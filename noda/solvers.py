# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Solve the diffusion equation."""

import time

import numpy as np
import scipy.interpolate as spi

import noda.lattice as la
import noda.utils as ut
from noda.utils import div

# =============================================================================
# Solver


def solver(thermo, mobility, space, init, BC, time_grid, lattice,
           show_completion, verbose, L_mean_kind, logger):
    """
    Solve diffusion equation.

    Explicit (forward Euler) solver with non-equilibrium vacancies and pores,
    finite sink strengths that represent dislocation climb and pore growth.

    Parameters
    ----------
    thermo : :class:`thermodynamics.Thermodynamics`
        Thermodynamic properties.
    mobility : :class:`mobility.Mobility`
        Mobility properties.
    space : :class:`space.SpaceGrid`
        Space grid.
    init : :class:`initial_conditions.InitialConditions`
        Initial conditions.
    BC : dict of :class:`boundary_conditions.BoundaryConditions`
        Boundary conditions.
    time_grid : :class:`time.TimeGrid`
        Time parameters.
    lattice : :class:`lattice.Lattice`
        Parameters related to sink strength.
    show_completion : bool
        Print completion rate while solver is running.
    verbose : str
        Verbosity level.
    L_mean_kind : str
        Kind of mean used to compute L values at nodes. See
        :func:`make_L_mean_fun`.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    :class:`utils.AtomFractionError`
        If any site fraction goes below 0 or above 1.

    Returns
    -------
    res: dict of dicts
        Dict of stored variables:

        ``{n: {var: array for var in variables} for n in saved_steps}``

        where variables are:

        * z : node positions, 1D array
        * c : concentrations, 2D array
        * fp : pore fractions, 1D array
        * mu : chemical potentials, 2D array
        * Jlat : fluxes in lattice-fixed frame, 2D array
        * yVa_eq : equilibrium vacancy fraction, 1D array
        * v : velocity field of lattice in laboratory frame, 1D array.
        * gamma_V : relative volume variation rate associated with molar volume
          variation, 1D array.
        * gamma_p : relative volume variation rate associated with pore growth,
          1D array.
        * alpha_d : sink term associated with dislocation climb, 1D array.
        * alpha_p : sink term associated with pore growth, 1D array.
        * L : Onsager coefficients, 3D array
        * deformation : cumulative deformation in diffusion direction (eps_zz),
          1D array.

        The first 3 variables (z, c, fp) are stored at the start of the time
        step, while the others are stored at the end of the step.
    dur : float
        Solver run time in s.

    """
    start_time = time.time()
    cvar = init.cvar
    comps = cvar.comps
    V_partial = cvar.V_partial
    Vk = cvar.Vk_arr
    z = space.z_init.copy()
    dz = np.diff(z)
    dz_min = dz.min()
    dz_max = dz.max()
    deformation = np.zeros(dz.size)
    nt = time_grid.nt
    dt = time_grid.dt
    k_dislo = lattice.k_dislo
    k_pores = lattice.k_pores
    rho_dislo = lattice.rho_dislo
    rho_pores = lattice.rho_pores
    L_mean_fun = make_L_mean_fun(L_mean_kind)
    res = {0: {'z': z.copy(), 'c': cvar.c.mid, 'fp': cvar.fp.mid}}
    completion_rates = np.arange(0.1, 1, step=0.1)
    completion_steps = [np.round((nt - 1)*r) for r in completion_rates]
    if show_completion:
        logger.info(f"{'Run':^10}{'':2}{'Completion':^10}")
        logger.info(f"{'time (s)':^10}{'':2}{'rate (%)':^10}")

    # The variables calculated in the loop (MU, Jlat, ...) are based on the
    # composition at the beginning of the loop, but are stored at the end of
    # the loop. Therefore the loop is run once more than needed for the
    # composition, in order to record the other variables at the last step.
    for n in range(nt):

        y = cvar.y.mid
        x = cvar.x.mid
        V = cvar.V.mid
        Vm = cvar.Vm.mid
        c = cvar.c.mid
        fm = cvar.fm.mid

        V0 = V_partial['Va'] if V_partial['Va'] != 'local' else Vm
        Vp = V_partial['pore'] if V_partial['pore'] != 'local' else Vm

        if n in time_grid.saved_steps:
            res[n] = {}
            res[n]['z'] = z.copy()
            res[n]['c'] = c.copy()
            res[n]['fm'] = fm.copy()
            res[n]['deformation'] = deformation.copy()

        if (y < 0).any() or (y > 1).any() or (y == np.nan).any():
            raise ut.AtomFractionError(n)

        # Compute chemical potentials
        MU = thermo.MU_funy(y[:-1])
        MU_diff = MU[1:] - MU[0]
        yVa_eq = thermo.yVa_fun(x[:-1])

        # Compute Onsager coefficients
        # Note: since c is the global concentration, L_eq includes fm. This is
        # only valid inasmuch as L_ij = 0 in the pores
        L_eq = mobility.L_fun(c[1:], x[:-1])
        # Factor in actual yVa
        L = L_eq * y[0]/yVa_eq
        # Select diagonal terms
        L_diag = np.array([L[k, k] for k in range(len(comps[1:]))])
        L_nod = L_mean_fun(L_diag)
        dz_nod = (dz[1:] + dz[:-1])/2
        Jbulk = -L_nod*np.diff(MU_diff)/dz_nod
        Jleft, Jright = compute_boundary_fluxes(n*dt, MU_diff, dz, BC,
                                                thermo.MU_funy,
                                                L, mobility.L_fun, Vk)
        Jlat = np.hstack((Jleft[None].T, Jbulk, Jright[None].T))

        if lattice.ideal:
            alphas = compute_alpha_ideal(z, cvar.c.nod(z), c, y[0], dt, Jlat,
                                         Vk, V0, Vp, V, Vm, thermo.yVa_fun,
                                         space.geometry)
        else:
            if rho_dislo or rho_pores:
                D0 = mobility.DVa_fun(y, x[:-1])
                if rho_dislo:
                    k_dislo = rho_dislo * D0
                if rho_pores:
                    k_pores = rho_pores * D0
            alphas = la.compute_alpha_nonideal(dt, y[0], yVa_eq, V, Vp, fm,
                                               k_dislo, k_pores)
        alpha_d, alpha_p = alphas

        c, v, gamma_V, gamma = solver_core(z, cvar.c.nod(z), c, y[0], dt, Jlat,
                                           Vk, V0, Vp, V, Vm, alpha_d, alpha_p,
                                           space.geometry)
        deformation += gamma*dt

        # New boundary positions
        z[0] = z[0] + v[0]*dt
        z[-1] = z[-1] + v[-1]*dt
        dz = np.diff(z)

        # Add/remove nods as needed
        if any(dz < dz_min/2) or any(dz > dz_max*3/2):
            z, c, deformation = remesh(n, z, dz_min, dz_max, c, deformation,
                                       verbose, logger)
            dz = np.diff(z)

        cvar.c.mid = c

        # Record variables
        if n in time_grid.saved_steps:
            res[n]['mu'] = MU
            res[n]['Jlat'] = Jlat
            res[n]['yVa_eq'] = yVa_eq
            res[n]['v'] = v
            res[n]['gamma_V'] = gamma_V
            res[n]['gamma_p'] = -alpha_p*Vp/V
            res[n]['alpha_d'] = alpha_d
            res[n]['alpha_p'] = alpha_p
            res[n]['L'] = L

        if show_completion:
            if n in completion_steps:
                dur = time.time() - start_time
                r = n/(nt - 1)
                logger.info(f'{dur:^10.3f}{"":2}{r*100:^10.0f}')

    dur = time.time() - start_time
    logger.info(f'Solver run time: {dur:.3f} s')

    return res


def solver_core(z, cnod, c, y_Va, dt, Jlat, Vk, V0, Vp, V, Vm,
                alpha_d, alpha_p, geometry):
    """
    Compute fluxes in laboratory frame and then new concentrations.

    See argument definitions in :func:`solver`.

    """
    alpha = alpha_d + alpha_p
    gamma_V = la.compute_gamma_V(Jlat, z, Vk, V0, Vm, V, y_Va, alpha, geometry)
    gamma_N = Vm/V*alpha
    gamma_p = -alpha_p*Vp/V
    gamma = gamma_V + gamma_N + gamma_p
    v = la.compute_velocity(gamma, z, Vk, V0, Jlat, y_Va, geometry)

    cnew = np.zeros(c.shape)

    Jref = Jlat + cnod[1:]*v
    derc = - div(Jref, z, geometry)
    cnew[1:] = c[1:] + dt*derc

    Jlat_Va = -sum(Jlat)
    Jref_Va = Jlat_Va + cnod[0]*v
    derc = - div(Jref_Va, z, geometry) + alpha/V
    cnew[0] = c[0] + dt*derc

    return cnew, v, gamma_V, gamma


def compute_alpha_ideal(z, cnod, c, y_Va, dt, Jlat, Vk, V0, Vp, V, Vm,
                        yVa_fun, geometry):
    """
    Compute ideal sink terms.

    See argument definitions in :func:`solver`.

    The term related to pore growth, alpha_p, is zero.

    The term related to dislocation climb, alpha_d, is computed from analytical
    expression if the equilibrium vacancy fraction is composition-fixed. If
    not, an iterative method is used instead, with the analytical version as a
    starting point. See Gheno 2022 [#Gheno_2022]_ for details.

    The error increases when partial molar volumes are different. The error is
    reduced by increasing the number of loops. Two loops is found to be
    sufficient in practice (y0 - y0_eq = 1e-13 - 1e-14).
    """
    Jlat_Va = -sum(Jlat)
    # Analytical expression
    alpha_d = V/(1 - y_Va)*div(Jlat_Va, z, geometry)
    alpha_p = np.zeros(alpha_d.size)

    # TODO : this can be made optional if yVa_eq is constant, ie when all L_kVa
    # interaction parameters are the same.
    for _ in range(2):
        c_, _, _, _ = solver_core(z, cnod, c, y_Va, dt, Jlat, Vk, V0,
                                  Vp, V, Vm, alpha_d, alpha_p, geometry)
        y_ = c_/sum(c_)
        x_ = y_[1:]/(1 - y_[0])
        yVa_eq_ = yVa_fun(x_[:-1])
        alpha_d += (yVa_eq_ - y_[0])/dt/(1 - y_[0])

    return alpha_d, alpha_p


# =============================================================================
# Auxiliary functions


def compute_boundary_fluxes(t, MU_diff, dz, BC, MU_fun, L, L_fun, Vk):
    """
    Compute fluxes on domain boundaries.

    Parameters
    ----------
    t : float
        Time (s).
    MU_diff : 2D array
        Diffusion potentials (MU_k - MU_0), initial shape (`ninds` + 1,
        `nz` - 1).
    dz : 1D array
        Space step (m).
    BC : dict
        Types and functions of left and right BC, see
        :class:`boundary_conditions.BoundaryConditions`.
    MU_fun : function
        Compute chemical potentials from site fractions.
    L : 3D array
        Onsager coefficients, initial shape (`ninds` + 1, `ninds` + 1,
        `nz` - 1).
    Vk : 1D array
        Partial molar volumes, shape (`ninds` + 2).

    Returns
    -------
    J['left'] : 1D array
        Flux on left-hand boundary, shape (`ninds` + 1,).
    J['right'] : 1D array
        Flux on right-hand boundary, shape (`ninds` + 1,).

    """
    J = {}

    for side in ['left', 'right']:
        sign = 1 if side == 'left' else -1
        i0 = 0 if side == 'left' else -1

        if BC[f'{side}'].type == 'Dirichlet':
            cvar_BC = BC[side].cvar_fun(t)
            y_BC = cvar_BC.y.mid
            MU_BC = MU_fun(y_BC[:-1])
            MU_diff_BC = MU_BC[1:] - MU_BC[0]
            gen_comps = range(MU_diff_BC.shape[0])
            L_BC_fullarray = L_fun(cvar_BC.c.mid[1:], cvar_BC.x.mid[:-1])
            L_BC = np.array([L_BC_fullarray[k, k, 0] for k in gen_comps])
            L_mid = np.array([L[k, k, i0] for k in gen_comps])
            L_side = (L_BC + L_mid)/2
            RL_side = 0.5*dz[i0]/L_side
            J[f'{side}'] = (-sign*(MU_diff[:, i0] - MU_diff_BC.T[0])
                            / RL_side)
            # Compensate volume change by setting flux of dependent constituent
            J[f'{side}'][-1] = (-sum(J[f'{side}'][:-1]*Vk[1:-1])/Vk[-1]).item()
        else:
            J[f'{side}'] = sign*BC[side].J_fun(t)

    return J['left'], J['right']


def make_L_mean_fun(L_mean_kind):
    r"""
    Make function that computes mean L values.

    L is defined in a volume. Fluxes are defined between volumes (at node
    points), and therefore require L evaluated at node points. The manner in
    which two neighboring L values are averaged to provide L at node points is
    determined by the 'L_mean_kind' parameter.

    Parameters
    ----------
    L_mean_kind : str
        Manner in which two neighboring Onsager coefficients are averaged.
        Possible values:

        * ``arithmetic``: :math:`\bar{L}_i = (L_i + L_{i - 1})/2`.
        * ``harmonic``: :math:`\bar{L}_i = \frac{2}{1/L_i + 1/L_{i - 1}}`.
        * ``geometric``: :math:`\bar{L}_i = \sqrt{L_iL_{i - 1}}`.
    
    Raises
    ------
    utils.UserInputError
        If L_mean_kind is invalid.

    Returns
    -------
    fun :
        Function that computes mean L at node points.

    """
    if L_mean_kind == 'arithmetic':
        fun = lambda L: (L[:, 1:] + L[:, :-1])/2
    elif L_mean_kind == 'harmonic':
        fun = lambda L: 2/(1/L[:, 1:] + 1/L[:, :-1])
    elif L_mean_kind == 'geometric':
        fun = lambda L: np.sqrt(L[:, 1:] * L[:, :-1])
    else:
        msg = f'L_mean_kind "{L_mean_kind}" not implemented'
        raise ut.UserInputError(msg) from None
    return fun

# =============================================================================
# Node adding/deleting algorithm


def remesh(n, z, dz_min, dz_max, cmid, deformation, verbose, logger):
    """
    Add or delete nodes.

    Parameters
    ----------
    n : int
        Time step.
    z : 1D array
        Node positions, initial shape (`nz`,).
    dz_min : float
        Minimum of initial dz.
    dz_max : float
        Maximum of initial dz.
    cmid : 2D array
        Concentrations, initial shape (`ninds` + 1, `nz` - 1).
    deformation : 1D array
        Cumulated deformation, initial shape (`nz` - 1,).
    verbose : str
        Verbosity level.

    Returns
    -------
    new_z : 1D array
        New node positions.
    new_cmid : 2D array
        New concentrations.
    new_deformation : 1D array
        New deformation.

    """
    dz = np.diff(z)
    zm = (z[1:] + z[:-1])/2

    delete_condition = dz < dz_min/2
    split_condition = dz > dz_max*3/2

    if any(delete_condition):
        idx_to_delete = np.nonzero(delete_condition)[0] + 1

        # avoid deleting last node
        if idx_to_delete[-1] == z.size - 1:
            idx_to_delete[-1] = z.size - 2

        if verbose > 0:
            text = f'Step {n}: deleting node(s) at index(ices) {idx_to_delete}'
            logger.info(text)

        idx_to_keep = [i for i in range(z.size) if i not in idx_to_delete]
        new_z = z[idx_to_keep]

    else:
        idx_to_split = np.nonzero(split_condition)[0]

        if verbose > 0:
            text = f'Step {n}: adding node(s) at index(ices) {idx_to_split}'
            logger.info(text)

        new_z = np.insert(z, idx_to_split + 1, zm[idx_to_split])

    new_zm = (new_z[1:] + new_z[:-1])/2
    new_cmid = np.zeros((len(cmid), new_z.size - 1))
    for i, cm in enumerate(cmid):
        f = spi.interp1d(zm, cm, fill_value='extrapolate')
        new_cmid[i] = f(new_zm)
    f_def = spi.interp1d(zm, deformation, fill_value='extrapolate')
    new_deformation = f_def(new_zm)

    return new_z, new_cmid, new_deformation
