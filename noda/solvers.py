# Copyright 2025 Onera
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


def solver(z_init, geometry, cvar, MU_fun, L_fun, yVa_fun, nt, dt,
           ideal_lattice, k_dislo, k_pores, rho_dislo, rho_pores, DVa_fun,
           Va_method, saved_steps, BC, show_completion, verbose, stencil,
           logger):
    """
    Solve diffusion equation.

    Explicit (forward Euler) solver with non-equilibrium vacancies and pores,
    finite sink strengths that represent dislocation climb and pore growth.

    Parameters
    ----------
    z_init : 1D array
        Initial node positions, shape (`nz`,).
    geometry : str
        Domain geometry (planar, cylindrical or spherical).
    cvar : :class:`composition_variables.CompositionVariables`
        Composition variables (x, y, c, Vm, ...).
    MU_fun : function
        Compute chemical potential from site fractions (see `MU_funy` in
        :class:`simu.System`).
    L_fun : func
        Same for phenomenological coefficients.
    yVa_fun : function
        Compute equilibrium vacancy site fraction, see
        :meth:`simu.System.add_Va_model`.
    nt : int
        Number of time steps, including time 0.
    dt : float
        Time step (s).
    ideal_lattice : bool
        Whether lattice is ideal, i.e., vacancy fraction is kept at
        equilibrium.
    k_dislo : float or str
        Sink strength associated with dislocation climb (str should be a valid
        file name in current job folder)
    k_pores : float or str
        Sink strength associated with pore growth (str should be a valid file
        name in current job folder)
    rho_dislo : float or str
        Dislocation density, used to compute k_dislo from local DVa.
    rho_pores : float
        Density used to compute k_pores from local DVa.
    DVa_fun : function
        Compute vacancy diffusion coefficient from atom fractions and site
        fractions, see :func:`thermo_functions.make_DVa_fun`.
    Va_method : str
        Type of function for equilibrium vacancy site fraction. Determined
        based on `thermo_method` in :meth:`simu.System.add_thermo_functions`.
    saved_steps : list
        Steps for which results are returned.
    BC : dict
        Types and functions of left and right BC, see
        :meth:`simu.Simulation.add_BC`.
    show_completion : bool
        Print completion rate while solver is running.
    verbose : str
        Verbosity level.
    stencil : str
        Name of discretization stencil. See :func:`compute_resistance`.

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

    comps = cvar.comps
    V_partial = cvar.V_partial
    Vk = cvar.Vk_arr

    z = z_init.copy()
    dz = np.diff(z)
    dz_min = dz.min()
    dz_max = dz.max()
    deformation = np.zeros(dz.size)

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

        if n in saved_steps:
            res[n] = {}
            res[n]['z'] = z.copy()
            res[n]['c'] = c.copy()
            res[n]['fm'] = fm.copy()

        if (y < 0).any() or (y > 1).any() or (y == np.nan).any():
            raise ut.AtomFractionError(n)

        # Compute chemical potentials
        MU = MU_fun(y[:-1])
        MU_diff = MU[1:] - MU[0]
        yVa_eq = yVa_fun(x[:-1])

        # Compute Onsager coeffs, factor in yVa
        # Note: since c is the global concentration, L_eq includes fm. This is
        # only valid inasmuch as L_ij = 0 in the pores
        L_eq = L_fun(c[1:], x[:-1])
        L = L_eq * y[0]/yVa_eq

        # Compute diffusion fluxes
        RL = compute_resistance(comps[1:], L, dz, stencil)

        Jbulk = -np.diff(MU_diff)/RL
        Jleft, Jright = compute_boundary_fluxes(n*dt, MU_diff, dz, BC, MU_fun,
                                                L, L_fun, Vk)
        Jlat = np.hstack((Jleft[None].T, Jbulk, Jright[None].T))

        if ideal_lattice:
            alphas = compute_alpha_ideal(z, cvar.c.nod(z), c, y[0], dt, Jlat,
                                         Vk, V0, Vp, V, Vm, yVa_fun, Va_method,
                                         geometry)
        else:
            if rho_dislo is not None or rho_pores is not None:
                D0 = DVa_fun(y, x[:-1])
                if rho_dislo is not None:
                    k_dislo = rho_dislo * D0
                if rho_pores is not None:
                    k_pores = rho_pores * D0
            alphas = la.compute_alpha_nonideal(dt, y[0], yVa_eq, V, Vp, fm,
                                               k_dislo, k_pores)
        alpha_d, alpha_p = alphas

        c, v, gamma_V, gamma = solver_core(z, cvar.c.nod(z), c, y[0], dt, Jlat,
                                           Vk, V0, Vp, V, Vm, alpha_d, alpha_p,
                                           geometry)
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
        if n in saved_steps:
            res[n]['mu'] = MU
            res[n]['Jlat'] = Jlat
            res[n]['yVa_eq'] = yVa_eq
            res[n]['v'] = v
            res[n]['gamma_V'] = gamma_V
            res[n]['gamma_p'] = -alpha_p*Vp/V
            res[n]['alpha_d'] = alpha_d
            res[n]['alpha_p'] = alpha_p
            res[n]['L'] = L
            res[n]['deformation'] = deformation

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
                        yVa_fun, Va_method, geometry):
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

    if Va_method != 'cst':
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
        :meth:`simu.Simulation.add_BC`.
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

        if BC[f'{side}'] == 'Dirichlet':
            cvar_BC = BC[f'c_{side}'](t)
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
            J[f'{side}'] = sign*BC[f'J_{side}'](t)

    return J['left'], J['right']


def compute_resistance(comps, L, dz, stencil):
    """
    Compute diffusion resistance.

    Defined as dz/L, with L the Onsager coefficients.

    Parameters
    ----------
    comps : list of str
        System constituents.
    L : 3D array
        Onsager coefficients, initial shape (`ninds` + 1, `ninds` + 1,
        `nz` - 1).
    dz : 1D array
        Space step, initial shape (`nz` - 1,).
    stencil : str
        Name of discretization stencil. Possible values: 'H', 'A', 'G'.

    Raises
    ------
    Exception
        If stencil is not supported.

    Returns
    -------
    RL : 2D array
        Diffusion resistance, initial shape (`ninds` + 1, `nz` - 2).

    """
    gen_comps = range(len(comps))
    RL = None
    if stencil == 'H':
        RL_mid = np.array([dz/L[k, k] for k in gen_comps])
        RL = (RL_mid[:, 1:] + RL_mid[:, :-1])/2

    elif stencil == 'A':
        L_mid = np.array([L[k, k] for k in gen_comps])
        L_bar = (L_mid[:, 1:] + L_mid[:, :-1])/2
        dz_bar = (dz[1:] + dz[:-1])/2
        RL = dz_bar/L_bar

    elif stencil == 'G':
        L_mid = np.array([L[k, k] for k in gen_comps])
        L_bar = np.sqrt(L_mid[:, 1:] * L_mid[:, :-1])
        dz_bar = (dz[1:] + dz[:-1])/2
        RL = dz_bar/L_bar
    return RL

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
