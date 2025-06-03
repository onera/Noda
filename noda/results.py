# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Organize simulation results."""

import numpy as np

from noda import plots
import noda.composition_variables as cv
import noda.utils as ut


class SimulationResults:
    """
    Organize composition variables from a simulation.

    Attributes
    ----------
    saved_steps : list
        Steps for which simulation results are stored in resdict.
    saved_times : 1D array
        Times in h (rounded) that correspond to the saved steps.
    comps : list of str
        System constituents, with 'Va' first and dependent constituent
        last.
    inds : list of str
        Independent constituents, with 'Va' first.
    V_partial : dict
        Partial molar volumes.
    title : str
        Default plot title.
    results : dict
        Simulation results stored as :class:`UnitResult` instances.

    """

    def __init__(self, simu_parameters):
        """
        Class constructor.

        Parameters
        ----------
        resdict : dict
            Simulation results (dict).
        simu_parameters : dict
            Simulation parameters, see class attributes.

        """
        self.saved_times = simu_parameters['saved_times']
        self.comps = simu_parameters['comps']
        self.inds = self.comps[1:-1]
        self.V_partial = simu_parameters['V_partial']
        self.title = simu_parameters['title']
        self.results = {}
        self.saved_steps = None

    def add_results(self, resdict):
        """
        Populate results dictionary.

        Generate :class:`UnitResult` instances and set up their
        :class:`plots.StaticProfile` attribute.
        """
        self.saved_steps = list(resdict.keys())
        for n in self.saved_steps:
            res = UnitResult(resdict[n], self.V_partial)
            th = self.saved_times[self.saved_steps.index(n)]
            res.static_prof = plots.StaticProfile(res, n, th, self.title)
            self.results[n] = res
        self.results[-1] = self.results[self.saved_steps[-1]]

    def result(self, step_index=-1, th=None):
        """
        Access results at required step index or nearest to asked time.

        Parameters
        ----------
        step_index : int
            Index of required time step in saved_steps/saved_times.
        th : float
            Required time in hour.

        Returns
        -------
        :class:`UnitResult`
            Simulation results.

        Raises
        ------
        Exception
            If the user specifies both step_index and th.
            If step_index is not an integer or is out of range.

        """
        if self.saved_steps is None:
            msg = ("The results cannot be accessed because the simulation has "
                   "not been run or read.")
            raise ut.ResultsError(msg)
        if th is not None and step_index != -1:
            msg = "Cannot specify both 'th' and 'step'."
            raise ut.ResultsError(msg)
        if not isinstance(step_index, int):
            msg = "'step_index' must be an integer."
            raise ut.ResultsError(msg)
        if th is not None:
            step_index = np.argmin(abs(self.saved_times - th))
        try:
            step = self.saved_steps[step_index]
        except IndexError as exc:
            msg = ("'step_index' is out of range.\n"
                    "Must be lower than len of saved_steps "
                    f"({len(self.saved_steps)}).")
            raise ut.ResultsError(msg) from exc

        return self.results[step]

    def plot(self, varname='x', step_index=-1, th=None, **kwargs):
        """
        Plot profile of variable varname at required time.

        Select results at time step with :meth:`result` (see doc of
        :meth:`result` for handling of step_index and th). Then call
        :meth:`UnitResult.plot`.

        Parameters
        ----------
        varname : str
            Name of variable to be plotted.
        step_index : int
            Index of required time step in saved_steps/saved_times.
        th : float
            Required time in hour.

        Returns
        -------
        fig, ax : matplotlib figure and axis

        """
        res = self.result(step_index=step_index, th=th)
        fig, ax = res.plot(varname=varname, **kwargs)
        return fig, ax

    def plot_quartet(self, step_index=-1, th=None, **kwargs):
        """
        Plot 2x2 subplots with profiles of simulation results.

        Select results at time step with :meth:`result` (see doc of
        :meth:`result` for handling of step_index and th). Then call
        :meth:`UnitResult.plot_quartet`.

        Parameters
        ----------
        step_index : int
            Index of required time step in saved_steps/saved_times.
        th : float
            Required time in hour.

        Returns
        -------
        fig, axes, lines : matplotlib figure, axes and lines

        """
        res = self.result(step_index=step_index, th=th)
        fig, axes, lines = res.plot_quartet(**kwargs)
        return fig, axes, lines

    def interactive_plot(self, variable='x'):
        """Interactive plot (call :class:`plots.InteractivePlot`)."""
        if variable == 'quartet':
            ip = plots.InteractivePlotQuartet(
                            self.comps, self.results,
                            self.saved_times, self.saved_steps,
                            self.title)
        else:
            ip = plots.InteractivePlot(variable, self.comps, self.results,
                                       self.saved_times, self.saved_steps,
                                       self.title)
        return ip

    def add_mass_balance(self):
        """
        Add mass balance.

        For each time step in steps, integrate concentrations and volume along
        profile, and compute difference with initial values.

        Variables:

        * N : amount of vacancies and atom species
        * Np : amount of pores
        * NW : amount of vacancies annihilated due to volume variation since
          initial step.
        * Nv : amount of vacancies in simulation box, in solution and in
          pores, plus amount annihilated due to volume variations.

        Type of operation:

        * INx : integral amount in mol/m2
        * dINx : difference with initial value, mol/m2
        * rdINx : relative difference with initial value

        """
        r0 = self.results[0]
        L0 = r0.z[-1] - r0.z[0]

        for r in self.results.values():
            V0 = (self.V_partial['Va'] if self.V_partial['Va'] != 'local'
                  else r.Vm_mid)
            Vp = (self.V_partial['pore'] if self.V_partial['pore'] != 'local'
                  else r.Vm_mid)
            r.IN = np.sum(r.c_mid_arr*np.diff(r.z), axis=1)
            r.INp = np.sum(r.fp_mid/Vp*np.diff(r.z))
            L = r.z[-1] - r.z[0]
            r.INW = (L0 - L)/np.mean(V0)
            r.INv = r.IN[0] + r.INp + r.INW

            r.dIN = r.IN - r0.IN
            r.dINp = r.INp - r0.INp
            r.dINW = r.INW - r0.INW
            r.dINv = r.INv - r0.INv

            r.rdIN = r.dIN/r0.IN
            r.rdINv = r.dINv/r0.INv

    def plot_mass_balance(self):
        """Call :func:`plots.plot_mass_balance`."""
        fig, ax = plots.plot_mass_balance(self)
        return fig, ax


class UnitResult:
    """
    Gather composition variables from one time step of a simulation.
    
    Most attributes below are available in three versions:
    * x_mid: variable x evaluated on mid points
    * x_nod: variable x evaluated on nodes
    * x: shorthand for x_nod

    Attributes
    ----------
    z_nod : 1D array
        Node positions.
    z_mid : 1D array
        Midpoint positions.
    c : dict of 1D arrays
        System concentrations.
    V : 1D array
        Average molar volume of system.
    y : dict of 1D arrays
        Metal site fractions.
    yVa_eq : 1D array
        Equilibrium vacancy site fraction.
    dyVa : 1D array
        Difference between actual and equilibrium vacancy fraction.
    ryVa : 1D array
        Relative difference between actual and eq. vacancy fraction.
    x : dict of 1D arrays
        Metal atom fractions.
    Vm : 1D array
        Average molar volume of metal.
    mu : dict of 1D arrays
        Chemical potentials.
    Jlat : dict of 1D arrays
        Fluxes in lattice-fixed frame.
    v : 1D array
        Velocity field of lattice relative to laboratory frame.
    Jref : dict of 1D arrays
        Fluxes in laboratory frame.
    alpha_dislo : 1D array
        Lattice sink rate due to dislocations.
    alpha_pores : 1D array
        Lattice sink rate due to pore growth.
    alpha : 1D array
        Total lattice sink rate.
    L : 2D array
        Onsager coefficients.
    fm : 1D array
        Metal volume fraction.
    fp : 1D array
        Pore volume fraction.
    gamma_V : 1D array
        Relative volume variation rate due to molar volume variations.
    gamma_N : 1D array
        Relative volume variation rate due to variations in the number of
        lattice sites.
    gamma_p : 1D array
        Relative volume variation rate due to pore growth.
    gamma : 1D array
        Total relative volume variation.
    deformation : 1D array
        Relative length variation.
    static_prof : :class:`plots.StaticProfile`
        Configure and plot static profiles.
    """

    def __init__(self, result, V_partial):
        """
        Class constructor.

        Parameters
        ----------
        comps : list of str
            System constituents, list with 'Va' first and dependent constituent
            last.
        result : dict
            Simulation results at one time step.
        V_partial : dict
            Partial molar volumes.

        """
        comps = [k for k in V_partial if k != 'pore']
        assert comps[0] == 'Va'
        self.comps = comps
        self.z_nod = result['z']
        self.z_mid = (self.z_nod[1:] + self.z_nod[:-1])/2
        self.z = self.z_nod

        c_mid_arr = result['c']
        self.c_mid_arr = c_mid_arr  # Record for use in mass balance
        c_nod_arr = cv.arr_to_nod(c_mid_arr, self.z)
        self.c_mid = {k: c_mid_arr[i] for i, k in enumerate(comps)}
        self.c_nod = {k: c_nod_arr[i] for i, k in enumerate(comps)}
        self.V_mid = 1/sum(self.c_mid_arr)

        self.y_mid = {k: self.c_mid[k]*self.V_mid for k in comps}
        self.y_nod = {k: self.to_nod(v) for k, v in self.y_mid.items()}
        self.yVa_eq_mid = result['yVa_eq']
        self.dyVa_mid = self.y_mid['Va'] - self.yVa_eq_mid
        self.ryVa_mid = self.dyVa_mid / self.yVa_eq_mid

        self.x_mid = {k: self.y_mid[k]/(1 - self.y_mid['Va'])
                      for k in comps[1:]}
        self.x_nod = {k: self.to_nod(v) for k, v in self.x_mid.items()}
        self.Vm_mid = self.compute_Vm(V_partial, self.x_mid, self.y_mid)

        mu_arr = result['mu']
        self.mu_mid = {k: mu_arr[i] for i, k in enumerate(comps)}
        self.mu_nod = {k: self.to_nod(v) for k, v in self.mu_mid.items()}

        Jlat_atoms = result['Jlat']
        Jlat_Va = - sum(Jlat_atoms)
        Jlat_arr = np.vstack((Jlat_Va, Jlat_atoms))
        self.Jlat = {k: Jlat_arr[i] for i, k in enumerate(comps)}

        self.v = result['v']
        Jref_arr = Jlat_arr + self.v*c_nod_arr
        self.Jref = {k: Jref_arr[i] for i, k in enumerate(comps)}

        try:
            self.alpha_dislo_mid = result['alpha_d']
        except KeyError:
            self.alpha_dislo_mid = result['alpha_h']  # for back-compatibility

        self.alpha_pores_mid = result['alpha_p']
        self.alpha_mid = self.alpha_dislo_mid + self.alpha_pores_mid

        self.L_mid = result['L']
        self.L_nod = self.compute_L_nod(self.L_mid)

        self.fm_mid = result['fm']
        self.fp_mid = 1 - self.fm_mid

        self.gamma_V_mid = result.get('gamma_V', None)
        self.gamma_N_mid = self.alpha_mid*self.fm_mid
        self.gamma_p_mid = result.get('gamma_p', None)
        try:
            self.gamma_mid = (self.gamma_V_mid
                              + self.gamma_N_mid
                              + self.gamma_p_mid)
        except TypeError:
            self.gamma_mid = None  # for backward compatibility

        self.deformation_mid = result.get('deformation', None)
        self.make_nod_attributes()
        self.make_shorthands()

        self.static_prof = None

    def to_nod(self, vec):
        """Interpolate vector on nodes."""
        return cv.vec_to_nod(vec, self.z)

    def make_nod_attributes(self):
        """
        Make attributes on nodes.

        For every variable x in the list, create attribute x_nod by
        interpolating x_mid on nodes.
        """
        variables = ['V', 'yVa_eq', 'dyVa', 'ryVa', 'Vm', 'alpha_dislo',
                     'alpha_pores', 'alpha', 'fm', 'fp', 'gamma_N', 'gamma_V',
                     'gamma_p', 'gamma', 'deformation']
        for x in variables:
            if hasattr(self, x + '_mid'):
                on_mid = getattr(self, x + '_mid')
                if on_mid is not None:
                    on_nod = self.to_nod(on_mid)
                    setattr(self, x + '_nod', on_nod)

    def make_shorthands(self):
        """
        Make attribute shorthands.

        For every variable x in the list, create attribute with shorthand x and
        assign existing x_nod value.
        """
        variables = ['c', 'V', 'y', 'yVa_eq', 'dyVa', 'ryVa', 'x', 'Vm', 'mu',
                     'alpha_dislo', 'alpha_pores', 'alpha', 'L', 'fm', 'fp',
                     'gamma_N', 'gamma_V', 'gamma_p', 'gamma', 'deformation']
        for x in variables:
            if hasattr(self, x + '_nod'):
                setattr(self, x, getattr(self, x + '_nod'))

    def compute_Vm(self, V_partial, x_mid, y_mid):
        """Compute average molar volume of the metal phase."""
        if V_partial['Va'] == 'local':
            arr = np.array([V_partial[k] for k in self.comps[1:]])
            Vk_arr = np.hstack((np.nan, arr))
            Vk_arr = Vk_arr[np.newaxis].T
            x_arr = np.array(list(x_mid.values()))
            Vm = sum(x_arr*Vk_arr[1:])
        else:
            Vk_arr = np.array([V_partial[k] for k in self.comps])
            Vk_arr = Vk_arr[np.newaxis].T
            y_arr = np.array(list(y_mid.values()))
            Vm = sum(y_arr*Vk_arr)
        return Vm

    def compute_L_nod(self, L_mid):
        """Interpolate Onsager coefficients on nodes."""
        nind = L_mid.shape[0]
        L_nod = np.zeros((nind, nind, self.z_nod.size))
        for i, arr in enumerate(L_mid):
            L_nod[i] = cv.arr_to_nod(arr, self.z_nod)
        return L_nod

    def plot(self, varname='x', title=None, **kwargs):
        """
        Plot profile of variable varname.

        Call :meth:`plots.StaticProfile.single`. See list of variable names in
        class documentation.

        """
        fig, ax = self.static_prof.single(varname=varname, title=title, **kwargs)
        return fig, ax

    def plot_quartet(self, **kwargs):
        """
        Plot 2x2 grid plot with profiles of simulation results.

        Call :meth:`plots.StaticProfile.quartet`. Plot x, yVa, Jlat, fp.

        """
        fig, axes, lines = self.static_prof.quartet(**kwargs)
        return fig, axes, lines
