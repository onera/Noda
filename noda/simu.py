# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define and load diffusion simulation."""

from pathlib import Path
import datetime

import numpy as np

from noda.alloy_system import AlloySystem, intro_msg
from noda.log_utils import get_and_log
import noda.para_io as pa
import noda.composition_variables as cv
import noda.meshing as mesh
import noda.solvers as so
import noda.results_io as rs
import noda.utils as ut
from noda import results


class Simulation(AlloySystem):
    """
    Thermodynamic system, simulation conditions, simulation results.

    This is a base class that does some pre-processing and provides methods to
    set up simulation conditions, run simulation and plot the results.

    It is possible to create a new simulation directly with this class, in
    which case the setup methods have to be called by the user. This provides
    the possibility of changing the system properties and simulation conditions
    in a user script, with functions not provided by this class (or the parent
    :class:`alloy_system.AlloySystem`). Note that this should be done with
    great care, as the custom modifications may cause unexpected issues and
    will not be recorded in the log file.

    The recommended way to create a simulation is to use :class:`NewSimulation`
    (new simulation) or :class:`ReadSimulation` (read from log file), where the
    default setup methods will be called.

    Also note that using :class:`ReadSimulation` to load a simulation first
    created with :class:`Simulation` and custom setup functions (instead of
    :class:`NewSimulation`) may produce an inconsistent state.

    This is a 1D version, with time-invariant space grid and time step.

    Attributes
    ----------
    x_profile : str
        Type of initial atom fraction profiles. See accepted values in
        :func:`make_profile`.
    zstep : float
        Position of step along step profile in m (see doc in
        :func:`meshing.make_step_profile`).
    q : float
        Common ratio for geometric grid.
    x_left : dict of floats
        Initial atom fractions of independent constituents on left-hand side.
    x_right : dict of floats
        Initial atom fractions of independent constituents on right-hand side.
    yVa_profile : str
        Type of initial vacancy site fraction profile. See accepted values in
        :func:`make_profile`.
    yVa_left : float
        Initial vacancy site fraction on left-hand side.
    yVa_right : float
        Initial vacancy site fraction on right-hand side.
    fp_profile : str
        Type of initial pore volume fraction profile. See accepted values in
        :func:`make_profile`.
    fp_left : float
        Initial pore volume fraction on left-hand side.
    fp_right : float
        Initial pore volume fraction on right-hand side.
    stencil : str
        Discretization stencil. See accepted values in
        :func:`solvers.compute_resistance`.
    title : str
        Title to be used in plots. See :func:`plots.make_plot_title`.
    zmin : float
        Initial position of left-hand domain boundary (m).
    zmax : float
        Initial position of right-hand domain boundary (m).
    nz_init : int
        Number of steps in initial space grid.
    z_init : 1D array
        Initial node positions. Ranges from 0 to zmax, size nz. This is where
        fluxes are evaluated.
    dz_init : 1D array
        Initial space step, shape (nz - 1,).
    zm_init : 1D array
        Initial midpoint positions, shape (nz - 1,). This is where
        composition variables are evaluated.
    geometry : str
        Domain geometry (planar, cylindrical or spherical).
    th : float
        Simulation time in h.
    ts : float
        Simulation time in s.
    dt_multiplier : float
        Factor by which default time step is multiplied (see :meth:`add_time`).
    dt : float
        Time step size in s.
    nt : int
        Number of time steps (positions in time sequence), including 0.
    num_out : int
        Number of time steps where simulation results are stored and saved on
        file.
    saved_steps : list
        Steps for which simulation results are stored and saved on file.
    saved_times : list
        Times in h (rounded) for which simulation results are stored and saved
        on file.
    x_init : dict of 1D arrays
        Initial atom fraction profiles of independent constituents.
    yVa_init : 1D array
        Initial vacancy site fraction profile.
    fm_init : 1D array
        Initial pore volume fraction profile.
    cvar : :class:`composition_variables.CompositionVariables`
        Composition variables (x, y, c, Vm, ...).
    BC : dict
        Type of boundary condition (Dirichlet, Neumann) and corresponding
        values or expressions. See  :meth:`add_BC`.
    res : dict
        Simulation results, nested dicts with syntax res[th][var][k].
    simres : :class:`results.SimulationResults`
        Simulation results.
    """

    def __init__(self, ref):
        """
        Class constructor.

        Use :class:`alloy_system.AlloySystem` constructor to define system
        parameters, then pre-process some simulation conditions. See
        :func:`para_io.read_parameter_file` for how input file is parsed.

        Perform consistency checks: if lattice is ideal, custom initial vacancy
        and pore fraction profiles are not allowed (vacancy fraction profile
        will be set to equilibrium, pore fraction profile will be set to 0).

        """
        super().__init__(ref)
        self.logger.info("Creating '%s' simulation.", ref)

        self.x_profile = self.params['x_profile']
        self.x_left = self.params.get('x_left', None)
        self.x_right = self.params.get('x_right', None)

        self.yVa_profile = self.params.get('yVa_profile', None)
        self.yVa_left = self.params.get('yVa_left', None)
        self.yVa_right = self.params.get('yVa_right', None)
        if self.yVa_profile is not None and self.ideal_lattice is True:
            msg = "Custom initial yVa profile not compatible with ideal "
            msg += "lattice"
            raise ut.UserInputError(msg) from None

        self.fp_profile = self.params.get('fp_profile', None)
        self.fp_left = self.params.get('fp_left', None)
        self.fp_right = self.params.get('fp_right', None)
        if self.fp_profile is not None and self.ideal_lattice is True:
            msg = "Custom initial fp profile not compatible with ideal "
            msg += "lattice"
            raise ut.UserInputError(msg) from None

        self.stencil = self.params.get('stencil',
                                       self.default_parameters['stencil'])
        if self.stencil not in ['H', 'A', 'G']:
            msg = f'Stencil "{self.stencil}" not implemented'
            raise ut.UserInputError(msg) from None
        self.title = make_plot_title(self.dep, self.x_profile, self.x_left,
                                     self.x_right, self.TC)
        self.logger.data(f"title = {self.title}")
        self.check_boundary_conditions()
        self.geometry = get_and_log(self.params, 'geometry',
                                    self.default_parameters['geometry'],
                                    self.logger)
        if self.geometry not in ('planar', 'cylindrical', 'spherical'):
            msg = f"Unknown geometry parameter {self.geometry}."
            raise ut.UserInputError(msg) from None
        self.check_grid_config()

        # Initialize other attributes
        self.zmin = None
        self.zmax = None
        self.nz_init = None
        self.z_init = None
        self.dz_init = None
        self.zm_init = None
        self.q = None
        self.zstep = None
        self.th = None
        self.ts = None
        self.dt_multiplier = None
        self.dt = None
        self.nt = None
        self.num_out = get_and_log(self.params, 'num_out',
                                   self.default_parameters['num_out'],
                                   self.logger)
        self.saved_steps = None
        self.saved_times = None
        self.x_init = None
        self.yVa_init = None
        self.fm_init = None
        self.cvar = None
        self.BC = None
        self.simres = None

    def check_boundary_conditions(self):
        """
        Check that boundary conditions are consistent.

        Make sure all expected constituents are present in the declared
        boundary conditions.

        Call :meth:`unit_check_BC_constituents` for the 4 conditions.

        * Condition on x: expect to find all independent constituents.
        * Condition on J: expect to find all atom constituents

        """
        for var in ['xBC_left', 'xBC_right']:
            if var in self.params:
                self.unit_check_BC_constituents(var, set(self.inds))
                self.unit_check_BC_values(var)
        for var in ['JBC_left', 'JBC_right']:
            if var in self.params:
                self.unit_check_BC_constituents(var, set(self.comps[1:]))

    def unit_check_BC_constituents(self, var, expected):
        """Make sure all expected constituents are present in the BC."""
        try:
            found = set(list(self.params[var]))
            assert found == expected
        except AssertionError as exc:
            name = var[0] + var[3:]
            msg = (f"Missing or extra elements in '{name}' boundary "
                   "condition.\n"
                   f"expected {expected}\n"
                   f"found {found}")
            raise ut.UserInputError(msg) from exc

    def unit_check_BC_values(self, var):
        """Make sure atom fractions are within accepted bounds."""
        min_val = self.default_parameters['min_atom_fraction']
        max_val = 1 - self.default_parameters['min_atom_fraction']
        for k in self.params[var]:
            v = float(self.params[var][k])
            if v < min_val:
                msg = (f"Input {var} {k} = {v} replaced "
                       f"by minimum allowed {min_val}.")
                self.logger.info(msg)
                v = min_val
                self.params[var][k] = str(min_val)
            if v > max_val:
                msg = (f"Input {var} {k} = {v} replaced "
                       f"by maximum allowed {max_val}.")
                self.logger.info(msg)
                v = max_val
                self.params[var][k] = str(max_val)

    def check_grid_config(self):
        """Make sure grid-related input is consistent."""
        if 'grid' not in self.params:
            if 'zmax' not in self.params:
                msg = "Grid not specified. Parameter zmax is required."
                raise ut.UserInputError(msg) from None
        elif not ut.isfilename(self.params['grid']):
            if 'zmax' not in self.params:
                msg = ("Grid parameter is not a file name "
                       f"(grid = '{self.params['grid']}'). "
                       "Parameter zmax is required.")
                raise ut.UserInputError(msg) from None
        else:
            for x in ['nz', 'zmin', 'zmax']:
                if x in self.params:
                    msg = ("Found filename as grid parameter "
                           f"(grid = '{self.params['grid']}'). "
                           f"Cannot specify {x}.")
                    raise ut.UserInputError(msg) from None

    def add_grid(self):
        """
        Add initial space grid.

        Behavior controlled by the 'grid' parameter:

        * 'linear': linear grid from zmin to zmax with size nz (zmax and nz
          must be included in parameter input file, zmin defaults to 0).
        * 'geo': geometric grid from zmin to zmax with size nz and common ratio
          q (zmax and nz must be included in input file, zmin and q are
          optional).
        * filename: read from sdir/filename using np.genfromtxt (zmin, zmax and
          nz are inferred from the grid).

        """
        grid = get_and_log(self.params, 'grid', 'linear', self.logger)

        if grid in ('linear', 'geo'):
            self.zmin = self.params.get('zmin', 0)
            self.zmax = self.params['zmax']
            self.nz_init = get_and_log(self.params, 'nz',
                               self.default_parameters['number_space_steps'],
                               self.logger)
            if grid == 'linear':
                self.z_init = np.linspace(self.zmin, self.zmax,
                                          num=self.nz_init)
            elif grid == 'geo':
                self.q = get_and_log(self.params, 'q',
                               self.default_parameters['common_ratio'],
                               self.logger)
                self.z_init = mesh.geo_grid(self.zmin, self.zmax, self.nz_init,
                                            self.q)
        else:
            fpath = self.work_dir / grid
            self.z_init = np.genfromtxt(fpath)
            self.zmin = self.z_init[0]
            self.zmax = self.z_init[-1]
            self.nz_init = self.z_init.size

        self.dz_init = np.diff(self.z_init)
        self.zm_init = (self.z_init[:-1] + self.z_init[1:])/2

        if self.geometry in ['cylindrical', 'spherical']:
            if self.zmin < 0:
                msg = (f"Found strictly negative zmin in {self.geometry} "
                       f"geometry (zmin = '{self.zmin}'). "
                       "zmin should be positive.")
                raise ut.UserInputError(msg) from None
            if self.zmax <= 0:
                msg = (f"Found negative zmax in {self.geometry} geometry "
                       f"(zmax = '{self.zmax}'). "
                       "zmax should be strictly positive.")
                raise ut.UserInputError(msg) from None
            if self.zmin >= self.zmax:
                msg = (f"Found zmin >= zmax in {self.geometry} geometry "
                       f"(zmin = '{self.zmin}', zmax = '{self.zmax}'). "
                       "zmin should be strictly smaller than zmax.")
                raise ut.UserInputError(msg) from None

    def add_time(self):
        """
        Add time step and number of time steps.

        The time step dt is calculated to ensure the stability of explicit
        schemes, with a default Fourier number Fo = 0.4. The DT value is 
        calculated from the initial concentration profile. The time step
        can be made smaller using the user-specified dt_multiplier (or
        nt_multiplier, see :func:`para_io.read_parameter_string`).

        """
        self.th = self.params['th']
        self.ts = self.th*3600
        self.dt_multiplier = self.params.get('dt_multiplier', 1)

        x_arr = np.array(list(self.x_init.values()))
        DTmax = np.max(self.DT_fun(x_arr))
        Fo = self.default_parameters['Fourier_number']
        self.dt = Fo*self.dt_multiplier*self.dz_init.min()**2/DTmax
        self.nt = int(self.ts/self.dt) + 1
        # Adjust dt to recover ts (rounding nt introduces error)
        self.dt = self.ts/(self.nt - 1)

        # make sure there are at least 10 time steps
        if self.nt < self.default_parameters['min_number_time_steps']:
            self.nt = self.default_parameters['min_number_time_steps']
            self.dt = self.ts/(self.nt - 1)

    def add_saved_steps(self):
        """
        Add list of time steps where simulation results will be saved to file.

        The time steps are determined based on the num_out parameter.
        """
        if self.num_out == 'all':
            self.num_out = self.nt
        steps = np.linspace(0, self.nt - 1, num=self.num_out)
        self.saved_steps = np.around(steps).astype(int)
        self.saved_times = self.saved_steps*self.dt/3600

        if len(self.saved_steps) > 1e4:
            question = """
            The simulation will generate a large output file. Do you want to
            proceed ? [y/N]
            (Check "num_out" parameter in input file.)"""
            ans = input(question) or 'n'
            if ans != ('y' or 'Y'):
                raise ut.UserInputError('Canceled')

    def add_profiles(self):
        """
        Add initial composition profiles.

        Make atom fraction, vacancy site fraction and pore volume fraction
        profiles. These are then used to initialize cvar (see
        :class:`composition_variables.CompositionVariables`).

        """
        if ut.isfilename(self.x_profile):
            self.add_x_profile()
            self.check_x_profile_from_file()
        else:
            self.check_x_profile_from_dicts()
            self.add_x_profile()
        self.x_init = {k: self.x_init[k] for k in self.inds}
        self.add_yVa_profile()
        self.add_fm_profile()
        self.cvar = cv.CompositionVariables(self.comps, self.x_init,
                                            self.yVa_init, self.V_partial,
                                            self.fm_init)

    def add_x_profile(self):
        """Make initial atom fraction profile."""
        self.x_init, self.zstep = make_profile(self.x_profile,
                                               self.zmin, self.zmax,
                                               self.zm_init,
                                               val_left=self.x_left,
                                               val_right=self.x_right,
                                               sdir=self.work_dir)

    def check_x_profile_from_file(self):
        """
        Prepare initial atom fraction profile.

        Applies to profiles provided as file.

        * Make sure the constituents in the initial atom fraction profile match
          the list of independent constituents declared in the input file.
        * Enforce bounds on initial atom fractions and print warning.

        """
        try:
            assert list(self.x_init) == self.inds
        except AssertionError as exc:
            msg = ('Missing or extra element in initial atom fraction profile.'
                   f"Declared independent constituents: {self.inds}.\n"
                   f"\n{self.x_profile} contains profiles for "
                   f"{list(self.x_init)}.")
            raise ut.UserInputError(msg) from exc

        min_val = self.default_parameters['min_atom_fraction']
        max_val = 1 - self.default_parameters['min_atom_fraction']
        for k in self.x_init:
            prof = self.x_init[k]
            if any(prof < min_val) or any(prof > max_val):
                msg = (f"Some atom fractions in {self.x_profile} were replaced"
                       f" by minimum allowed {min_val}"
                       f" or maximum allowed {max_val}.")
                self.logger.info(msg)
                prof[:] = np.clip(prof, min_val, max_val)

    def check_x_profile_from_dicts(self):
        """
        Prepare initial atom fraction profile.

        Applies to profiles provided as predefined types (see list in
        :func:`make_profile`).

        * Make sure the constituents in the initial atom fraction profile match
          the list of independent constituents declared in the input file.
        * Enforce bounds on initial atom fractions and print warning.

        """
        min_val = self.default_parameters['min_atom_fraction']
        max_val = 1 - self.default_parameters['min_atom_fraction']
        varlist = ['x_left', 'x_right'] if self.x_right else ['x_left']
        for var in varlist:
            x_side = getattr(self, var)
            for k in x_side:
                if x_side[k] < min_val:
                    msg = (f"Input {var} {k} = {x_side[k]} replaced "
                           f"by minimum allowed {min_val}.")
                    self.logger.info(msg)
                    x_side[k] = min_val
                if x_side[k] > max_val:
                    msg = (f"Input {var} {k} = {x_side[k]} replaced "
                           f"by maximum allowed {max_val}.")
                    self.logger.info(msg)
                    x_side[k] = max_val

        try:
            assert set(self.x_left) == set(self.inds)
            if self.x_right:
                assert set(self.x_right) == set(self.inds)
        except AssertionError as exc:
            msg = ("Missing or extra element(s) in initial atom fraction "
                   "profile.\n"
                   f"Declared independent constituents: {self.inds}.\n"
                   "Constituents in initial profile:\n"
                   f"x_left: {list(self.x_left)}")
            if self.x_right:
                msg += f"\nx_right: {list(self.x_right)}"
            raise ut.UserInputError(msg) from exc

    def add_yVa_profile(self):
        """Make initial vacancy site fraction profile."""
        if self.yVa_profile is not None:
            yVa_left = {'yVa': self.yVa_left}
            yVa_right = {'yVa': self.yVa_right}
            yVa_init, _ = make_profile(self.yVa_profile, self.zmin, self.zmax,
                                       self.zm_init,
                                       val_left=yVa_left, val_right=yVa_right,
                                       sdir=self.work_dir)
            self.yVa_init = yVa_init['yVa']
        else:
            x_init_arr = np.array([self.x_init[k] for k in self.inds])
            self.yVa_init = self.yVa_fun(x_init_arr)

    def add_fm_profile(self):
        """Make initial metal volume fraction profile."""
        if self.fp_profile is not None:
            fp_left = {'fp': self.fp_left}
            fp_right = {'fp': self.fp_right}
            fp_init, _ = make_profile(self.fp_profile, self.zmin, self.zmax,
                                      self.zm_init,
                                      val_left=fp_left, val_right=fp_right,
                                      sdir=self.work_dir)
            self.fm_init = 1 - fp_init['fp']
        else:
            self.fm_init = np.ones(self.nz_init - 1)

    def add_BC(self):
        """
        Add boundary conditions.

        The keys of the BC attribute are 'left', 'right', 'cvar_left',
        'cvar_right', 'J_left', 'J_right'.

        'left' and 'right' indicate the type of BC and can be either:

        * 'Dirichlet': prescribed composition, if 'xBC_left' or 'xBC_right'
          is provided in params.
        * 'Neumann': prescribed flux, if 'JBC_left' or 'JBC_right' is
          provided in params.

        The other keys ('cvar_left', 'cvar_right', 'J_left', 'J_right') are
        functions of time (ex: (3*t + 2)**(1/2)) that return composition
        variables and flux arrays.

        Defaults to 0-flux BC if no BC is specified in the input file.

        """
        self.BC = {}
        info_list = []

        for side in ['left', 'right']:
            self.BC[f'{side}'] = self.params.get(f'BC_{side}_type', 'Neumann')
            if self.BC[f'{side}'] == 'Dirichlet':
                x_strings = self.params[f'xBC_{side}']
                self.BC[f'c_{side}'] = self.make_cBC_fun(x_strings)
            else:
                J_strings = self.params.get(f'JBC_{side}', {})
                for k in self.comps[1:]:
                    if k not in J_strings:
                        J_strings[k] = "0"
                        info_list.append((side, k))
                self.BC[f'J_{side}'] = self.make_JBC_fun(J_strings)

        if len(info_list) > 0:
            self.logger.info("Auto boundary conditions:")
            for side, k in info_list:
                text = f"* {side:5} BC for {k} set to 0 flux"
                self.logger.info(text)

    def make_JBC_fun(self, J_strings):
        """
        Make function that returns a boundary flux array.

        The inner function returns a 1D array of fluxes, which include
        all atom constituents.

        """
        funs = {k: lambda t, s=J_strings[k]: eval(s) for k in self.comps[1:]}

        def fun(t):
            J_dict = {k: funs[k](t) for k in self.comps[1:]}
            J_arr = np.array(list(J_dict.values()))
            return J_arr

        return fun

    def make_cBC_fun(self, x_strings):
        """
        Make function that returns a boundary composition variable.

        The inner function returns a
        :class:`composition_variables.CompositionVariables` instance. Two
        assumptions are made:

        * vacancies are at equilibrium,
        * the pore fraction is 0.

        """
        fm = 1

        def fun(t, s=x_strings):
            x_dict = {k: eval(s[k]) for k in self.inds}
            x_arr = np.array([x_dict[k] for k in self.inds])
            yVa = self.yVa_fun(x_arr)
            cvar = cv.CompositionVariables(self.comps, x_dict, yVa,
                                           self.V_partial, fm)
            return cvar

        return fun

    def prepare_simulation_log(self):
        """Add simulation info to log."""
        self.logger.data(f'nt = {self.nt}')
        self.logger.data(f'dt = {self.dt}')
        self.logger.data(f'saved_steps = {self.saved_steps.tolist()}')
        self.logger.data(f'saved_times = {self.saved_times.tolist()}')
        self.logger.data(f'zstep = {self.zstep}')
        self.logger.info('Running simulation')

    def run(self, show_completion=False, verbose=1):
        """
        Run diffusion simulation.

        * Call to :func:`solvers.solver`. The function returns variables at
          saved_steps.
        * These are printed to files, and stored in res dict. The keys of res
          are the time steps. The last time step can be accessed with -1 as
          key.
        * A mass balance is performed (see :meth:`add_mass_balance`).

        Parameters
        ----------
        show_completion : bool, optional
            Print completion rate while simulation is running.
        verbose : int, optional
            Verbosity level, sets amount of information printed while
            simulation is running. Recommended values: 0 (less verbose) and 1
            (more verbose). See :func:`solvers.solver` and
            :func:`solvers.remesh`.

        """
        self.prepare_simulation_log()

        resdict = so.solver(self.z_init, self.geometry, self.cvar,
                            self.MU_funy, self.L_fun, self.yVa_fun,
                            self.nt, self.dt,
                            self.ideal_lattice,
                            self.k_dislo, self.k_pores,
                            self.rho_dislo, self.rho_pores,
                            self.DVa_fun,
                            self.Va_method,
                            self.saved_steps, self.BC, show_completion,
                            verbose,
                            self.stencil,
                            self.logger)

        self.logger.results(resdict)
        self.simres.add_results(resdict)


class NewSimulation(Simulation):
    """
    Create new simulation from input file.

    This is the recommended way to create a new simulation. See
    :class:`Simulation` for documentation on attributes and methods.

    """

    def __init__(self, ref):
        """
        Class constructor.

        Call :class:`Simulation` constructor to define system parameters and
        pre-process simulation conditions. Then call :class:`Simulation`
        methods to complete simulation setup.

        """
        super().__init__(ref)
        self.logger.info("Adding space grid.")
        self.add_grid()
        self.logger.info("Adding initial composition profiles.")
        self.add_profiles()
        self.logger.info("Adding boundary conditions.")
        self.add_BC()
        self.logger.info("Adding time steps.")
        self.add_time()
        self.add_saved_steps()

        # Initialize results
        simu_parameters = {'comps': self.comps,
                           'V_partial': self.V_partial,
                           'title': self.title,
                           'saved_times': self.saved_times}
        self.simres = results.SimulationResults(simu_parameters)
        self.results = self.simres.results  # shortcut
        self.res = self.results  # backward compatibility
        self.result = self.simres.result  # shortcut

        # Shortcuts to useful plot methods
        self.plot = self.simres.plot
        self.plot_quartet = self.simres.plot_quartet
        self.interactive_plot = self.simres.interactive_plot


class ReadSimulationOld(Simulation):
    """
    Create simulation from results file.

    See :class:`Simulation` for documentation on attributes and methods.

    """

    def __init__(self, ref):
        """
        Class constructor.

        Read input parameters, log and simulation results from output file.
        Check that the parameters in the output file coincide with those in the
        input file. A difference would mean that the input file has been
        modified after the simulation was run, and would be error-inducing.
        Then add attributes to complete initialization.

        """
        super().__init__(ref)

        read_pstring, log, resdict = rs.read_res_from_file(ref, self.work_dir)

        try:
            assert read_pstring == self.pstring
        except AssertionError:
            info = "Warning: the set of parameter read in the result file "
            info += f"differs from that found in {ref}-input.txt"
            print(info)

        log_dict = pa.read_log_string(log)
        saved_steps = list(resdict.keys())
        dt = float(log_dict['dt'])
        saved_times = np.array([step*dt/3600 for step in saved_steps])
        simu_parameters = {'comps': self.comps,
                            'V_partial': self.V_partial,
                            'title': self.title,
                            'saved_times': saved_times}
        self.simres = results.SimulationResults(simu_parameters)
        self.simres.add_results(resdict)
        self.results = self.simres.results  # shortcut
        self.res = self.results  # backward compatibility
        self.result = self.simres.result  # shortcut

        # Shortcuts to useful plot methods
        self.plot = self.simres.plot
        self.plot_quartet = self.simres.plot_quartet
        self.interactive_plot = self.simres.interactive_plot


class ReadSimulation:
    """
    Create simulation from results file.

    See :class:`Simulation` for documentation on attributes and methods.

    """

    def __init__(self, ref):
        """
        Class constructor.

        Read log file of a previous simulation, extract useful attributes (see
        list in :class:`alloy_system.AlloySystem` and :class:`Simulation`) and
        make instance of :class:`results.SimulationResults`. The latter
        organizes the results and exposes plot methods. The most useful plot
        methods are added as attributes of this class for ease of access (and
        backward compatibility).

        """
        print(intro_msg(), flush=True)
        now = datetime.datetime.today()
        now = now.strftime('%Y-%m-%d %H:%M:%S')
        print(f"INFO   : {now} Reading '{ref}' results.")

        self.work_dir = Path.cwd()
        data, resdict = rs.parse_log(ref, self.work_dir)

        self.comps = data['comps']
        self.inds = self.comps[1:-1]
        self.V_partial = data['V_partial']
        self.title = data['title']
        self.saved_steps = list(resdict.keys())
        self.saved_times = np.array(data['saved_times'])
        self.nt = int(data['nt'])
        self.dt = data['dt']
        self.th = data['th']
        self.zstep = data['zstep']
        self.thermo_db = data['thermo_db']
        self.mob_db = data['mob_db']
        simu_parameters = {'comps': self.comps,
                           'V_partial': self.V_partial,
                           'title': self.title,
                           'saved_times': self.saved_times}
        self.simres = results.SimulationResults(simu_parameters)
        self.simres.add_results(resdict)
        self.results = self.simres.results  # shortcut
        self.res = self.results  # backward compatibility
        self.result = self.simres.result  # shortcut

        # Shortcuts to useful plot methods
        self.plot = self.simres.plot
        self.plot_quartet = self.simres.plot_quartet
        self.interactive_plot = self.simres.interactive_plot


def make_profile(profile, zmin, zmax, zm_init,
                 val_left=None, val_right=None, sdir=None):
    """
    Make profile along distance axis.

    Behavior determined by the 'profile' parameter:

    * 'step' followed by a float argument: step profile with step at zstep (see
      doc in :func:`meshing.make_step_profile`) and values from val_left and
      val_right.

      If the argument is > 0.1, it is interpreted as a fraction of zrange
      -> zstep = zmin + argument*(zmax - zmin)

      If the argument is < 0.1, it is interpreted as zstep in m.
    * 'smooth step' followed by a float argument: same with an error function
      instead of heaviside (see :func:`meshing.make_smooth_step_profile`)
    * 'flat': flat profile with values from val_left.
    * filename: read from sdir/filename using np.genfromtxt

    The input file must have the constituent profiles arranged by
    columns, with the constituent names on the first line. The size
    of the columns must match that of zm (i.e., `nz - 1`).

    The 'val_left' and 'val_right' parameters are dicts with the independent
    constituents as keys.
    """
    prof = None
    zstep = None
    if profile[0] in ('step', 'smooth step'):
        arg = profile[1]
        if arg > 0.1:
            zstep = zmin + arg * (zmax - zmin)
        else:
            zstep = arg
        if profile[0] == 'step':
            prof = mesh.make_step_profile(zm_init, zstep, val_left, val_right)
        elif profile[0] == 'smooth step':
            prof = mesh.make_smooth_step_profile(zm_init, zstep, val_left,
                                                 val_right)
    elif profile[0] == 'flat':
        prof = {k: np.ones(zm_init.shape)*v for k, v in val_left.items()}
    else:
        fpath = sdir / profile
        arr = np.genfromtxt(fpath, names=True, delimiter=',')
        try:
            assert arr.shape[0] == zm_init.size
        except AssertionError as exc:
            msg = f"Initial atom fraction profile provided in {profile} has "
            msg += f"size {arr.shape[0]}. It is not compatible with initial "
            msg += f"grid of size {zm_init.size} (nz = {zm_init.size + 1} "
            msg += f"nodes i.e. {zm_init.size} volumes)."
            raise ut.UserInputError(msg) from exc
        prof = {k: arr[k] for k in arr.dtype.names}
    return prof, zstep


def make_compo_string_from_xdict(x, dep):
    """
    Make label indicating an alloy composition.

    Parameters
    ----------
    x : dict
        Atom fractions of independent constituents.
    dep : str
        Name of dependent constituent.

    Returns
    -------
    A : str
        Label indicating an alloy composition.

    Examples
    --------
    >>> x = {'Cr': 0.2, 'Al': 0.062}
    >>> dep = 'Ni'
    >>> make_compo_string_from_xdict(x, dep)
    'Ni-20Cr-6.2Al'

    """
    s = {}
    for k in x:
        N = round(x[k]*100, 1)
        s[k] = f"{N:.0f}" if N.is_integer() else f"{N:.1f}"
    A = dep + ''.join(f'-{v}{k}' for k, v in s.items() if float(v) != 0)
    return A


def make_plot_title(dep, profile_type, x_left, x_right, TC):
    """
    Make plot title.

    Parameters
    ----------
    dep : str
        Name of dependent constituent.
    profile_type : str
        Type of initial atom fraction profile. See :func:`make_profile`.
    x_left : dict of floats
        Initial atom fraction of dependent constituents on left hand side.
    x_right : dict of floats
        Initial atom fraction of dependent constituents on right hand side.
    TC : float
        Temperature in Celsius.

    Returns
    -------
    title : str
        Plot title.

    """
    if profile_type[0] == 'step':
        A_left = make_compo_string_from_xdict(x_left, dep)
        A_right = make_compo_string_from_xdict(x_right, dep)

        title = rf'{A_left} vs. {A_right} at {int(TC)} $^\circ$C'

    elif profile_type[0] == 'flat':
        A_left = make_compo_string_from_xdict(x_left, dep)
        title = rf'{A_left} at {int(TC)} $^\circ$C'

    else:
        title = f'No title implemented for "{profile_type}" profile'

    return title
