# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define and load diffusion simulation."""

from pathlib import Path
import datetime
import tomllib
import textwrap as tw

from noda import paths
import noda.log_utils as lut
import noda.data_io as da
import noda.constants as co
import noda.solvers as so
import noda.utils as ut
from noda.thermodynamics import Thermodynamics
from noda.mobility import Mobility
from noda.space import SpaceGrid
from noda.initial_conditions import InitialConditions
from noda.boundary_conditions import BoundaryConditions
from noda.time_grid import TimeGrid
from noda.temperature import Temperature
from noda.lattice import Lattice
from noda import results
from noda import __version__


class Simulation:
    """
    Handle system properties, simulation conditions, simulation results.

    This is a base class that organizes all input parameter processing
    and prepares simulation results. It is meant to be used only through
    subclasses :class:`NewSimulation` (create new simulation) or
    :class:`ReadSimulation` (read simulation from log file), not directly by
    the user.

    This is a 1D version, with constant temperature, time-invariant space grid
    and time step.

    Attributes
    ----------
    work_dir : pathlib.Path
        Work directory.
    logger : :class:`log_utils.CustomLogger`
        Logger.
    data_dir : pathlib.Path
        Data directory, user- or package-provided (see
        :func:`paths.get_data_dir`).
    config : dict
        Simulation input parameters.
    db_register : dict
        Names and file paths of available databases.
    default_parameters : dict
        Parameters used when not specified in user input.
    comps : list of str
        System constituents (Va + atom components).
    inds : list of str
        Independent constituents (Va + atom components).
    phase : str
        System phase name.
    temperature : :class:`temperature.Temperature`
        Temperature program.
    TK : float
        Temperature in Kelvin
    TC : float
        Temperature in Celsius
    databases : dict
        Names of user-specified databases (thermo, mobility, optional :
        molar_volume, vacancy_formation_energy)
    volume_db : str
        Name of partial molar volume database.
    V_partial : dict
        Partial molar volumes.
    thermo : :class:`thermodynamics.Thermodynamics`
        Thermodynamic properties.
    mob : :class:`mobility.Mobility`
        Mobility properties.
    space : :class:`space.SpaceGrid`
        Space grid parameters.
    BC : :class:`boundary_conditions.BoundaryConditions`
        Boundary conditions.
    lattice : :class:`lattice.Lattice`
        Parameters relative to vacancy annihilation/creation.
    init : :class:`initial_conditions.InitialConditions`
        Initial conditions.
    time : :class:`time.TimeGrid`
        Time-related parameters.
    stencil : str
        Space discretization stencil.
    ready : bool
        Whether simulation is ready to run.
    simres : :class:`results.SimulationResults`
        Simulation results.
    results : dict
        Simulation results, nested dicts with syntax res[th][var][k], see
        :class:`results.SimulationResults`.
    res : dict
        Alias of results.
    result :
        Function to access simulation results, see
        :meth:`results.SimulationResults.result`.
    """

    def __init__(self, config, work_dir, logger):
        """Class constructor."""
        self.work_dir = work_dir
        self.logger = logger
        self.data_dir = paths.get_data_dir(work_dir, logger)
        (volume_databases, vacancy_databases,
         self.db_register, self.default_parameters) = self.get_user_data()
        min_atom_fraction = self.default_parameters['min_atom_fraction']
        self.config = config
        logger.input(config)
        self.check_required_parameters()
        self.comps = self.get_constituents()
        self.inds = self.comps[1:-1]
        self.phase = config['system']['phases']
        self.temperature = Temperature(config['temperature'])
        self.TC = self.temperature.TC
        self.TK = self.temperature.TK
        self.databases = config['databases']
        vacancy_db = self.get_and_log('vacancy', stream=False)
        self.volume_db = self.get_and_log('molar_volume')
        self.V_partial = da.get_volume_data(volume_databases,
                                            self.volume_db,
                                            self.comps,
                                            logger)
        self.thermo = self.get_thermo_handler(vacancy_databases, vacancy_db)
        self.mob = self.get_mob_handler()
        if 'space' in config:
            self.space = SpaceGrid(config['space'],
                                   self.default_parameters,
                                   work_dir,
                                   logger)
        if 'boundary_conditions' in config:
            self.BC = self.get_boundary_conditions(
                                  config['boundary_conditions'],
                                  min_atom_fraction)
        options = config.get('options', {})
        self.lattice = Lattice(options, work_dir, logger)
        if self.lattice.ideal is False:
            self.thermo.set_nonideal_lattice()

        # Need to handle after lattice because depends on thermo.ideal_lattice
        if 'initial_conditions' in config:
            self.init = InitialConditions(config['initial_conditions'],
                                          self.V_partial,
                                          self.space,
                                          work_dir,
                                          min_atom_fraction,
                                          self.thermo,
                                          logger)
        if 'time' in config:
            self.time = TimeGrid(config['time'],
                                 self.init.x,
                                 self.space.dz_init,
                                 self.mob.DT_fun,
                                 self.default_parameters,
                                 logger)

        self.stencil = options.get('stencil', self.default_parameters['stencil'])
        if self.stencil not in ['H', 'A', 'G']:
            msg = f'Stencil "{self.stencil}" not implemented'
            raise ut.UserInputError(msg) from None

        # TODO : check that all input config entries are valid. Difficulty:
        # nested structure both in config and in instance attributes

        # Initialize results if simulation is ready to run
        cats = ['space', 'init', 'time', 'BC']
        self.ready = all(hasattr(self, cat) for cat in cats)
        if self.ready:
            params = {'comps': self.comps,
                      'V_partial': self.V_partial,
                      'saved_th': self.time.saved_th}
            self.simres = results.SimulationResults(params)
            self.results = self.simres.results  # shortcut
            self.res = self.results  # shortcut shortcut
            self.result = self.simres.result  # shortcut

            # Shortcuts to useful plot methods
            self.plot = self.simres.plot
            self.plot_quartet = self.simres.plot_quartet
            self.interactive_plot = self.simres.interactive_plot

    def check_required_parameters(self):
        """Make sure required parameters are present in input dict."""
        for cat in ['databases', 'system', 'temperature']:
            if cat not in self.config:
                msg = f"Missing required '{cat}' entry in input dict."
                raise ut.UserInputError(msg) from None
        for cat in ['thermo', 'mobility']:
            if cat not in self.config['databases']:
                msg = (f"Missing required '{cat}' entry in input dict, in "
                       "'databases' subdict.")
                raise ut.UserInputError(msg) from None
        for cat in ['components', 'phases']:
            if cat not in self.config['system']:
                msg = (f"Missing required '{cat}' entry in input dict, in "
                       "'system' subdict.")
                raise ut.UserInputError(msg) from None

    def get_constituents(self):
        """Get atom components from input dict and add vacancies."""
        raw = self.config['system']['components']
        comps = [x.strip() for x in raw.split(',')]
        comps = [ut.format_element_symbol(x) for x in comps]
        inds = comps[1:]
        return ['Va'] + inds + [comps[0]]

    def get_user_data(self):
        """
        Get data from 'user_data.toml' file.

        * Molar volume databases : required
        * Vacancy formation energy databases : required
        * Register of thermodynamics and mobility databases : required
        * Default parameters : optional : for each parameter, defaults to
          value in package-provided `co.factory_default_parameters`.

        """
        user_data = da.get_user_data(self.data_dir, self.logger)
        volume_databases = ut.get_or_raise(user_data, 'molar_volume')
        vacancy_databases = ut.get_or_raise(user_data,
                                            'vacancy_formation_energy')
        thermo_register = ut.get_or_raise(user_data, 'thermodynamics')
        mobility_register = ut.get_or_raise(user_data, 'mobility')
        # Merge and make all keys lowercase
        db_register = thermo_register | mobility_register
        db_register = {k.lower(): val for k, val in db_register.items()}
        # Get numerical parameters from user data or from constants module
        factory = co.factory_default_parameters
        if 'default_parameters' in user_data:
            user_dct = user_data['default_parameters']
            default_parameters = {k: user_dct.get(k, val)
                                  for k, val in factory.items()}
        else:
            default_parameters = factory
        return (volume_databases, vacancy_databases,
                db_register, default_parameters)

    def get_and_log(self, key, stream=True):
        """
        Get item from self.databases, defaults to default_parameters.

        Wrapper around :func:`log_utils.get_and_log`.

        Parameters
        ----------
        key : str
            Key of item in self.databases dict.
        stream : bool, optional
            Log to screen in addition to file. The default is True.

        Returns
        -------
        res :
            Item of interest.

        """
        res = lut.get_and_log(self.databases,
                              key,
                              self.default_parameters[f'{key}_database'],
                              self.logger,
                              stream=stream)
        return res

    def get_thermo_handler(self, vacancy_databases, vacancy_db):
        """
        Process input parameters and make thermodynamic properties handler.

        Parameters
        ----------
        vacancy_databases : dict
            Vacancy formation energy databases.
        vacancy_db : str
            Name of vacancy formation energy database.

        Returns
        -------
        :class:`thermodynamics.Thermodynamics`
            Thermodynamic properties handler.

        """
        db_name = self.databases['thermo'].lower()
        db = ut.get_or_raise(self.db_register, db_name)
        fpath = self.data_dir / db["file"]
        params = da.get_thermo_from_file(fpath, self.phase,
                                         self.comps[1:],
                                         self.TK,
                                         self.logger)
        db_file = ut.get_or_raise(self.db_register[db_name], 'file')
        msg = ("Thermodynamic functions built from database file "
               f"('{db_file}').")
        self.logger.info(msg)
        return Thermodynamics(params, self.comps, self.phase, self.TK,
                              vacancy_databases, vacancy_db, self.logger)

    def get_mob_handler(self):
        """
        Process input parameters and make mobility properties handler.

        Returns
        -------
        :class:`mobility.Mobility`
            Mobility properties handler.

        """
        db_name = self.databases['mobility'].lower()
        db = ut.get_or_raise(self.db_register, db_name)
        fpath = self.data_dir / db["file"]
        params = da.get_mob_from_file(fpath, self.comps[1:], self.TK,
                                      self.logger)
        db_file = ut.get_or_raise(self.db_register[db_name], 'file')
        msg = f"Mobility functions built from database file ('{db_file}')."
        self.logger.info(msg)
        return Mobility(params, self.comps[1:], self.TK)

    def get_boundary_conditions(self, params, min_atom_fraction):
        """
        Make dict of :class:`boundary_conditions.BoundaryConditions` instances.

        Parameters
        ----------
        params : dict
            Boundary conditions-related input parameters.
        min_atom_fraction : min_atom_fraction : float
            Minimum atom fraction accepted.

        Returns
        -------
        dct : dict
            :class:`boundary_conditions.BoundaryConditions` instances.

        """
        dct = {}
        for side in ['left', 'right']:
            side_params = params.get(side, {})
            dct[side] = BoundaryConditions(side_params,
                                           self.thermo,
                                           self.V_partial,
                                           min_atom_fraction,
                                           self.logger,
                                           side)
        return dct


class NewSimulation(Simulation):
    """
    Create new simulation from input file or dict.

    Process user input and make :class:`Simulation` instance. The user is
    expected to provide either a file or a dict, using the 'file' or 'config'
    keyword argument.

    See :class:`Simulation` for documentation on attributes and methods.

    """
    def __init__(self, file=None, config=None, ref=None):
        """Class constructor."""
        print(intro_msg(), flush=True)
        now = get_timedate()
        if ( (config is None and file is None)
            or (config is not None and file is not None) ):
            msg = ("Expecting either a file path ('file') or a dictionary "
                   "('config').")
            raise ut.UserInputError(msg)
        if file:
            input_file = Path(file).resolve()
            try:
                with open(input_file, mode="rb") as fp:
                    config = tomllib.load(fp)
            except FileNotFoundError as exc:
                msg = f"Input file '{input_file}' not found."
                raise ut.UserInputError(msg) from exc
            init_msg = f"Reading input file '{input_file}'."
            work_dir = input_file.parent
            if ref is None:
                ref = input_file.stem
        else:
            work_dir = Path().cwd().resolve()
            if ref is None:
                ref = get_timedate(file_format=True)
            init_msg = "Reading input dictionary."
        logger = lut.CustomLogger(work_dir, ref)
        logger.info(f"{now} Creating '{ref}' simulation.")
        log_fpath = logger.handlers[1].baseFilename
        logger.info(f"Log saved in '{log_fpath}'.")
        logger.info(init_msg)
        for x in config:
            # TODO : this should search through nested dicts
            if list(config).count(x) > 1:
                msg = f"Parameter {x} found more than once in input file."
                raise ut.UserInputError(msg) from None
        super().__init__(config, work_dir, logger)

    def prepare_simulation_log(self):
        """Add simulation info to log."""
        self.logger.data(f'nt = {self.time.nt}')
        self.logger.data(f'dt = {self.time.dt}')
        self.logger.data(f'saved_steps = {self.time.saved_steps.tolist()}')
        self.logger.data(f'saved_th = {self.time.saved_th.tolist()}')
        self.logger.info('Running simulation')

    def run(self, show_completion=False, verbose=1):
        """
        Prepare log and run diffusion simulation.

        * Call to :func:`solvers.solver`. The function returns variables at
          saved_steps.
        * These are logged to file, and stored in simres attribute
          (:class:`results.SimulationResults` instance).

        Parameters
        ----------
        show_completion : bool, optional
            Print completion rate to screen while simulation is running.
        verbose : int, optional
            Verbosity level, sets amount of information printed while
            simulation is running. Valid values: 0 (less verbose) and 1
            (more verbose). See :func:`solvers.solver` and
            :func:`solvers.remesh`.

        Raises
        ------
        ut.UserInputError
            If simulation is not ready.

        """
        if not self.ready:
            for cat in ['space', 'initial_conditions', 'time']:
                if cat not in self.config:
                    msg = (f"Simulation is not ready to run. Required '{cat}' "
                           "parameters are missing from input.")
                    raise ut.UserInputError(msg)

        self.prepare_simulation_log()
        resdict = so.solver(self.thermo, self.mob, self.space, self.init,
                            self.BC, self.time, self.lattice, show_completion,
                            verbose, self.stencil, self.logger)
        self.logger.results(resdict)
        self.simres.add_results(resdict)


class ReadSimulation(Simulation):
    """
    Create simulation from log file.

    See :class:`Simulation` for documentation on attributes and methods.

    """

    def __init__(self, file):
        """
        Class constructor.

        Read log file of a previous simulation, get input dict and simulation
        results and make instance of :class:`results.SimulationResults`.

        """
        print(intro_msg(), flush=True)
        fpath = Path(file).resolve()
        work_dir = fpath.parent
        ref = fpath.stem
        logger = lut.CustomLogger(work_dir, ref, log=False)
        now = get_timedate()
        print(f"INFO   : {now} Reading '{ref}' simulation.")
        try:
            config, resdict = lut.parse_log(fpath)
        except FileNotFoundError as exc:
            msg = f"Input file '{fpath}' not found."
            raise ut.UserInputError(msg) from exc
        msg = f"Reading log file '{fpath}'."
        for line in tw.wrap(msg, lut.MESSAGE_WIDTH):
            print("INFO   : " + line)
        super().__init__(config, work_dir, logger)
        self.simres.add_results(resdict)


def intro_msg():
    """Prepare intro message."""
    msg = '-'*35 + '  Noda  ' + '-'*36 + '\n'
    msg += '| A Python package for simulating diffusion in multicomponent '
    msg += 'alloys.' + ' '*9 + '|\n'
    msg += f'| Version {__version__}' + ' '*63 + '|\n'
    msg += '-'*79
    return msg

def get_timedate(file_format=False):
    """Get current timedate str."""
    now = datetime.datetime.today()
    if file_format:
        res = now.strftime('%Y-%m-%d_%Hh%Mm%S.%f')[:-3] + 's'
    else:
        res = now.strftime('%Y-%m-%d %H:%M:%S')
    return res
