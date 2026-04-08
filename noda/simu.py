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
        Names and file paths of databases in data folder.
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
        Databases used in the simulation (thermo, mobility,
        molar_volume, vacancy_formation_energy)
    V_partial : dict
        Partial molar volumes.
    thermo : :class:`thermodynamics.Thermodynamics`
        Thermodynamic properties.
    mobility : :class:`mobility.Mobility`
        Mobility properties.
    space : :class:`space.SpaceGrid`
        Space grid parameters.
    boundary_conditions : dict
        Instances of :class:`boundary_conditions.BoundaryConditions`, with
        'left' and 'right' as keys.
    lattice : :class:`lattice.Lattice`
        Parameters relative to vacancy annihilation/creation.
    initial_conditions : :class:`initial_conditions.InitialConditions`
        Initial conditions.
    time : :class:`time.TimeGrid`
        Time-related parameters.
    L_mean_kind : str
        Kind of mean used to compute L values at nodes.
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
        self.user_data = da.get_user_data(self.data_dir, self.logger)
        self.db_register = self.get_database_register()
        self.default_parameters = self.get_default_parameters()
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
        self.V_partial = da.get_partial_molar_volume(
                            self.databases,
                            self.db_register['partial_molar_volume'],
                            self.comps,
                            self.default_parameters['partial_molar_volume'],
                            logger)
        self.thermo = self.get_thermo_handler()
        self.mobility = self.get_mob_handler()
        if 'space' in config:
            self.space = SpaceGrid(config['space'],
                                   self.default_parameters,
                                   work_dir,
                                   logger)
        options = config.get('options', {})
        self.lattice = Lattice(options, work_dir, logger)
        if self.lattice.ideal is False:
            self.thermo.set_nonideal_lattice()

        # Need to handle after lattice because depends on thermo.ideal_lattice
        if 'initial_conditions' in config:
            self.initial_conditions = InitialConditions(
                                          config['initial_conditions'],
                                          self.V_partial,
                                          self.space,
                                          work_dir,
                                          min_atom_fraction,
                                          self.thermo,
                                          logger)
        if 'time' in config:
            self.time = TimeGrid(config['time'],
                                 self.initial_conditions.x,
                                 self.space.dz_init,
                                 self.mobility.DT_funx,
                                 self.default_parameters,
                                 logger)
        
        # Flag ready if tables required for run are present in config.
        # This is used in BoundaryConditions to determine whether auto-boundary
        # conditions INFO should be streamed.
        cats = ['space', 'initial_conditions', 'time']
        self.ready = all(hasattr(self, cat) for cat in cats)
        self.boundary_conditions = self.get_boundary_conditions(min_atom_fraction)
        self.L_mean_kind = options.get('L_mean_kind',
                                       self.default_parameters['L_mean_kind'])

        # TODO : check that all input config entries are valid. Difficulty:
        # nested structure both in config and in instance attributes

        # Initialize results if simulation is ready to run
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
        possible_keys = ['thermo', 'thermodynamics']
        if not any(k in self.config['databases'] for k in possible_keys):
            msg = ("Missing required 'thermo' or 'thermodynamics' entry in "
                   "input dict, in 'databases' subdict.")
            raise ut.UserInputError(msg) from None
        if 'mobility' not in self.config['databases']:
            msg = ("Missing required 'mobility' entry in input dict, in "
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

    def get_database_register(self):
        """
        Get database register from 'user_data.toml'.

        Four categories of databases are handled:

        * Partial molar volume
        * Vacancy formation energy
        * Thermodynamics
        * Mobility

        All are optional : if no entry in 'user_data.toml', default to empty
        dict.

        """
        db_register = {}
        for cat in ['partial_molar_volume',
                    'vacancy_formation_energy',
                    'thermodynamics',
                    'mobility']:
            dct = self.user_data.get(cat, {})
            db_register[cat] = {k.lower(): val for k, val in dct.items()}
        return db_register

    def get_default_parameters(self):
        """
        Get default parameters.

        For each parameter defined in
        :data:`constants.factory_default_parameters`, look for entry in
        'user_data.toml', and default to entry in
        :data:`constants.factory_default_parameters`.

        """
        factory = co.factory_default_parameters
        user = self.user_data.get('default_parameters', {})
        default_parameters = {k: user.get(k, val)
                              for k, val in factory.items()}
        return default_parameters

    def get_thermo_handler(self):
        """
        Process input parameters and make thermodynamic properties handler.

        Returns
        -------
        :class:`thermodynamics.Thermodynamics`
            Thermodynamic properties handler.

        """
        possible_keys = ['thermo', 'thermodynamics']
        key = (self.databases.keys() & possible_keys).pop()
        name = self.databases[key].lower()
        if Path(name).is_file():
            fpath = Path(name)
        else:
            db = ut.get_or_raise(self.db_register['thermodynamics'], name)
            fpath = self.data_dir / db["file"]
        params = da.get_thermo_from_file(fpath,
                                         self.phase,
                                         self.comps[1:],
                                         self.TK,
                                         self.logger)
        msg = f"Reading thermodynamic data in '{fpath.resolve()}'."
        self.logger.info(msg)
        GfV = da.get_vacancy_formation_energy(
                        self.databases,
                        self.db_register['vacancy_formation_energy'],
                        self.comps,
                        self.default_parameters['vacancy_formation_energy'],
                        self.logger)
        return Thermodynamics(params,
                              self.comps,
                              self.phase,
                              self.TK,
                              GfV,
                              self.logger)

    def get_mob_handler(self):
        """
        Process input parameters and make mobility properties handler.

        Returns
        -------
        :class:`mobility.Mobility`
            Mobility properties handler.

        """
        name = self.databases['mobility'].lower()
        if Path(name).is_file():
            fpath = Path(name)
        else:
            db = ut.get_or_raise(self.db_register['mobility'], name)
            fpath = self.data_dir / db["file"]
        params = da.get_mob_from_file(fpath, self.comps[1:], self.TK,
                                      self.logger)
        msg = f"Reading mobility data in '{fpath.resolve()}'."
        self.logger.info(msg)
        return Mobility(params, self.comps[1:], self.TK)

    def get_boundary_conditions(self, min_atom_fraction):
        """
        Make dict of :class:`boundary_conditions.BoundaryConditions` instances.

        Parameters
        ----------
        min_atom_fraction : min_atom_fraction : float
            Minimum atom fraction accepted.

        Returns
        -------
        dct : dict
            :class:`boundary_conditions.BoundaryConditions` instances.

        """
        params = self.config.get('boundary_conditions', {})
        dct = {}
        for side in ['left', 'right']:
            side_params = params.get(side, {})
            dct[side] = BoundaryConditions(side_params,
                                           self.thermo,
                                           self.V_partial,
                                           min_atom_fraction,
                                           self.logger,
                                           side,
                                           self.ready)
        return dct


class NewSimulation(Simulation):
    """
    Create new simulation from input file or dict.

    Process user input and make :class:`Simulation` instance. The user is
    expected to provide either a file or a dict, using the 'file' or 'config'
    keyword argument.

    See :class:`Simulation` for documentation on attributes and methods.

    """
    def __init__(self, file=None, config=None, ref=None, log=True):
        """Class constructor."""
        if log:
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
        logger = lut.CustomLogger(work_dir, ref, log)
        logger.info(f"{now} Creating '{ref}' simulation.")
        log_fpath = logger.handlers[1].baseFilename if log else None
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
        utils.UserInputError
            If simulation is not ready.

        """
        if not self.ready:
            for cat in ['space', 'initial_conditions', 'time']:
                if cat not in self.config:
                    msg = (f"Simulation is not ready to run. Required '{cat}' "
                           "parameters are missing from input.")
                    raise ut.UserInputError(msg)

        self.prepare_simulation_log()
        resdict = so.solver(self.thermo,
                            self.mobility,
                            self.space,
                            self.initial_conditions,
                            self.boundary_conditions,
                            self.time,
                            self.lattice,
                            show_completion,
                            verbose,
                            self.L_mean_kind,
                            self.logger)
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
