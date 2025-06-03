# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Define thermodynamic and mobility properties of alloy system.
"""

from pathlib import Path
import datetime

import numpy as np

from noda.log_utils import get_and_log, CustomLogger
import noda.data_io as da
import noda.para_io as pa
import noda.thermo_functions as tfu
import noda.thermo_wrappers as wrap
import noda.utils as ut
import noda.constants as co
from noda import paths
from noda import __version__


class AlloySystem:
    """
    Thermodynamic and mobility properties of alloy system.

    This class mainly serves as a base class for :class:`simu.Simulation`. It
    can also be used on its own to examine thermodynamics and mobility data.

    This is an isothermal, single-phase version.

    Attributes
    ----------
    ref: str
        System reference, name of sub-directory with input and output files.
    work_dir: pathlib.Path
        Path to directory with system-related input and output files.
    database_register : dict
        Register of available thermodynamic and mobility databases.
    volume_databases : dict
        Available molar volume data.
    vacancy_databases : dict
        Available vacancy formation energy data.
    pstring : str
        Content of system input file, set of user-defined parameters.
    params: dict
        Set of system parameters, processed by
        :func:`para_io.read_parameter_string` from pstring.
    inds: list
        Independent constituents.
    dep: str
        Dependent constituent.
    comps: list
        Conservative (atomic) species, with dependent constituent at the end.
    phase: str
        Name of phase.
    thermo_db: str
        Name of thermodynamic database.
    mob_db: str
        Name of mobility database.
    TC: float
        Temperature in Celsius.
    TK: float
        Temperature in Kelvin.
    alloy: str
        Alloy system str representation, with dependent constituent first.
    stamp: str
        Label with syntax phase-alloy-TC.
    volume_db : str
        Name of partial molar volume database.
    V_partial : dict
        Partial molar volumes. See See :func:`data_io.get_volume_data`.
    vacancy_db : str
        Name of vacancy formation energy database.
    GfV : dict
        Vacancy formation energy in pure elements, given as [enthalpy, entropy]
        in J/mol and J/mol/K. See :func:`data_io.get_vacancy_formation_energy`.
    Va_method : str
        Method used to compute equilibrium vacancy fraction: 'cst'
        (composition-fixed) or 'RK' (composition-dependant via Redlich-Kister
        polynomial).
    thermo_params : dict of floats
        Thermodynamic parameters arranged as follows:

        | ``A: G_A for A in endmembers``
        | ``AB: [L0, L1] for AB in binary subsystems``

    mob_params : dict of dicts
        Mobility parameters arranged as follows:

        ``{i: subdict for i in comps}``

        subdict: ``{j: val for j in subsystems}``.

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
    ideal_lattice : bool
        Whether lattice is ideal, i.e., vacancy fraction is kept at
        equilibrium.
    G_funx : function
        Compute Gibbs free energy at composition provided as array of atom
        fractions of the dependent constituents.
    G_fun : function
        Alias of G_funx
    G_funy : function
        Compute Gibbs free energy at composition provided as array of site
        fractions of vacancies and of the dependent constituents. See
        :meth:`add_thermo_model`.
    MU_funx : function
        Compute chemical potential at composition provided as array of atom
        fractions of the dependent constituents.
    MU_fun : function
        Alias of MU_funx
    MU_funy: function
        Compute chemical potentials at composition provided as array of site
        fractions of vacancies and of the dependent constituents. See
        :meth:`add_thermo_model`.
    yVa_fun : function
        Compute equilibrium vacancy fraction at composition provided as array
        of atom fractions of the dependent constituents. See
        :meth:`add_Va_model`.
    DT_fun: function
        Compute tracer diffusion coefficients at composition provided as array
        of atom fractions of the dependent constituents.
    L_fun: function
        Compute phenomenological coefficients at composition provided as array
        of atom fractions of the dependent constituents. Derived from `DT_fun`
        (see :func:`thermo_functions.make_Lfun`).
    DVa_fun : function
        Compute vacancy diffusion coefficient from atom fractions and site
        fractions, see :func:`thermo_functions.make_DVa_fun`.
    """

    def __init__(self, ref):
        """
        Class constructor.

        Read input file and define system parameters. See
        :func:`para_io.read_parameter_string` for how input file is parsed.
        """
        self.logger = CustomLogger(ref)
        print(intro_msg(), flush=True)
        now = datetime.datetime.today()
        now = now.strftime('%Y-%m-%d %H:%M:%S')
        self.logger.info(f"{now} Creating '{ref}' system.")

        self.ref = ref
        self.work_dir = Path.cwd()
        self.data_dir = paths.get_data_dir(self.work_dir, self.logger)
        self.get_user_data()
        para_fpath = self.work_dir / (ref + '-input.txt')
        if not para_fpath.exists():
            msg = (f"Input file '{ref}-input.txt' not found in "
                   f"'{self.work_dir}'.")
            raise ut.UserInputError(msg) from None
        with open(para_fpath, encoding='utf-8') as file:
            raw_input = file.read()
        self.logger.info(f"Reading input file '{ref}-input.txt'.")
        self.pstring, self.params = pa.read_parameter_string(raw_input,
                                                             self.logger)
        self.logger.data('Begin input file.')
        for line in self.pstring.split('\n'):
            self.logger.data(line)
        self.logger.data('End input file.')
        self.inds = self.params['inds']
        self.dep = self.params['dep']
        self.comps = ['Va'] + self.inds + [self.dep]
        self.phase = self.params['phase']
        self.thermo_db = self.params['thermo_db'].lower()
        self.mob_db = self.params['mob_db'].lower()
        self.TC = self.params['TC']
        self.TK = self.TC + 273.15
        self.alloy = self.dep + ''.join(self.inds)
        self.stamp = f'{self.phase}-{self.alloy}-{round(self.TC)}C'

        self.volume_db = get_and_log(self.params, 'volume_db',
                                     self.default_parameters['volume_db'],
                                     self.logger)
        self.V_partial = da.get_volume_data(self.volume_databases,
                                            self.volume_db,
                                            self.comps,
                                            self.logger)
        self.vacancy_db = get_and_log(self.params, 'vacancy_db',
                                      self.default_parameters['vacancy_db'],
                                      self.logger)
        self.GfV = da.get_vacancy_formation_energy(self.vacancy_databases,
                                                   self.vacancy_db,
                                                   self.phase,
                                                   self.comps[1:],
                                                   self.logger)
        self.thermo_params = None
        self.mob_params = None

        self.k_dislo = self.params.get('k_dislo', None)
        self.k_pores = self.params.get('k_pores', None)
        self.rho_dislo = self.params.get('rho_dislo', None)
        self.rho_pores = self.params.get('rho_pores', None)

        # Log some attributes needed to read simulation results
        self.logger.data(f'comps = {self.comps}')
        self.logger.data(f'V_partial = {self.V_partial}')

        # Initialize attributes
        self.G_fun = None
        self.G_funx = None
        self.G_funy = None
        self.MU_fun = None
        self.MU_funx = None
        self.MU_funy = None
        self.yVa_fun = None
        self.DT_fun = None
        self.DVa_fun = None
        self.L_fun = None

        # Define attributes
        self.logger.info("Processing thermodynamic and mobility data.")
        self.ideal_lattice = self.is_lattice_ideal()
        self.add_thermo_functions()
        self.add_mob_functions()
        self.logger.info("Processing lattice sink terms.")
        if self.ideal_lattice is False:
            self.process_sink_terms()

    def get_user_data(self):
        """
        Get data from user.

        * Molar volume databases
        * Vacancy formation energy databases
        * Register of thermodynamics and mobility databases
        * Default parameters

        """
        user_data = da.get_user_data(self.data_dir, self.logger)
        self.volume_databases = self.get_or_raise(user_data, 'molar_volume')
        self.vacancy_databases = self.get_or_raise(user_data,
                                                   'vacancy_formation_energy')
        thermo_register = self.get_or_raise(user_data, 'thermodynamics')
        mobility_register = self.get_or_raise(user_data, 'mobility')
        # Merge and make all keys lowercase
        database_register = thermo_register | mobility_register
        self.database_register = {k.lower(): val for k, val in
                                  database_register.items()}
        # Get numerical parameters from user data or from constants module
        self.default_parameters = {}
        factory = co.factory_default_parameters
        if 'default_parameters' in user_data:
            user_dct = user_data['default_parameters']
            for k in factory:
                if k in user_dct:
                    self.default_parameters[k] = user_dct[k]
                else:
                    self.default_parameters[k] = factory[k]
        else:
            self.default_parameters = factory

    def get_or_raise(self, dictionary, key):
        """Get value or raise exception if key not present."""
        try:
            val = dictionary[key]
        except KeyError:
            msg = f"Entry '{key}' not found in 'user_data.toml' file."
            raise ut.UserInputError(msg) from None
        return val

    def is_lattice_ideal(self):
        """
        Check whether lattice is ideal and print info messages.

        If the input files contains no value for either k_dislo, k_pores,
        rho_dislo or rho_pores, the lattice is ideal. If either parameter is
        given, the lattice is non-ideal.

        """
        if (self.k_dislo is None and self.k_pores is None
            and self.rho_dislo is None and self.rho_pores is None):
            msg = ("Ideal lattice: sink term set to maintain vacancy fraction "
                   "at equilibrium. No pore will form.")
            self.logger.info(msg)
            res = True
        else:
            msg = ("Non-ideal lattice: the Kirkendall effect may produce "
                   "non-equilibrium vacancy fractions and porosity.")
            self.logger.info(msg)
            res = False
        return res

    def add_thermo_functions(self):
        """Add thermodynamics functions."""
        self.add_thermo_model()
        self.Va_method = 'RK'
        # TIP Va_method is always 'RK' -> no longer useful currently. Could
        # either suppress or reimplement composition-fixed equilibrium ('cst')
        # vacancy fraction.
        self.yVa_fun = wrap.yVa_fun_from_params(self.comps,
                                                self.thermo_params,
                                                self.TK)
        msg = "Equilibrium vacancy fraction computed from G model."
        self.logger.info(msg)

        self.G_fun = self.G_funx
        self.MU_fun = self.MU_funx

    def add_thermo_model(self):
        """Add analytical models for G and MU."""
        db = self.get_or_raise(self.database_register, self.thermo_db)
        fpath = self.data_dir / db["file"]
        self.thermo_params = da.get_thermo_from_file(fpath, self.phase,
                                                     self.comps[1:],
                                                     self.TK,
                                                     self.logger)
        db = self.get_or_raise(self.database_register[self.thermo_db], 'file')
        self.thermo_params['Va'] = co.GV0
        for k in self.comps[1:]:
            solvent = 'Va' + k
            p0 = tfu.LkV_fun(self.GfV, k, self.TK)
            self.thermo_params[solvent] = [p0, 0]
        for i in range(1, len(self.comps) - 1):
            for j in range(i + 1, len(self.comps)):
                A = self.comps[i]
                B = self.comps[j]
                self.thermo_params['Va' + A + B] = [0, 0]
        self.G_funx = wrap.Gfun_from_params(self.comps[1:], self.thermo_params,
                                            self.TK)
        self.MU_funx = wrap.MUfun_from_params(self.comps[1:],
                                              self.thermo_params, self.TK)
        if self.ideal_lattice is True:
            self.MU_funy = self.extend_MU_funx()
        else:
            self.G_funy = wrap.Gfun_from_params(self.comps, self.thermo_params,
                                                self.TK)
            self.MU_funy = wrap.MUfun_from_params(self.comps,
                                                  self.thermo_params,
                                                  self.TK)
        msg = f"Thermodynamic functions built from parameter file ({db})."
        self.logger.info(msg)

    def extend_MU_funx(self):
        """
        Make function that computes MU from site fractions.

        Here the vacancy site fraction is at its equilibrium value (MU_Va = 0).

        """
        def fun(y):
            x = y[1:]/(1 - y[0])
            MU = self.MU_funx(x)
            res = np.vstack((np.zeros(MU.shape[1]), MU))
            return res
        return fun

    def add_mob_functions(self):
        """Add mobility functions."""
        self.add_mob_model()
        self.L_fun = tfu.make_Lfun(self.DT_fun, self.TK)
        self.DVa_fun = tfu.make_DVa_fun(self.DT_fun)

    def add_mob_model(self):
        """Add analytical model for DT."""
        db = self.get_or_raise(self.database_register, self.mob_db)
        fpath = self.data_dir / db["file"]
        self.mob_params = da.get_mob_from_file(fpath, self.comps[1:], self.TK,
                                               self.logger)
        db = self.get_or_raise(self.database_register[self.mob_db], 'file')
        self.DT_fun = wrap.DTfun_from_params(self.comps[1:], self.mob_params,
                                             self.TK)
        msg = f"Mobility functions built from parameter file ({db})."
        self.logger.info(msg)

    def process_sink_terms(self):
        """
        Process user input for non-ideal lattice parameters.

        Parameters of interest: k_dislo, k_pores, rho_dislo, rho_pores. Input
        can be a float or a string indicating a file name. If neither k_dislo
        or rho_dislo is given, k_dislo defaults to 0. Same for k_pores and
        rho_pores.

        """
        for name in ['dislo', 'pores']:
            k = getattr(self, f"k_{name}")
            rho = getattr(self, f"rho_{name}")
            if k is None and rho is None:
                setattr(self, k, 0)
                msg = f"No value found for k_{name}, defaults to 0."
                self.logger.warning(msg)
            if isinstance(k, str):
                fpath = self.work_dir / k
                setattr(self, f"k_{name}", np.genfromtxt(fpath))
            if isinstance(rho, str):
                fpath = self.work_dir / rho
                setattr(self, f"rho_{name}", np.genfromtxt(fpath))


def intro_msg():
    """Prepare intro message."""
    msg = '-'*35 + '  Noda  ' + '-'*36 + '\n'
    msg += '| A Python package for simulating diffusion in multicomponent alloys.'
    msg += ' '*9 + '|\n'
    msg += f'| Version {__version__}' + ' '*63 + '|\n'
    msg += '-'*79
    return msg
