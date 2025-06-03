# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Load thermodynamics and mobility data."""

import json
import tomllib

import pandas as pd

import noda.thermo_functions as tfu
import noda.utils as ut
import noda.constants as co
from noda.constants import factory_default_parameters as factory
from noda.paths import pkg_data_dir


def get_user_data(data_dir, logger):
    """
    Get user data from toml file.
    
    If toml file is not found in data folder, 

    Parameters
    ----------
    data_dir : pathlib.Path
        Path of data folder.
    
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Returns
    -------
    res : dict
        Content of file.

    """
    fpath = data_dir / 'user_data.toml'
    if fpath.exists():
        with open(fpath, 'rb') as file:
            res = tomllib.load(file)
    else:
        with open(pkg_data_dir / 'user_data.toml', 'rb') as file:
            res = tomllib.load(file)
        msg = (f"No 'user_data.toml' file found in {data_dir}. Using "
               "'user_data.toml' file from package installation directory "
               f"instead ({pkg_data_dir}).")
        logger.warning(msg)
    return res


def get_volume_data(volume_databases, volume_db, comps, logger):
    """
    Get partial molar volumes from specified database.

    Parameters
    ----------
    volume_databases : dict
        Available molar volume data.
    volume_db : str
        Name of partial molar volume database.
    comps : list of str
        System constituents.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    Exception
        If database is not present in databases dict.
    Exception
        If database entry is formatted incorrectly.

    Returns
    -------
    res : dict
        ``{k: V_k for k in comps}``

    """
    try:
        di = volume_databases[volume_db]
    except KeyError:
        msg = (f"Unknown molar volume database '{volume_db}'. "
               "Please check the 'molar_volume' table in your "
               "'user_data.toml' file")
        raise ut.UserInputError(msg) from None
    res = {}
    sorted_keys = comps + ['pore']
    for k in sorted_keys:
        try:
            res[k] = di[k]
        except KeyError:
            if 'default' in di:
                res[k] = di['default']
                msg = (f"No entry for {k} in molar volume database "
                       f"'{volume_db}'. Using default entry in the database.")
                if k in ('Va', 'pore'):
                    logger.data(msg)
                else:
                    logger.info(msg)
            else:
                res[k] = factory['partial_molar_volume']
                msg = (f"Molar volume database '{volume_db}' contains no data "
                       f"for {k}, and no default value. Using system-wide "
                       f"default value (Vm = {res[k]} m3/mol) instead.")
                if k in ('Va', 'pore'):
                    logger.data(msg)
                else:
                    logger.info(msg)
        if k in ('Va', 'pore'):
            if not (isinstance(res[k], float) or res[k] == 'local'):
                msg = (f"Invalid entry for species {k} in molar volume "
                       f"database '{volume_db}' (found '{res[k]}', should be "
                       "a float or string 'local').")
                raise ut.UserInputError(msg) from None
        else:
            if not isinstance(res[k], float):
                msg = (f"Invalid entry for species {k} in database "
                       f"'{volume_db}' (found '{res[k]}', should be a float).")
                raise ut.UserInputError(msg) from None

    return res


def get_vacancy_formation_energy(vacancy_databases, vacancy_db, phase, comps,
                                 logger):
    """
    Get vacancy formation energy in pure metals.

    Parameters
    ----------
    volume_databases : dict
        Available vacancy formation energy data.
    vacancy_db : str
        Name of database with vacancy formation energy in pure metals.
    phase : str
        Name of metal phase.
    comps : list of str
        System constituents.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    Exception
        If database is not included in databases dict.

    Returns
    -------
    res : dict
        ``{k: [enthalpy, entropy] for k in comps}``

    """
    try:
        di = vacancy_databases[vacancy_db]
    except KeyError:
        msg = (f"Unknown vacancy formation energy database '{vacancy_db}'. "
               "Please check the 'vacancy_formation_energy' table in your "
               "'user_data.toml' file")
        raise ut.UserInputError(msg) from None

    di_phase = {}
    for k in comps:
        try:
            di_phase[k] = di[f"{phase}-{k}"]
        except KeyError:
            if 'default' in di:
                di_phase[k] = di['default']
                msg = (f"No entry for {phase}-{k} in vacancy formation energy "
                       f"database '{vacancy_db}'. Using default entry in the "
                       "database.")
                logger.info(msg)
            else:
                di_phase[k] = factory['vacancy_formation_energy']
                msg = (f"Vacancy formation energy database '{vacancy_db}' "
                       f"contains no data for {k}, and no default value. "
                       f"Using system-wide default value ({di_phase[k]}) "
                       "instead.")
                logger.info(msg)

    res = {k: [v*co.EV*co.NA for v in di_phase[k]] for k in comps}
    return res


# =============================================================================
# Begin code from OPTIMOB
# =============================================================================

# =============================================================================
# Thermodynamic parameters

def get_thermo_from_file(fpath, phase, comps, TK, logger):
    """
    Get parameters needed to calculate Gibbs free energy from spreadsheet file.

    Parameters
    ----------
    fpath : pathlib.Path
        Path of file with thermodynamic database.
    phase : str
        Name of metal phase.
    comps : list of str
        System constituents.
    TK : float
        Temperature in Kelvin.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    Exception
        If file is not found.

    Returns
    -------
    p : dict of floats
        Thermodynamic parameters arranged as follows:

        | ``A: G_A for A in endmembers``
        | ``AB: [L0, L1] for AB in binary subsystems``
        | ``ABC: [L0, L1] for ABC in ternary subsystems``

    """
    if not fpath.exists():
        fname = fpath.name
        data_dir = fpath.parents[0]
        pkg_fpath = pkg_data_dir / fname
        if pkg_fpath.exists():
            fpath = pkg_fpath
            msg = (f"Thermodynamic database file '{fname}' not found in "
                   f"{data_dir}. Using file from package installation "
                   f"directory instead ({pkg_data_dir}).")
            logger.warning(msg)
        else:
            msg = (f"Thermodynamic database file '{fname}' not found in "
                   f"{data_dir} or in package installation directory. \n"
                   "Please provide a thermodynamic database file.")
            raise ut.UserInputError(msg) from None

    G0_para = get_G0_parameters(fpath)
    G0 = {k: tfu.G0_fun(G0_para[k], TK, phase) for k in comps}
    solvents = ut.make_combinations(comps)['mix']
    interactions = get_thermo_interaction_parameters(fpath, solvents, logger)
    L_para = interactions['L']
    L = make_L_isotherm(L_para, TK)
    p = {**G0, **L}
    return p


def get_G0_parameters(fpath):
    """
    Get parameters needed to compute Gibbs free energy of endmembers.

    Data retrieved from spreadsheet file. Requires an external dependency:

    ======  ============
    format  package name
    ======  ============
    xls     xlrd
    xlsx    openpyxl
    ods     odfpy
    ======  ============

    Data is in the form of G - H_SER. See Dinsdale 1991 [#Dinsdale_1991]_.

    Parameters
    ----------
    fpath : pathlib.Path
        Path of file with thermodynamic database.

    Returns
    -------
    dict
        Parameters, see in file.
    """
    df = pd.read_excel(fpath, sheet_name='Pure elements',
                       skiprows=2, index_col=0)
    di = df.to_dict()
    return {ut.format_element_symbol(k): v for k, v in di.items()}


def get_thermo_interaction_parameters(fpath, solvents, logger):
    """
    Get thermodynamic interaction parameters from spreadsheet file.

    Requires an external dependency:

    ======  ============
    format  package name
    ======  ============
    xls     xlrd
    xlsx    openpyxl
    ods     odfpy
    ======  ============

    The parameters belong to the following categories (variables) depending
    on the quantity they are related to:

    * L : Gibbs free energy
    * Tc : critical temperature
    * beta : magnetism

    The parameters correspond to:

    * binary interactions (orders 0 and 1)
    * ternary interactions (order 0, with `L1 = 0` for compatibility)

    Parameters are given as:

    * order 0: A and B in `L0 = A + B*T`
    * order 1: C and D in `L1 = C + D*T`

    If a variable is not included in the input file, all parameters are set to
    0 for this variable.

    Parameters
    ----------
    fpath : pathlib.Path
        Path to thermodynamic database file.
    solvents : list of str
        Binary and ternary subsystems, concatenated to strings.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Returns
    -------
    di: dict
        Thermodynamic interaction parameters,

        ``{var: subdi for var in ['L', 'Tc', 'beta']}``

        where subdi is result of :func:`process_interaction_parameters`.
    """
    df = pd.read_excel(fpath, sheet_name='Interactions', skiprows=2)
    di = {}
    for k in ['L', 'Tc', 'beta']:
        sub_df = df.loc[df['variable'] == k]
        if sub_df.empty:
            di[k] = {k: 0 for k in solvents}
        else:
            sub_df = sub_df.set_index('solvent')
            di[k] = process_interaction_parameters(sub_df, fpath, solvents,
                                                   logger)
    return di


def process_interaction_parameters(df, fpath, solvents, logger):
    """
    Process dataframe with thermodynamic interaction parameters.

    Apply sanitary checks and convert from dataframe to dict. For each
    subsystem in solvents, operation depends on number of matching keys in the
    dataframe:

    * if 0, set all interaction parameters to 0
    * if 1, get interaction parameters
    * if more than 1, raise exception.

    Parameters
    ----------
    df : dataframe
        Thermodynamic interaction parameters for one variable. Columns:

        * variable : either of L, Tc or beta (see
          :func:`get_thermo_interaction_parameters`).
        * solvent : constituents of subsystem concatenated to string.
        * A, B, C, D : interaction parameters, with

            * order 0 = A + B*T
            * order 1 = C + D*T

    fpath : pathlib.Path
        Path of file with thermodynamic database.
    solvents : list of str
        Binary and ternary subsystems, concatenated to strings.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    Exception
        If several equivalent subsystems (ie permutations of the same
        subsystem) are present in the database.

    Returns
    -------
    di_reduced : dict of dicts
        Thermodynamic interaction parameters,

        ``{k: {letter: val for letter in 'ABCD'} for k in solvents}``.

    """
    variable = df['variable'].values[0]

    for k in solvents:
        possible = ut.make_permutations_samesize(k)
        kfile_list = [s for s in df.index if s in possible]

        if len(kfile_list) == 0:
            df.loc[k] = [variable] + [0 for letter in 'ABCD']
            # Interactions parameters for Tc and beta are not always available,
            # and are not used in Noda -> do not issue any warning for these
            # variables.
            if variable not in ['Tc', 'beta']:
                msg = f"{k} interaction parameters for variable '{variable}' "
                msg += "are missing. Using 0 as default."
                logger.warning(msg)

        elif len(kfile_list) == 1:
            # In binary subsystems, if the 2 elements are reversed, the order 1
            # parameters must be changed because
            # G = ... + x_i*x_j*(L0_ij + L1_ij*(x_i - x_j)) + ...
            # with L1_ij = C + D*T
            # In ternary subsystems, order 1 parameters are 0 for the moment
            # Introducing order 1 parameters would be more complex -> order of
            # elements does not matter

            kfile = kfile_list[0]
            if k != kfile:
                df = df.rename(index={kfile: k})
                if len(possible) == 2:  # binary subsystems
                    df.loc[k, 'C'] *= -1
                    df.loc[k, 'D'] *= -1

        else:
            msg = (f'Input file: "{fpath.name}"\n'
                   f"{k} interaction parameters for {variable}\n"
                   f"Several equivalent solvents given: {kfile_list}")
            raise ut.UserInputError(msg) from None
    di = df.to_dict(orient='index')
    di_reduced = {k: di[k] for k in solvents}
    return di_reduced


def make_L_isotherm(L, T):
    """
    Evaluate interaction parameters at given temperature.

    Parameters
    ----------
    L : dict
        Interaction parameters,

        ``{k: {letter: val for letter in 'ABCD'} for k in solvents}``.
    T : float or int
        Temperature in Kelvin.

    Returns
    -------
    res : dict of lists
        Interaction parameters, ``{k: [L0, L1] for k in solvents}``.

    """
    res = {k: [L[k]['A'] + L[k]['B']*T,
               L[k]['C'] + L[k]['D']*T]
           for k in L}
    return res


# =============================================================================
# Mobility parameters

def get_mob_from_file(fpath, comps, TK, logger):
    """
    Get mobility parameters from input file.

    File formats currently supported:

    * xls, xlsx or ods: database from literature
    * json: database from OPTIMOB.

    Parameters
    ----------
    fpath : pathlib.Path
        Path of file with mobility database.
    comps : list of str
        System constituents.
    TK : float
        Temperature in Kelvin.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    Raises
    ------
    Exception
        If file is not found or file format not accepted.

    Returns
    -------
    p : dict of dicts
        ``{i: subdict for i in comps}``

        subdict: ``{j: val for j in subsystems}``.

    """
    if not fpath.exists():
        fname = fpath.name
        data_dir = fpath.parents[0]
        pkg_fpath = pkg_data_dir / fname
        if pkg_fpath.exists():
            fpath = pkg_fpath
            msg = (f"Mobility database file '{fname}' not found in "
                   f"{data_dir}. Using file from package installation "
                   f"directory instead ({pkg_data_dir}).")
            logger.warning(msg)
        else:
            msg = (f"Mobility database file '{fname}' not found in "
                   f"{data_dir} or in package installation directory. \n"
                   "Please provide a mobility database file.")
            raise ut.UserInputError(msg) from None

    ext = fpath.name.split('.')[-1]
    if ext.startswith('xls') or ext == 'ods':
        p = get_mob_from_spreadsheet(fpath, comps, TK)
    elif ext == 'json':
        p = get_mob_from_json(fpath, comps, TK)
    else:
        msg = f'Input file format (.{ext}) not accepted.'
        raise ut.UserInputError(msg) from None

    return p


def get_mob_from_spreadsheet(fpath, comps, TK):
    """
    Get mobility parameters from spreadsheet file.

    Requires an external dependency:

    ======  ============
    format  package name
    ======  ============
    xls     xlrd
    xlsx    openpyxl
    ods     odfpy
    ======  ============

    Parameters
    ----------
    fpath : pathlib.Path
        Path of file with mobility database.
    comps : list of str
        System constituents.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    p : dict of dicts
        ``{i: subdict for i in comps}``

        subdict: ``{j: val for j in subsystems}``

    """
    di = pd.read_excel(fpath, sheet_name=None, header=3)
    df = di[list(di)[0]]
    df.solute = df.solute.apply(ut.format_element_symbol)

    solvents = ut.make_combinations(comps)['all']

    p = {}
    for i in comps:
        p[i] = {}
        for j in solvents:
            redf = get_reduced_df(df, j, i)
            A = redf['A'].values[0]
            B = redf['B'].values[0]
            p[i][j] = A + B*TK

    return p


def get_mob_from_json(fpath, comps, TK):
    """
    Get mobility parameters from json file.

    Parameters
    ----------
    fpath : pathlib.Path
        Path of file with mobility database.
    comps : list of str
        System constituents.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    p : dict of dicts
        ``{i: subdict for i in comps}``

        subdict: ``{j: val for j in subsystems}``

    """
    # pylint: disable=too-many-locals
    with open(fpath, 'r', encoding='utf-8') as f:
        di_full = json.load(f)

    # Check that temperature is correct
    TC_file_list = []
    for k in di_full:
        if 'exp ' in k:
            di = di_full[k]
            TC_file_list.append(di_full[k]['TC'])
    TC_file = TC_file_list[0]
    assert all(TC == TC_file for TC in TC_file_list)
    TC = TK - 273.15
    if TC != round(TC_file, 0):
        msg = f'Optimized parameters are valid at {TC_file} *C, '
        msg += f'not compatible with simulation at {TC} *C.'
        raise ut.UserInputError(msg) from None

    # Get nested dicts with parameters
    di = di_full["popt"]

    # Keep values for required solute-solvent combinations
    solvents = ut.make_combinations(comps)['all']
    p = {}
    for i in comps:
        p[i] = {}
        for j in solvents:
            p[i][j] = di[i][j]['value']

    return p


def get_reduced_df(df, solvent, solute):
    """
    Filter dataframe to keep mobility data for solute in solvent.

    Parameters
    ----------
    df : dataframe
        Parameters to compute mobility of solutes in solvents.
    solvent : str
        Solvent of interest (constituents concatenated to string).
    solute : str
        Solute of interest.

    Raises
    ------
    Exception
        If several equivalent solvents (ie permutations of the same solvent)
        are present in the df.

    Returns
    -------
    res : dataframe
        Reduced dataframe.

    """
    possible_solvents = ut.make_permutations_samesize(solvent)
    possible_lower_key = [k.lower() for k in possible_solvents]

    res = df[(df.solute == solute)
             & (df.solvent.apply(str.lower).isin(possible_lower_key))]

    if len(res.index) == 0:
        msg = f"Mobility parameters for {solute} in {solvent} are missing. "
        raise ut.UserInputError(msg) from None

    if len(res.index) > 1:
        msg = (f"Mobility parameters for {solute} in {solvent}. "
               f"Several equivalent solvents given: {res.solvent.values}")
        raise ut.UserInputError(msg) from None

    return res

# =============================================================================
# End code from OPTIMOB
# =============================================================================
