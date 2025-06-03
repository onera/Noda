# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Parse the input file that defines the system and various options."""

import re
import io
from pathlib import Path
import itertools as it

import pandas as pd

import noda.utils as ut


alias_db = {'k_dislo': ['kin_homo', 'km', 'nu_dislo'], # backward compatibility
            'k_pores': ['kin_pore', 'kp', 'nu_pores'], # backward compatibility
            'x_profile': 'profile_type',
            'grid': 'grid_type',
            'dep': 'dependent',
            'inds': 'independent',
            'nt_multiplier': 'nt_mult'}

extra_entries = ['nz', 'grid_type', 'nt_multiplier',
                 'xBC_left', 'xBC_right',
                 'JBC_left', 'JBC_right',
                 'xrange', 'independent',
                 'tim', 'mim',  # backward compatibility
                 ]

required_entries = ['dep', 'phase', 'thermo_db', 'mob_db', 'TC']

phase_alias_db = {'fcc': ['fcc_a1']}


def read_log_string(s):
    """
    Read the log created by a diffusion simulation.

    Parameters
    ----------
    s : str
        Simulation log.

    Returns
    -------
    di : dict
        Parameters and options defining the simulation.

    """
    sio = io.StringIO(s)
    df = pd.read_csv(sio, sep=r'\s+=\s+', engine='python', header=None,
                     comment='#', index_col=0)
    di = df.T.to_dict(orient='list')
    for var in di:
        di[var] = di[var][0]
    return di


def scrap_simu_attributes():
    """
    Generate list of valid attributes.

    Get list of all attributes of :class:`alloy_system.AlloySystem` and
    :class:`simu.Simulation`.

    """
    attributes = []
    for module in ['alloy_system', 'simu']:
        fpath = Path(__file__).parents[0] / f'{module}.py'
        with open(fpath, encoding='utf-8') as f:
            raw = f.read()
        attr = re.findall(r'self.(\w+) = (?:self.params|get_and_log)', raw)
        attr = [a for a in attr if "func" not in a]
        attributes += attr
    return list(set(attributes))


def find_parameter_name(ref_name, di):
    """
    Find key used for the parameter ref_name in the input dict.

    For some parameters in the user input file, there exists one prefered name
    and one or more aliases. This function recognizes the name used as an alias
    of the reference name as key of the input dictionary and returns it. This
    was implemented for backward compatibility reasons (aliases are former
    parameter names).

    Parameters
    ----------
    ref_name : str
        Reference name of the parameter.
    di : dict
        Input parameters.

    Raises
    ------
    Exception
        If input dict contains two or more keys that correspond to the same
        parameter.

    Returns
    -------
    par_name : str
        Key used for the parameter in the input dict.

    """
    aliases = alias_db[ref_name]
    if not isinstance(aliases, list):
        aliases = [aliases]

    possible_keys = aliases + [ref_name]

    for k1, k2 in it.combinations(possible_keys, r=2):
        if k1 in di and k2 in di:
            msg = f"Use either {k1} or {k2}, not both."
            raise ut.UserInputError(msg) from None

    used_aliases = list(di.keys() & aliases)
    if len(used_aliases) == 1:
        par_name = used_aliases[0]
    else:
        par_name = ref_name
    return par_name


def handle_aliases(di, logger):
    """
    Replace aliases by their prefered equivalent name.

    Use :func:`find_parameter_name` to find all possible aliases in input dict,
    and then replace them by the prefered parameter name.

    """
    for ref_name in alias_db:
        par_name = find_parameter_name(ref_name, di)
        if par_name != ref_name:
            di[ref_name] = di[par_name]
            del di[par_name]
            msg = (f"Parameter name '{par_name}' replaced by prefered name"
                   f" '{ref_name}'.")
            logger.data(msg)
    return di


def handle_phase_aliases(phase, logger):
    """
    Replace phase alias by the prefered phase name.

    Main purpose is backward compatibility.

    """
    for ref_name, aliases in phase_alias_db.items():
        if phase in aliases:
            msg = (f"Phase entry '{phase}' replaced by prefered name"
                   f" '{ref_name}'.")
            logger.info(msg)
            return ref_name
    return phase


def process_xrange(xr_raw_string):
    """Make dict of atom fraction ranges from input string."""
    xr_strings = xr_raw_string.split(',')
    xr_lists = [re.split(':|-', xr) for xr in xr_strings]
    try:
        xr_raw = {xr[0].strip(' '): [float(xr[1]), float(xr[2])]
                  for xr in xr_lists}
    except ValueError:
        msg = (f"Invalid entry '{xr_raw_string}' in input file. "
               "Required format for range of allowed atom fractions: "
               "'element1: min - max, element2: min - max, ...'.")
        raise ut.UserInputError(msg) from None
    res = {ut.format_element_symbol(k): v for k, v in xr_raw.items()}
    return res


def get_independent_constituents(inds_raw):
    """Make list of independent constituents from string."""
    res = inds_raw.split(',')
    res = [k.strip() for k in res]
    res = [ut.format_element_symbol(k) for k in res]
    return res


def preprocess_initial_profile(profile_string):
    """Separate initial profile input into keyword and option."""
    try:
        groups = profile_string.split(' ')
        words = [' '.join(x for x in groups if x.isalpha())]
        option = [float(x) for x in groups if not x.isalpha()]
        res = words + option
    except ValueError:
        res = profile_string
    return res


def check_profile_config(di, var):
    """Make sure input related to initial profile of var is consistent."""
    profile = di[f'{var}_profile']
    par_name = f"{var}_profile"
    original_input = ' '.join(str(x) for x in profile)
    msg = f"Found {par_name} = {original_input} in input. "
    if ut.isfilename(profile):
        if f'{var}_left' in di or f'{var}_right' in di:
            msg += f"{par_name} parameter is a filename. "
            msg += f"Cannot specify {var}_left or {var}_right."
            raise ut.UserInputError(msg) from None
    elif profile[0] == 'step' or profile[0] == 'smooth step':
        if f'{var}_left' not in di or f'{var}_right' not in di:
            msg += f"Parameters {var}_left and {var}_right are required."
            raise ut.UserInputError(msg) from None
        if len(profile) < 2:
            msg += "Step position is missing."
            raise ut.UserInputError(msg) from None
        if len(profile) > 2:
            msg += "Step type profiles take only one argument."
            raise ut.UserInputError(msg) from None
        if not isinstance(profile[1], (float, int)):
            msg += "Step argument should be a float (or int)."
            raise ut.UserInputError(msg) from None
    elif profile[0] == 'flat':
        if f'{var}_left' not in di:
            msg += f"Parameter {var}_left is required."
            raise ut.UserInputError(msg) from None
        if len(profile) > 1:
            msg += "Flat profile takes no argument."
            raise ut.UserInputError(msg) from None
    else:
        msg += "Invalid profile parameter."
        raise ut.UserInputError(msg) from None


def process_x_profile(x_string):
    """Make dict of atom fractions from initial condition input string."""
    side_strings = x_string.split(',')
    side_list = [re.split(':', x) for x in side_strings]
    try:
        res = {x[0].strip(' '): float(x[1]) for x in side_list}
    except ValueError:
        msg = (f"Invalid entry '{x_string}' in input file. "
               "Required format for initial atom fractions: "
               "'element1: value, element2: value, ...'.")
        raise ut.UserInputError(msg) from None
    res = {ut.format_element_symbol(k): v for k, v in res.items()}
    return res


def process_var_profile(var_string, side):
    """Make dict of var values from initial condition input string."""
    try:
        res = float(var_string)
        try:
            assert 0 <= res <= 1
        except AssertionError as exc:
            raise f"Value of {side} must be between 0 and 1." from exc
    except ValueError:
        res = var_string
    return res


def make_BC_type(di, side):
    """Guess BC type from input dict and make sure input is consistent."""
    if f'xBC_{side}' in di:
        res = 'Dirichlet'
        if f'JBC_{side}' in di:
            msg = f"{side} boundary condition: cannot fix both "
            msg += "composition and flux"
            raise ut.UserInputError(msg)
    elif f'JBC_{side}' in di:
        res = 'Neumann'
    else:
        res = None
    return res


def make_BC_dict(di, var):
    """Convert raw BC input into BC dict."""
    var_strings = di[var].split(',')
    var_list = [re.split(':', x) for x in var_strings]
    res = {}
    for x in var_list:
        constituent = x[0].strip(' ')
        constituent = ut.format_element_symbol(constituent)
        res[constituent] = x[1].strip(' ')
    return res


def check_num_out(di):
    """Make sure num_out is greater than 2 if present."""
    if 'num_out' in di:
        if di['num_out'] < 2:
            msg = ("'num_out' must be greater than 2. "
                   f"Found 'num_out = {di['num_out']}'")
            raise ut.UserInputError(msg) from None


def check_time_multiplier_config(di):
    """Make sure input related to time multipliers is consistent."""
    if 'nt_multiplier' in di:
        if 'dt_multiplier' in di:
            msg = "Use either dt_multiplier or nt_multiplier, not both."
            raise ut.UserInputError(msg) from None
        di['dt_multiplier'] = 1/di['nt_multiplier']


def check_sink_config(di):
    """Make sure input related to sink strength is consistent."""
    for name in ['dislo', 'pores']:
        if f'rho_{name}' in di and 'k_{name}' in di:
            msg = f"Cannot specify both k_{name} and rho_{name}."
            raise ut.UserInputError(msg) from None


def compact(s):
    """Remove empty lines and commented entries in input string."""
    lines = s.split('\n')
    res = '\n'.join([line for line in lines
                     if line != '' and not line.startswith('#')])
    return res


def input_string_to_dict(s):
    """Convert input string into a dictionary."""
    buffer = io.StringIO(s)
    df = pd.read_csv(buffer, sep=r'\s+=\s+', engine='python', header=None,
                     comment='#', index_col=0)
    entries = list(df.T.columns)
    for x in entries:
        if entries.count(x) > 1:
            msg = f"Parameter {x} found more than once in input file."
            raise ut.UserInputError(msg) from None
    di = df.T.to_dict(orient='list')
    for var in di:
        di[var] = di[var][0].rstrip()
    return di


def make_valid_entries():
    """List all recognized parameters for user input file."""
    valid_entries = scrap_simu_attributes()
    for k, v in alias_db.items():
        valid_entries.append(k)
        if isinstance(v, str):
            valid_entries.append(v)
        elif isinstance(v, list):
            valid_entries += v
    valid_entries += extra_entries
    valid_entries += required_entries
    return list(set(valid_entries))


def check_valid_entries(di):
    """Make sure all entries in input dict are valid."""
    valid_entries = make_valid_entries()
    for key in di:
        if key not in valid_entries:
            msg = f"Unknown entry '{key}' in input file."
            raise ut.UserInputError(msg) from None


def check_required_entries(di):
    """Make sure all required entries are present in dict."""
    for x in required_entries:
        if x not in di:
            msg = f"Missing required parameter '{x}' in input file."
            raise ut.UserInputError(msg) from None


def read_parameter_string(s, logger):
    """
    Parse input file to get set of parameters.

    In the file, parameter names and values must be separated by '='.
    Commented lines (#) and blank lines are ignored.
    Blank spaces around '=' delimiters are ignored.
    Blank spaces in future dictionaries (x_left, x_right) are ignored.
    See :class:`simu.System` for role of parameters.

    Parameters
    ----------
    s: str
        Input file given as str.

    Returns
    -------
    compacted_string : str
        Input string without unnecessary lines.
    di: dict
        Set of parameters that define a system, simulation conditions and
        various options.

    """
    # Pre-processing
    compacted_string = compact(s)
    di = input_string_to_dict(s)
    di = handle_aliases(di, logger)
    check_valid_entries(di)
    check_required_entries(di)
    di['phase'] = handle_phase_aliases(di['phase'], logger)

    # System constituents
    try:
        inds = get_independent_constituents(di['inds'])
    except KeyError:
        xrange = process_xrange(di['xrange'])
        inds = list(xrange.keys())
    di['inds'] = sorted(inds)
    di['dep'] = ut.format_element_symbol(di['dep'])

    # Initial profiles
    init_profile_variables = ['x', 'yVa', 'fp']
    for var in init_profile_variables:
        key = var + '_profile'
        if key in di:
            di[key] = preprocess_initial_profile(di[key])
            check_profile_config(di, var)

    for side in ['x_left', 'x_right']:
        if side in di:
            di[side] = process_x_profile(di[side])

    for var in [k for k in init_profile_variables if k != 'x']:
        for side in [var + '_left', var + '_right']:
            if side in di:
                di[side] = process_var_profile(di[side], side)

    # Boundary conditions
    for side in ['left', 'right']:
        di[f'BC_{side}_type'] = make_BC_type(di, side)

    BC_dicts = ['xBC_left', 'xBC_right', 'JBC_left', 'JBC_right']
    for var in BC_dicts:
        if var in di:
            di[var] = make_BC_dict(di, var)

    # Convert numerical values to float if present
    for var in di:
        try:
            di[var] = float(di[var])
        except (TypeError, ValueError):
            pass

    # Convert numerical values to int if present
    for var in ['nz', 'num_out']:
        try:
            di[var] = int(di[var])
        except KeyError:
            pass

    # Other sanity checks
    check_num_out(di)
    check_time_multiplier_config(di)
    check_sink_config(di)

    return compacted_string, di


if __name__ == "__main__":
    accepted = make_valid_entries()
    print("\n".join(accepted))
