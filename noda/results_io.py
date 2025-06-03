# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Load/save simulation results."""

import re
import json
import ast

import numpy as np


def print_res_to_file(ref, sdir, pstring, log, results):
    """
    Print simulation results to file.

    Note: this is obsolete, kept for backward compatibility.

    Parameters
    ----------
    ref : str
        Reference of the simulation, prepended to file name.
    sdir: Path or str
        Path to sub-directory where file will be created.
    pstring : str
        Set of simulation parameters.
    log : str
        Simulation log.
    results : dict
        Simulation results (see :func:`simu.Simulation.run` and
        :func:`solvers.solver`).

    Returns
    -------
    None.

    """
    res_lists = {n: {var: d[var].tolist() for var in d}
                 for n, d in results.items()}

    jdict = {'parameters': pstring.split('\n'),
             'log': log.split('\n'),
             'results': res_lists}

    fpath = sdir / f'{ref}-results.json'

    with open(fpath, 'w', encoding='utf-8') as file:
        json.dump(jdict, file, indent=1)


def read_res_from_file(ref, sdir):
    """
    Get simulation results from file.

    The file format is expected to conform to :func:`print_res_to_file`.

    Note: this is obsolete, kept for backward compatibility.

    Parameters
    ----------
    ref: str
        Reference at start of file name.
    sdir: Path or str
        Path to sub-directory where file is to be found.

    Returns
    -------
    pstring : str
        Set of simulation parameters.
    log : str
        Simulation log.
    results : dict
        Simulation results (see :func:`simu.Simulation.run` and
        :func:`solvers.solver`).

    """
    fpath = sdir / f'{ref}-results.json'

    with open(fpath, 'r', encoding='utf-8') as file:
        raw = file.read()

    jdict = json.loads(raw)

    pstring = '\n'.join(jdict['parameters'])
    log = '\n'.join(jdict['log'])
    res_lists = jdict['results']
    results = {int(n): {var: np.array(d[var]) for var in d}
               for n, d in res_lists.items()}

    return pstring, log, results


def parse_log(ref, sdir):
    """
    Get simulation results from log file.

    Parameters
    ----------
    ref: str
        Simulation reference.
    sdir: Path or str
        Path to sub-directory where file is to be found.

    Returns
    -------
    data : dict
        Simulation parameters (from input file + some logged in
        :class:`simu.System` and :class:`simu.Simulation`).
    results : dict
        Simulation results (see :func:`simu.Simulation.run` and
        :func:`solvers.solver`).

    """
    fpath = sdir / f'{ref}.nod'

    with open(fpath, 'r', encoding='utf-8') as file:
        raw = file.read()

    data = re.findall(r"DATA   : (\w+)\s*=\s*(.*)\n", raw)
    data = {x[0]: x[1] for x in data}
    for x in ['comps', 'V_partial', 'saved_steps', 'saved_times']:
        data[x] = ast.literal_eval(data[x])
    for x in ['nt', 'dt', 'th']:
        data[x] = float(data[x].split('#')[0].rstrip())
    try:  # back-compatibility: older versions did not have zstep saved
        data['zstep'] = float(data['zstep'].split('#')[0].rstrip())
    except KeyError:
        data['zstep'] = None

    res_match = re.findall(r"RESULTS: (.*)", raw)
    res_di = json.loads(res_match[0].replace("'", '"'))
    results = {int(n): {var: np.array(d[var]) for var in d}
               for n, d in res_di.items()}

    return data, results
