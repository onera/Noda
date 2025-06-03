# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Thermodynamics and mobility functions based on analytical models."""

import numpy as np

import noda.thermo_functions as tfu
from noda.utils import make_combinations
from noda.constants import R


def Gfun_from_params(comps, pdict, TK):
    """
    Make function to compute Gibbs free energy.

    Process parameters and call :func:`thermo_functions.G_model`.

    Parameters
    ----------
    comps : list of str
        System constituents.
    pdict : dict
        Thermodynamic parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    fun : function
        Function that takes composition array (shape (`n_inds`, `nz`)) as
        argument, and returns Gibbs free energy evaluated on this composition
        grid (shape (`nz`,)).

    """
    solvents = make_combinations(comps)

    G0 = [pdict[k] for k in comps]
    L = []
    for s in solvents['mix']:
        L.append(pdict[s][0])
        L.append(pdict[s][1])
    L = np.array(L)

    def fun(x):
        return tfu.G_model(x, G0, L, TK)
    return fun


def MUfun_from_params(comps, pdict, TK):
    """
    Make function that computes chemical potentials.

    Process parameters and call :func:`thermo_functions.MU_model`.
    Use a Redlich-Kister polynomial with binary interactions of order 0
    and 1 (see :func:`thermo_functions.MU_model`).

    Parameters
    ----------
    comps : list of str
        System constituents.
    pdict : dict of floats
        Thermodynamic parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    fun : function
        Function that takes composition array (shape (`n_inds`, `nz`)) as
        argument, and returns chemical potentials evaluated on this composition
        grid (shape (`n_comps`, `nz`)).

    """
    solvents = make_combinations(comps)
    parr = [pdict[k] for k in comps]
    for s in solvents['mix']:
        parr.append(pdict[s][0])
        parr.append(pdict[s][1])

    def fun(x):
        return tfu.MU_model(x, parr, TK)
    return fun


def yVa_fun_from_params(comps, pdict, TK):
    """
    Make function to compute equilibrium vacancy site fractions.

    Process parameters and call :func:`thermo_functions.yVa_model`.

    Parameters
    ----------
    comps : list of str
        System constituents, including Va.
    pdict : dict
        Thermodynamic parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    fun : function
        Function that takes atom fractions (shape (`n_inds`, `nz`)) as
        argument, and returns the equilibrium vacancy fraction evaluated on
        this composition array (shape (`nz`,)).

    """
    solvents = make_combinations(comps)
    G0 = [pdict[k] for k in solvents['unaries']]
    L = np.array([pdict[k] for k in solvents['mix']])

    def fun(x):
        return tfu.yVa_model(x, G0, L, TK)
    return fun


def yVa_fun_compfixed(GfV, TK):
    """
    Make function to compute equilibrium vacancy site fractions.

    Here, the equilibrium vacancy fraction depends on the temperature but not
    on the composition. The parameters provided in GfV normally correspond to
    the dependent constituent.

    Parameters
    ----------
    GfV : list
        Gibbs free energy of vacancy formation in a metal.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    fun : function
        Function that takes atom fractions (shape (`n_inds`, `nz`)) as
        argument, and returns the equilibrium vacancy fraction evaluated on
        this composition array (shape (`nz`,)).

    """
    def fun(x):
        Gf = GfV[0] - TK*GfV[1]
        val = np.exp(-Gf/(R*TK))
        return np.ones(x.shape[-1])*val
    return fun


def DTfun_from_params(comps, pdict, TK):
    """
    Generate DT function based on :func:`thermo_functions.lnDT_model`.

    This uses a Redlich-Kister polynomial with binary and ternary interactions
    of order 0.

    Parameters
    ----------
    comps : list of str
        System constituents.
    pdict : dict of floats
        Mobility parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    fun : function
        Function that takes composition array (shape (`n_inds`, `nz`)) as
        argument, and returns variable evaluated on this composition grid
        (shape (`n_comps`, `nz`)).

    """
    solvents = make_combinations(comps)['all']
    parr = [pdict[k][s] for k in comps for s in solvents]

    def fun(x):
        return np.exp(tfu.lnDT_model(x, parr, TK))
    return fun
