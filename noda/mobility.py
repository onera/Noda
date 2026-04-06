# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Handle mobility properties."""

import numpy as np

import noda.thermo_functions as tfu
from noda.utils import make_combinations


class Mobility:
    """
    Methods to compute mobility-related quantities.

    Attributes
    ----------
    params : dict of dicts
        Mobility parameters arranged as follows:

        ``{i: subdict for i in comps}``

        subdict: ``{j: val for j in subsystems}``.
    comps : list of str
            Atom components.

    Methods
    -------
    DT_funx(x) :
        Tracer diffusion coefficients, see :func:`get_DT_funx`.
    L_fun(x):
        Onsager coefficients, see :func:`thermo_functions.make_Lfun`.
    DVa_fun(y, x)
        Vacancy diffusion coefficient, see
        :func:`thermo_functions.make_DVa_fun`.

    """
    def __init__(self, params, comps, TK):
        """
        Class constructor.

        Parameters
        ----------
        params : dict of dicts
            Mobility parameters arranged as follows:

            ``{i: subdict for i in comps}``

            subdict: ``{j: val for j in subsystems}``.

        comps : list of str
            Atom components.
        TK : float
            Temperature in Kelvin.

        """
        self.params = params
        self.comps = comps
        self.DT_funx = get_DT_funx(comps, params, TK)
        self.L_fun = tfu.make_Lfun(self.DT_funx, TK)
        self.DVa_fun = tfu.make_DVa_fun(self.DT_funx)
    
    def DT_fun(self, x):
        """
        Calculate tracer diffusion coefficients of atom components.

        Parameters
        ----------
        x : np.array (shape (`n_inds`, `nz`))
            Atom fractions.

        Returns
        -------
        dict
            Tracer diffusion coefficients, shape (`nz`)
            ``{k: DT_k for k in components}``.

        """
        return {k: self.DT_funx(x)[i] for i, k in enumerate(self.comps)}


def get_DT_funx(comps, pdict, TK):
    """
    Generate DT function based on :func:`thermo_functions.lnDT_model`.

    This uses a Redlich-Kister polynomial with binary and ternary
    interactions of order 0.

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
        Function that takes ayom fraction array (shape (`n_inds`, `nz`)) as
        argument, and returns tracer diffusion coefficients evaluated on this
        composition grid (shape (`n_comps`, `nz`)).

    """
    solvents = make_combinations(comps)['all']
    parr = [pdict[k][s] for k in comps for s in solvents]

    def fun(x):
        return np.exp(tfu.lnDT_model(x, parr, TK))
    return fun
