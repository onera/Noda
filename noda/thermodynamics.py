# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Handle thermodynamic properties."""

import numpy as np

import noda.thermo_functions as tfu
import noda.constants as co
import noda.utils as ut
import noda.data_io as da


class Thermodynamics:
    """
    Provide methods to compute thermodynamic properties.

    Methods that compute the Gibbs free energy and chemical potentials as a
    function of atom fraction array (G_funx and MU_funx) assume vacancies are
    at equilibrium.

    Methods that compute the Gibbs free energy and chemical potentials as a
    function of site fraction array (G_funy and MU_funy) are initialized based
    on the same assumption, but the class provides a method to update these
    methods in the case of a non-ideal lattice (set_nonideal_lattice).

    Attributes
    ----------
    comps : list of str
        System components.
    TK : float
        Temperature in Kelvin.
    GfV : dict of list
        Vacancy formation energy in pure metals, given as [enthalpy, entropy]
        in J/mol and J/mol/K. See :func:`data_io.get_vacancy_formation_energy`.
    params : dict of floats
        Thermodynamic parameters arranged as follows:

        | ``A: G_A for A in endmembers``
        | ``AB: [L0, L1] for AB in binary subsystems``

    ideal_lattice : bool
        Whether lattice is ideal, in the sense that the vacancy fraction is
        maintained at equilibrium.

    Methods
    -------
    G_funx(x) :
        Returns Gibbs free energy at composition given by atom fraction array.
    MU_funx(x) :
        Returns chemical potentials at composition given by atom fraction
        array.
    yVa_fun(x) :
        Returns equilibrium vacancy fraction at composition given by atom
        fraction array.
    G_fun(x) :
        Alias of G_fun.x
    MU_fun(x):
        Alias of MU_funx.
    G_funy(y):
        Returns Gibbs free energy at composition given by site fraction array.
    MU_funy :
        Returns chemical potentials at composition given by site fraction
        array.

    """
    def __init__(self, params, comps, phase, TK, vacancy_databases,
                 vacancy_db, logger):
        """
        Class constructor.

        Parameters
        ----------
        params : dict of floats
            Thermodynamic parameters arranged as follows:

            | ``A: G_A for A in endmembers``
            | ``AB: [L0, L1] for AB in binary subsystems``

        comps : list of str
            System components.
        phase : str
            Phase name.
        TK : float
            Temperature in Kelvin.
        vacancy_databases : dict
            Register of vacancy formation energy databases.
        vacancy_db : str
            Name of database with vacancy formation energy in pure metals.
        logger : :class:`log_utils.CustomLogger`
            Logger.

        """
        self.comps = comps
        self.TK = TK
        self.GfV = da.get_vacancy_formation_energy(vacancy_databases,
                                                   vacancy_db,
                                                   phase,
                                                   comps[1:],
                                                   logger)
        self.params = self.process_parameters(params)
        self.G_funx = make_G_fun(comps[1:], self.params, TK)
        self.MU_funx = make_MU_fun(comps[1:], self.params, TK)
        self.yVa_fun = make_yVa_fun(comps, self.params, TK)
        self.G_fun = self.G_funx
        self.MU_fun = self.MU_funx
        self.ideal_lattice = True
        self.G_funy = self.extend_G_funx()
        self.MU_funy = self.extend_MU_funx()

    def process_parameters(self, params):
        """Add vacancy-related interactions to thermodynamic parameters."""
        params['Va'] = co.GV0
        for k in self.comps[1:]:
            solvent = 'Va' + k
            p0 = tfu.LkV_fun(self.GfV, k, self.TK)
            params[solvent] = [p0, 0]
        for i in range(1, len(self.comps) - 1):
            for j in range(i + 1, len(self.comps)):
                A = self.comps[i]
                B = self.comps[j]
                params['Va' + A + B] = [0, 0]
        return params

    def extend_G_funx(self):
        """
        Make function that computes G from site fractions.

        Here the vacancy site fraction is at its equilibrium value.

        """
        def fun(y):
            x = y[1:]/(1 - y[0])
            return self.G_funx(x)
        return fun

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

    def set_nonideal_lattice(self):
        """
        Update G_funy and MU_funy functions in cas of non-ideal lattice.

        Update functions that compute Gibbs free energy and chemical potentials
        from site fraction array when lattice is non-ideal (vacancy site
        fraction not maintained at equilibrium).

        """
        self.ideal_lattice = False
        self.G_funy = make_G_fun(self.comps, self.params, self.TK)
        self.MU_funy = make_MU_fun(self.comps, self.params, self.TK)


def make_G_fun(comps, pdict, TK):
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
        argument, and returns Gibbs free energy evaluated on this
        composition grid (shape (`nz`,)).

    """
    solvents = ut.make_combinations(comps)

    G0 = [pdict[k] for k in comps]
    L = []
    for s in solvents['mix']:
        L.append(pdict[s][0])
        L.append(pdict[s][1])
    L = np.array(L)

    def fun(x):
        return tfu.G_model(x, G0, L, TK)
    return fun


def make_MU_fun(comps, pdict, TK):
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
        argument, and returns chemical potentials evaluated on this
        composition grid (shape (`n_comps`, `nz`)).

    """
    solvents = ut.make_combinations(comps)
    parr = [pdict[k] for k in comps]
    for s in solvents['mix']:
        parr.append(pdict[s][0])
        parr.append(pdict[s][1])

    def fun(x):
        return tfu.MU_model(x, parr, TK)
    return fun

def make_yVa_fun(comps, pdict, TK):
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
    solvents = ut.make_combinations(comps)
    G0 = [pdict[k] for k in solvents['unaries']]
    L = np.array([pdict[k] for k in solvents['mix']])

    def fun(x):
        return tfu.yVa_model(x, G0, L, TK)
    return fun
