# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Functions related to thermodynamics and mobility models."""

import numpy as np
import numdifftools as nd

from noda import utils as ut
from noda.constants import R, GV0


# =============================================================================
# Gibbs free energy in pure metal

def Gmag_fun(T, Tcritical, beta, phase):
    """
    Compute magnetic part of Gibbs free energy.

    See Dinsdale 1991 [#Dinsdale_1991]_.
    """
    tau = T/Tcritical

    p = 0.28
    if phase.upper() == 'BCC_A2':
        p = 0.4

    D = 518/1125 + 11692/15975 * (1/p - 1)

    if tau <= 1:
        g = 1 - (79/(140*p*tau)
                 + 474/497*(1/p - 1)*(tau**3/6 + tau**9/135 + tau**15/600)
                 )/D
    else:
        g = -(tau**(-5)/10 + tau**(-15)/315 + tau**(-25)/1500)/D

    return R*T*np.log(beta + 1)*g


def G0_fun(Gpara, T, phase):
    """
    Compute Gibbs free energy in pure metal.

    Data in the form of G - H_SER. See Dinsdale 1991 [#Dinsdale_1991]_.
    """
    a = Gpara['-']
    b = Gpara['T']
    c = Gpara['T*ln(T)']
    d2 = Gpara['T^2']
    d3 = Gpara['T^3']
    e = Gpara['1/T']
    Tc = Gpara['Tcritical']
    beta = Gpara['beta']

    Gmain = a + b*T + c*T*np.log(T) + d2*T**2 + d3*T**3 + e/T

    if np.isnan(Tc) or np.isnan(beta):
        Gmag = 0
    else:
        Gmag = Gmag_fun(T, Tc, beta, phase)

    return Gmain + Gmag


# =============================================================================
# Gibbs free energy in alloy

def G_model(x, G0, L, T):
    r"""
    Compute molar Gibbs free energy at given composition.

    Use excess term with binary interactions of order 0 and 1, and ternary
    interactions of order 0.

    Interactions are represented with a Redlich-Kister polynomial:

    .. math::

       &G_{m} = \sum_{i = 1}^{n}{x_{i}G_{i}}
                + RT\sum_{i = 1}^{n}{x_{i}\ln x_{i}}
                + ^{xs}G_{m}

       &_{}^{xs}G_{m} = \sum_{i = 1}^{n - 1}{
                          \sum_{j = i + 1\ }^{n}{
                        x_{i}x_{j}\left(^0L_{ij} + ^1L_{ij}(x_i - x_j)\right)}}
       + \sum_{i = 1}^{n - 2}{
           \sum_{j = i + 1\ }^{n - 1}{
               \sum_{k = j + 1\ }^{n}{x_{i}x_{j}x_{k}L_{ijk}}}}

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    G0 : list of floats
        Gibbs free energy of the endmembers, len `ncomps`.
    L : 1D array
        Interaction parameters, given as a flat array.

        | Number of terms in a n-component system:
        | binary  : 2*C(n, 2) = n*(n - 1)
        | ternary : C(n, 3) = n*(n - 1)*(n - 2)/6
        |
        | In a 2-comp. system, :math:`L = [L_{12}]`
        | In a 3-comp. system, :math:`L = [L_{12}, L_{13}, L_{23}, L_{123}]`
        | where binary terms are couples :math:`L_{ij} = (^0L_{ij}, ^1L_{ij})`

    T : float or int
        Temperature in Kelvin.

    Returns
    -------
    1D array
        Gibbs free energy.

    """
    n = len(G0)
    xn = 1 - sum(x)
    x_ = np.vstack((x, xn))

    L_ = L.reshape(-1, 2)
    n_bin = ut.count_combinations(n)['binaries']
    L_bin = L_[:n_bin]
    L_ter = L_[n_bin:]

    mech = sum(G0[i]*x_[i] for i in range(n))
    smix = R*T*sum(x_[i]*np.log(x_[i]) for i in range(n))

    binary = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            idx = ut.index_binary_combination(i, j, n)
            binary += x_[i]*x_[j]*(L_bin[idx, 0] +
                                   L_bin[idx, 1]*(x_[i] - x_[j]))

    ternary = 0
    for i in range(n - 2):
        for j in range(i + 1, n - 1):
            for k in range(j + 1, n):
                idx = ut.index_ternary_combination(i, j, k, n)
                ternary += x_[i]*x_[j]*x_[k]*L_ter[idx, 0]

    return mech + smix + binary + ternary


# =============================================================================
# Chemical potentials

def partial_derivative(fun, x, i, args=(), dx=1e-10):
    """
    Compute the derivative of fun with respect to i-th row of x.

    Adapted from https://stackoverflow.com/a/20708578.

    """
    X = list(x)

    def wraps(y):
        X[i] = y
        return fun(np.array(X), *args)
    return nd.Derivative(wraps, step=dx)(x[i])


def MU_model(x, p, TK):
    r"""
    Compute chemical potentials.

    .. math:: \mu_k = G + \frac{\partial G}{\partial x_k}
                      - \sum_{i=0}^n {x_i \frac{\partial G}{\partial x_i}}

    The G model uses excess term with binary interactions of order 0 and 1 and
    ternary interactions of order 0 (see :func:`G_model`).

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    p : list of floats
        Thermodynamic parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    MU : 2D array
        Chemical potentials, shape (`ncomps`, `nz`).

    """
    nc = x.shape[0] + 1
    nz = x.shape[-1]

    xn = 1 - sum(x)
    x_ = np.vstack((x, xn))

    G0 = p[:nc]
    L = np.array(p[nc:])

    G = G_model(x, G0, L, TK)
    dG = np.zeros((nc, nz))
    MU = np.zeros((nc, nz))

    for k in range(nc - 1):
        dG[k] = partial_derivative(G_model, x, k, args=(G0, L, TK))

    for k in range(nc):
        MU[k] = G + dG[k] - sum(x_[i]*dG[i] for i in range(nc))

    return MU


# =============================================================================
# Equilibrium vacancy fraction

def LkV_fun(GfV, k, TK):
    """
    Compute thermodynamic interaction parameter between vacancy and atom k.

    Parameters
    ----------
    GfV : dict of lists
        Gibbs free energy of vacancy formation in pure metals.

        ``{k: [enthalpy, entropy] for k in comps}``
    k : str
        Atom.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Interaction parameter.

    """
    GfV_at_T = GfV[k][0] - TK*GfV[k][1]
    return GfV_at_T - GV0


def yVa_model(x, G0, L, TK):
    r"""
    Compute composition-dependent equilibrium vacancy fractions.

    Use a Redlich-Kister expansion with binary interactions of order 0 and 1
    and ternary interactions of order 0:

    .. math::

       &y_{0}^{eq} = \exp \left( - \frac{G_{f,Va}^{\text{alloy}}}{RT} \right)

       &G_{f,Va}^{\text{alloy}} = G_{0}
       + \sum_{i = 1}^{n}{y_i \left(^0L_{0i} - y_i\ ^1L_{0i} \right)}
       - \sum_{i = 1}^{n - 1}{
           \sum_{j = i + 1\ }^{n}{
               y_iy_j \left(^0L_{ij} + 2\ ^1L_{ij}(y_i - y_j) \right)}}

       &\qquad{}+ \sum_{i = 1}^{n - 1}{\sum_{j = i + 1\ }^{n}{y_i y_j L_{0ij}}}
               - 2 \sum_{i = 1}^{n - 2}{
                     \sum_{j = i + 1\ }^{n - 1}{
                       \sum_{k = j + 1\ }^{n}{y_i y_j y_k L_{ijk}}}}

    See Gheno 2022 [#Gheno_2022]_ for details.

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    G0 : list of floats
        Gibbs free energy of the endmembers, len `ncomps`.
    L : 1D array
        Interaction parameters, given as a flat array.

        | Number of terms in a n-component system:
        | binary  : 2*C(n, 2) = n*(n - 1)
        | ternary : C(n, 3) = n*(n - 1)*(n - 2)/6
        |
        | In a 2-comp. system, :math:`L = [L_{12}]`
        | In a 3-comp. system, :math:`L = [L_{12}, L_{13}, L_{23}, L_{123}]`
        | where binary terms are couples :math:`L_{ij} = (^0L_{ij}, ^1L_{ij})`

    T : float or int
        Temperature in Kelvin.

    Returns
    -------
    1D array
        Equilibrium vacancy site fraction.

    """
    n = len(G0)
    xn = 1 - sum(x)
    y_ = np.vstack((np.zeros(xn.size), x, xn))

    n_bin = ut.count_combinations(n)['binaries']
    L_bin = L[:n_bin]
    L_ter = L[n_bin:]

    binary = 0
    for i in range(1, n):
        binary += y_[i]*(L_bin[i - 1, 0] - y_[i]*L_bin[i - 1, 1])
    for i in range(1, n - 1):
        for j in range(i + 1, n):
            idx = ut.index_binary_combination(i, j, n)
            binary -= y_[i]*y_[j]*(L_bin[idx, 0]
                                   + 2*(y_[i] - y_[j])*L_bin[idx, 1])

    ternary = 0
    for i in range(1, n - 1):
        for j in range(i + 1, n):
            idx = ut.index_ternary_combination(0, i, j, n)
            ternary += y_[i]*y_[j]*L_ter[idx, 0]
    for i in range(1, n - 2):
        for j in range(i + 1, n - 1):
            for k in range(j + 1, n):
                idx = ut.index_ternary_combination(i, j, k, n)
                ternary += -2*y_[i]*y_[j]*y_[k]*L_ter[idx, 0]

    Gf = G0[0] + binary + ternary

    return np.exp(-Gf/(R*TK))


# =============================================================================
# Vacancy tracer diffusion coefficient

def make_DVa_fun(DT_fun):
    r"""
    Make function to compute vacancy tracer diffusion coefficient.

    .. math::

        D_{0}^{*} = \frac{1}{y_{0}}\sum_{k = 1}^{n}{y_{k}D_{k}^{*}}

    Parameters
    ----------
    DT_fun : function
        Compute tracer diffusion coefficient of atom constituents.

    Returns
    -------
    function
        Compute vacancy tracer diffusion coefficient.

    """
    def fun(y, x):
        DT = DT_fun(x)
        return 1/y[0]*sum(y[1:]*DT)

    return fun


# =============================================================================
# Onsager coefficients

def make_Lfun(DT_fun, TK):
    """
    Generate function to evaluate Onsager coefficients.

    See background in Andersson 1992 [#Andersson_1992]_.

    Arguments
    ---------
    DT_fun: function
        Function evaluating tracer diffusion coefficients on composition array.
    TK: float
        Temperature in Kelvin.

    Returns
    -------
    function
        Function that takes composition array (typical shape (`n_comps`, `nz`))
        as argument, and returns array (shape (`n_comps`, `n_comps`, `nz`)) of
        phenomenological coefficients evaluated on this array (ex:
        ``L[i, j, k]``, where component i diffuses in a gradient of j, and k is
        a space index).

    """
    def fun(c, x):

        DT = DT_fun(x)

        nc, nz = c.shape
        res = np.zeros((nc, nc, nz))
        for k in range(nc):
            for j in range(nc):
                res[k, j, :] = ut.kd(k, j)*c[k]*DT[k]/(R*TK)
        return res

    return fun


# =============================================================================
# Mobility

def RK_pol(x, p):
    r"""
    Redlich-Kister polynomial with binary and ternary interactions of order 0.

    .. math::

       MQ = \sum_{i = 1}^{n}{x_{i}L_{i}}
       + \sum_{i = 1}^{n - 1}{
           \sum_{j = i + 1\ }^{n}{x_{i}x_{j}L_{ij}}}
       + \sum_{i = 1}^{n - 2}{
           \sum_{j = i + 1\ }^{n - 1}{
             \sum_{k = j + 1\ }^{n}{x_{i}x_{j}x_{k}L_{ijk}}}}

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    p : 1D array
        | Parameters of the polynomial, given as a flat array.
        | Number of terms in a n-component system:
        | unary   : C(n, 1) = n
        | binary  : C(n, 2) = n*(n - 1)/2
        | ternary : C(n, 3) = n*(n - 1)*(n - 2)/6
        |
        | In a 2-component system, :math:`L = [L_1, L_2, L_{12}]`
        | In a 3-component system, :math:`L = [L_1, L_2, L_3, L_{12}, L_{13},
                                               L_{23}, L_{123}]`

    Returns
    -------
    1D array
        Polynomial evaluated along composition profile, shape (`nz`,).

    """
    n = x.shape[0] + 1
    xn = 1 - sum(x)
    x_ = np.vstack((x, xn))
    n_bin = ut.count_combinations(n)['binaries']
    a = p[:n]  # unary terms
    b = p[n:n + n_bin]  # binary interactions
    c = p[n + n_bin:]  # ternary interactions

    unary = sum(a[i]*x_[i] for i in range(n))

    binary = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            idx = ut.index_binary_combination(i, j, n)
            binary += x_[i]*x_[j]*b[idx]

    ternary = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                idx = ut.index_ternary_combination(i, j, k, n)
                ternary += x_[i]*x_[j]*x_[k]*c[idx]

    return unary + binary + ternary


def lnDT_model(x, p, TK):
    """
    Compute log of tracer diffusion coefficients.

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    p : list of floats
        Mobility parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    2D array
        Log of DT, shape (`ncomps`, `nz`).

    """
    nc = x.shape[0] + 1
    nx = x.shape[1]
    lnDT = np.zeros((nc, nx))
    npa = int(len(p)/nc)
    for i in range(nc):
        p_i = p[i*npa:(i + 1)*npa]
        lnDT[i] = unit_lnDT(x, p_i, TK)
    return lnDT


def unit_lnDT(x, p, TK):
    """
    Compute log of tracer diffusion coefficient of one component.

    Use Redlich-Kister polynomial with binary and ternary interactions of order
    0 (see :func:`RK_pol`).

    Parameters
    ----------
    x : 2D array
        Atom fraction of the independent constituents, shape (`ninds`, `nz`).
    p : list of floats
        Mobility parameters.
    TK : float
        Temperature in Kelvin.

    Returns
    -------
    1D array
        Log of DT, shape (`nz`,).

    """
    MQ = RK_pol(x, p)
    return MQ/(R*TK)
