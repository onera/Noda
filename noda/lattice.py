# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compute variables describing the evolution of the crystal lattice."""

import numpy as np

from noda.utils import div, integrate


def compute_alpha_nonideal(dt, yVa, yVa_eq, V, Vp, fm,
                           k_dislo, k_pores):
    """
    Compute lattice sink rate in non-ideal case.

    Parameters
    ----------
    dt : float
        Time step.
    yVa : 1D array
        Vacancy site fraction, shape (`nz` - 1,).
    yVa_eq : 1D array
        Equilibrium vacancy site fraction, shape (`nz` - 1,).
    V : 1D array
        System average molar volume, shape (`nz` - 1,).
    Vp : dict of floats
        Partial molar volumes.
    fm : 1D array
        Metal volume fraction, shape (`nz` - 1,).
    k_dislo : 1D array
        Sink strength associated with dislocation climb, shape (`nz` - 1,).
    k_pores : 1D array
        Sink strength associated with pore growth, shape (`nz` - 1,).

    Returns
    -------
    alpha_d : 1D array
        Sink term associated with dislocation climb, shape (`nz` - 1,).
    alpha_p : 1D array
        Sink term associated with pore growth, shape (`nz` - 1,).

    """
    y_xs = yVa - yVa_eq

    # Creation/annihilation via bulk lattice activity
    alpha_d = -k_dislo*y_xs

    # Creation/annihilation via exchange with pores
    sursat = y_xs > 0
    fp = 1 - fm
    inpore = fp > 0
    contribute_to_fp = sursat | inpore
    alpha_p = np.where(contribute_to_fp, -k_pores*y_xs, 0)

    # This is used to prevent fp < 0
    # Negative fp can occur because of small oscillations in fp: in regions
    # with y_xs < 0, if fp is slightly positive instead of 0, contribute_to_fp
    # will be True and the positive alpha_p will lead to a decrease of fp,
    # possibly below 0.
    # The expression for alpha_p_max is an approximation: it neglects
    # V/Vp*div(fp*v), which cannot be evaluated because v is not known (but
    # this term is small).
    alpha_p_max = V/Vp*fp/dt
    alpha_p = np.clip(alpha_p, -np.inf, alpha_p_max)

    return alpha_d, alpha_p


def compute_gamma(Jlat, z, Vk, V0, Vp, V, alpha_d, alpha_p, y_Va, geometry):
    """
    Compute relative volume variation rate.

    Parameters
    ----------
    Jlat : 1D array
        Fluxes in the lattice frame, shape (`ninds` + 1, `nz`).
    z : 1D array
        Node positions, shape (`nz`,).
    Vk : 1D array
        Partial molar volumes.
    V0 : float or 1D array
        Partial molar volume of the vacancy.
    Vp : float or 1D array
        Partial molar volume of the pore.
    V : 1D array
        System average molar volume, shape (`nz` - 1,).
    alpha_d : 1D array
        Sink term associated with dislocation climb, shape (`nz` - 1,).
    alpha_p : 1D array
        Sink term associated with pore growth, shape (`nz` - 1,).
    y_Va : 1D array
        Vacancy site fraction, shape (`nz` - 1,).
    geometry : str
        Domain geometry (planar, cylindrical or spherical).

    Returns
    -------
    gamma : 1D array
        Relative volume variation rate, shape (`nz` - 1,).

    """
    # TIP: this function is not used anywhere
    alpha = alpha_d + alpha_p
    sum_V_divJ = sum((V0 - Vk[1:])*div(Jlat, z, geometry))
    if not isinstance(V0, float):
        gamma = sum_V_divJ/(1 - y_Va) + (alpha*V0 - alpha_p*Vp)/V
    else:
        gamma = sum_V_divJ + (alpha*V0 - alpha_p*Vp)/V
    return gamma


def compute_gamma_V(Jlat, z, Vk, V0, Vm, V, y_Va, alpha, geometry):
    """
    Compute relative volume variation rate due to molar volume variation.

    Parameters
    ----------
    Jlat : 1D array
        Fluxes in the lattice frame, shape (`ninds` + 1, `nz`).
    z : 1D array
        Node positions, shape (`nz`,).
    Vk : 1D array
        Partial molar volumes.
    V0 : float or 1D array
        Partial molar volume of the vacancy.
    Vm : 1D array
        Average molar volume of the metal phase, shape (`nz` - 1,).
    V : 1D array
        System average molar volume, shape (`nz` - 1,).
    y_Va : 1D array
        Vacancy site fraction, shape (`nz` - 1,).
    alpha : 1D array
        Sink term, shape (`nz` - 1,).
    geometry : str
        Domain geometry (planar, cylindrical or spherical).

    Returns
    -------
    gamma_V : 1D array
        Relative volume variation rate due to molar volume variation, shape
        (`nz` - 1,).

    """
    sum_V_divJ = sum((V0 - Vk[1:])*div(Jlat, z, geometry))
    if not isinstance(V0, float):
        gamma_V = sum_V_divJ/(1 - y_Va)
    else:
        gamma_V = sum_V_divJ + alpha*(V0 - Vm)/V
    return gamma_V


def compute_velocity(gamma, z, Vk, V0, Jlat, yVa, geometry):
    """
    Compute velocity field of the lattice in the laboratory frame.

    `v` is calculated via an integral from 0 to `z` -> need to add value at 0.
    `v_left` is obtained assuming ideal lattice activity on left boundary
    (and neglecting the gradient of `Vm/(1 - yVa)` in getting the primitive).

    Parameters
    ----------
    gamma : 1D array
        Relative volume variation rate, shape (`nz` - 1,).
    z : 1D array
        Node positions, shape (`nz`,).
    Vk : 1D array
        Partial molar volumes.
    V0 : float or 1D array
        Partial molar volume of the vacancy.
    Jlat : 1D array
        Fluxes in the lattice frame, shape (`ninds` + 1, `nz`).
    yVa : 1D array
        Vacancy site fraction, shape (`nz` - 1,).
    geometry : str
        Domain geometry (planar, cylindrical or spherical).

    Returns
    -------
    v : 1D array
        Lattice velocity, shape (`nz`,)

    """
    V0_left = V0 if isinstance(V0, float) else V0[0]
    v_left = -sum((Vk[1:] + V0_left*yVa[0]/(1 - yVa[0]))*Jlat[:, [0]])
    v_left = v_left.item()
    v = integrate(gamma, z, v_left, geometry)
    return v
