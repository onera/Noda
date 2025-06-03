# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define space grids and composition profiles."""

import numpy as np
import scipy.special as sp

import noda.utils as ut


def geo_grid(zmin, zmax, nz, q):
    """
    Make space grid according to geometric progression.

    Parameters
    ----------
    zmin : float
        Start position.
    zmax : float
        End position.
    nz : int
        Number of elements.
    q : float
        Common ratio.

    Returns
    -------
    1D array
        Geometric grid.

    """
    dz0 = (zmax - zmin)*(1 - q)/(1 - q**(nz - 1))
    idx = np.arange(nz)
    return zmin + dz0*(1 - q**idx)/(1 - q)


def make_step_profile(z, zstep, x_left, x_right):
    """
    Make step profiles.

    Parameters
    ----------
    z: 1D array
        Node positions, shape (`nz`,).
    zstep: float
        Position of the step. If `zstep` falls on the `z` grid, the value at
        `zstep` is the average of the end-values.
    x_left: dict of floats
        Atom fractions on the left-hand side.
    x_right: dict of floats
        Atom fractions on the right-hand side.

    Returns
    -------
    x: dict of 1D arrays
        Atom fraction profiles on `z` grid.
    """
    x = {}
    for i in x_left:
        xl = x_left[i]
        xr = x_right[i]
        # Round to avoid floating point precision error, which can place zstep
        # off the z grid when it should in fact be on a grid point.
        z_rel = np.around(z - zstep, decimals=16)
        x[i] = xl + (xr - xl)*np.heaviside(z_rel, 0.5)
    return x


def make_smooth_step_profile(z, zstep, x_left, x_right):
    """
    Make smooth step profiles.

    Parameters
    ----------
    z: 1D array
        Node positions, shape (`nz`,).
    zstep: float
        Position of the step. If `zstep` falls on the `z` grid, the value at
        `zstep` is the average of the end-values.
    x_left: dict of floats
        Atom fractions on the left-hand side.
    x_right: dict of floats
        Atom fractions on the right-hand side.

    Returns
    -------
    x: dict of 1D arrays
        Atom fraction profiles on `z` grid.
    """
    zspan = z[-1] - z[1]
    x = {}
    for i in x_left:
        xl = x_left[i]
        xr = x_right[i]
        # Round to avoid floating point precision error, which can place zstep
        # off the z grid when it should in fact be on a grid point.
        z_rel = np.around(z - zstep, decimals=16)
        x[i] = (xl + xr)/2 + (xr - xl)/2 * sp.erf(z_rel/zspan*z.size)
    return x


def make_x_full(x_sides, nx=20):
    """
    Map the composition space.

    Only supports binary and ternary systems.

    Parameters
    ----------
    x_sides : dict of dicts
        End values of the atom fractions of each constituent.
    nx : int, optional
        Number of steps in each dimension of the composition grid. The default
        is 20.

    Raises
    ------
    Exception
        If the system contains more than two independent constituents.

    Returns
    -------
    x_full : 2D array
        Atom fraction grid, shape (`n_inds`, `nx**n_inds`).

    """
    # TIP: this function is no longer used. Max DT is calculated from initial x
    # array (see simu.add_time). It would be useful to extend make_x_full to
    # systems of any size to better scan possible compositions (?).

    x_range = {}
    for k in x_sides['left']:
        x_range[k] = np.linspace(x_sides['left'][k], x_sides['right'][k],
                                 num=nx)

    n_inds = len(x_sides['left'])

    if n_inds == 1:
        x_full = np.array(list(x_range.values()))

    elif n_inds == 2:

        x_full = np.zeros((2, nx**2))

        A = list(x_range)[0]
        B = list(x_range)[1]

        x_full[0] = np.hstack([np.ones(nx)*x_range[A][i] for i in range(nx)])
        x_full[1] = np.hstack([x_range[B] for i in range(nx)])

    else:
        msg = f'"make_x_full" not implemented for n_inds = {n_inds}'
        raise ut.UserInputError(msg)

    return x_full


def make_xtab(x):
    """
    Make 2D array from dict of 1D arrays.

    Parameters
    ----------
    x : dict of 1D arrays
        `n_inds` composition arrays, shape (`nx_i`,).

    Returns
    -------
    points : 2D array
        Composition array, shape (`product(nx_i)`, `n_inds`).

    """
    grids = np.meshgrid(*list(x.values()), indexing='ij')
    points = np.array(grids).T.reshape(-1, len(x))
    return points
