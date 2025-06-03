# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gather and organize the variables that describe a system composition."""

import numpy as np
import scipy.interpolate as spi

def vec_to_nod(x, z):
    """
    Evaluate `x` (defined on midpoints) on nodes by linear interpolation.

    Parameters
    ----------
    x : 1D array
        Variable, shape (`nz` - 1).
    z : 1D array
        Node positions, shape (`nz`,).

    Returns
    -------
    arr : 1D array
        Variable evaluated on nodes, shape (`nind`, `nz`).

    """
    zm = (z[:-1] + z[1:])/2
    f = spi.interp1d(zm, x, kind='linear', fill_value='extrapolate')
    return f(z)


def arr_to_nod(x, z, interpkind='cubic'):
    """
    Evaluate `x` (defined on midpoints) on nodes by linear interpolation.

    Parameters
    ----------
    x : 2D array
        Variable, shape (`nind`, `nz` - 1).
    z : 1D array
        Node positions, shape (`nz`,).

    Returns
    -------
    arr : 2D array
        Variable evaluated on nodes, shape (`nind`, `nz`).

    """
    nind = x.shape[0]
    arr = np.zeros((nind, z.size))
    zm = (z[:-1] + z[1:])/2
    for i, row in enumerate(x):
        f = spi.interp1d(zm, row, kind='linear', fill_value='extrapolate')
        arr[i] = f(z)
    return arr


class MultiVariable:
    """
    Vector variable evaluated on midpoints `zm`.

    One row per independent constituent.
    """

    def __init__(self, v_init):
        """
        Class constructor.

        Parameters
        ----------
        v_init : dict
            Initial value, dict of 1D arrays of shape (`nz` - 1,), with
            independent constituents as keys.

        Attributes
        ----------
        mid : 2D array
            Variable evaluated on midpoints, shape (`ninds`, `nz` - 1)

        """
        self.mid = np.array(list(v_init.values()))

    def nod(self, z):
        """Evaluate on node, see :func:`arr_to_nod`."""
        return arr_to_nod(self.mid, z)


class UniVariable:
    """Scalar variable evaluated on midpoints `zm`."""

    def __init__(self, v_init):
        """
        Class constructor.

        Parameters
        ----------
        v_init : dict
            Initial value, 1D array of shape (`nz` - 1,).

        Attributes
        ----------
        mid : 1D array
            Variable evaluated on midpoints, shape (`nz` - 1,)

        """
        self.mid = v_init

    def nod(self, z):
        """
        Evaluate on nodes (`z`) by interpolation.

        Parameters
        ----------
        z : 1D array
            Node positions, shape (`nz`,).

        Returns
        -------
        1D array
            Variable evaluated on nodes, shape (`nz`,).

        """
        zm = (z[:-1] + z[1:])/2
        f = spi.interp1d(zm, self.mid, fill_value='extrapolate')
        return f(z)


class CompositionVariables:
    """
    Contain and organize composition variables.

    This describes the composition of a system with one metal phase and one
    pore phase.

    The system composition is initialized with metal atom fractions. Once all
    variables are built, the reference variable is the system concentration
    (in mol/m3). Only the system concentration will be updated in the
    calculations. The other variables will derive from the system
    concentrations --- the other variables are 'read-only'.

    Attributes
    ----------
    x : 2D array
        Metal atom fractions.
    y : 2D array
        Metal site fractions.
    V : 1D array
        Average system molar volume.
    Vm : 1D array
        Average metal molar volume.
    fm : 1D array
        Metal volume fraction.
    fp : 1D array
        Pore volume fraction.
    c : 2D array
        System concentrations.

    """

    # pylint: disable=too-many-instance-attributes

    def __init__(self, comps, x_init, yVa_init, V_partial, fm):
        """
        Class constructor.

        | Build metal atom fractions.
        | Build metal site fractions.
        | Build average molar volume of metal (`Vm`), then of system (`V`).
        | Build system concentration.

        Parameters
        ----------
        comps : list of str
            System constituents, ordered: ['Va'] + inds + [dep].
        x_init : dict
            Initial atom fractions, dict of 1D arrays of shape (`nz` - 1,).
        yVa_init : 1D array
            Initial vacancy site fraction, shape (`nz` - 1,).
        V_partial : dict
            Partial molar volumes.
        fm : 1D array
            Initial metal volume fraction, shape (`nz` - 1,).

        """
        inds = comps[1:-1]
        dep = comps[-1]
        self.V_partial = V_partial

        x_dep = 1 - sum(x_init.values())
        x_full = {k: x_init[k] for k in inds}
        x_full[dep] = x_dep

        y_full = {k: x_full[k]*(1 - yVa_init) for k in comps[1:]}
        y_full['Va'] = yVa_init

        if V_partial['Va'] == 'local':
            Vm = sum(x_full[k]*V_partial[k] for k in comps[1:])
            Vk_arr = np.hstack((np.nan,
                                np.array([V_partial[k] for k in comps[1:]])
                                ))
        else:
            Vm = sum(y_full[k]*V_partial[k] for k in comps)
            Vk_arr = np.array([V_partial[k] for k in comps])

        self.Vk_arr = Vk_arr[np.newaxis].T
        V = Vm/fm
        c_full = {k: y_full[k]/V for k in comps}

        self.comps = comps
        self._x = MultiVariable(x_full)
        self._y = MultiVariable(y_full)
        self.c = MultiVariable(c_full)
        self._Vm = UniVariable(Vm)
        self._fm = UniVariable(fm)
        self._fp = UniVariable(1 - fm)
        self._V = UniVariable(V)

    @property
    def x(self):
        """
        Atom fractions evaluated on midpoints.

        Returns
        -------
        2D array
            Shape (`ninds` + 1, `nz` - 1).
        """
        self._x.mid = self.y.mid[1:]/(1 - self.y.mid[0])
        return self._x

    @property
    def y(self):
        """
        Site fractions of all metal constituents evaluated on midpoints.

        Returns
        -------
        2D array
            Shape (`ninds` + 2, `nz` - 1).

            | Order of rows:
            |   0 : vacancies
            |   1.. : independent atom constituents
            |   n : dependent constituent
        """
        self._y.mid = self.c.mid/sum(self.c.mid)
        return self._y

    @property
    def V(self):
        """Average molar volume of the system, 1D array of shape (`nz` - 1)."""
        self._V.mid = 1/sum(self.c.mid)
        return self._V

    @property
    def Vm(self):
        """
        Average molar volume of the metal, 1D array of shape (`nz` - 1).

        Choice of two definitions depending on whether the vacancy molar volume
        is independent or that of the local environment.
        """
        if self.V_partial['Va'] == 'local':
            self._Vm.mid = sum(self.x.mid*self.Vk_arr[1:])
        else:
            self._Vm.mid = sum(self.y.mid*self.Vk_arr)
        return self._Vm

    @property
    def fm(self):
        """Metal volume fraction, 1D array of shape (`nz` - 1)."""
        self._fm.mid = self.Vm.mid/self.V.mid
        return self._fm

    @property
    def fp(self):
        """Pore volume fraction, 1D array of shape (`nz` - 1)."""
        self._fp.mid = 1 - self.fm.mid
        return self._fp
