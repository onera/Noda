# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Miscellaneous functions."""

import re
import itertools as it
from math import comb

import numpy as np


def kd(a, b):
    """
    Kronecker delta.

    Returns 1 if a = b, 0 otherwise.

    """
    if a == b:
        res = 1
    else:
        res = 0
    return res


def format_element_symbol(s):
    """Format string to match valid element symbol."""
    if len(s) == 1:
        return s.upper()
    if len(s) == 2:
        return s[0].upper() + s[1].lower()
    msg = f"Invalid element symbol '{s}' in input file."
    raise UserInputError(msg) from None


def get_or_raise(dictionary, key):
    """Get value or raise exception if key not present in dictionary."""
    try:
        val = dictionary[key]
    except KeyError:
        msg = f"Entry '{key}' not found in 'user_data.toml' file."
        raise UserInputError(msg) from None
    return val

# =============================================================================
# Combinations

def count_combinations(n):
    """
    Count p-size combinations in n-size set, with p up to 3.

    Parameters
    ----------
    n : int
        Size of set.

    Returns
    -------
    res : dict of ints
        ``{k: val for k in keys}`` with the following keys:

        * 'unaries' : n.
        * 'binaries' : C(n, 2).
        * 'ternaries' : C(n, 3).
        * 'mix' : binaries + ternaries.
        * 'all' : unaries + binaries + ternaries.
    """
    res = {}
    res['unaries'] = n
    res['binaries'] = comb(n, 2)
    res['ternaries'] = comb(n, 3)
    res['mix'] = res['binaries'] + res['ternaries']
    res['all'] = res['unaries'] + res['mix']
    return res


def make_combinations(constituents):
    """
    Make combinations of constituents up to size 3 (unary, binary, ternary).

    The combinations are stored as concatenated strings, for example "AuAg".

    Parameters
    ----------
    constituents : list of str
        Constituents.

    Returns
    -------
    res : dict of lists
        ``{k: list_of_strings for k in keys}`` with the following keys:

        * 'unaries' : Constituents.
        * 'binaries' : All binary subsystems.
        * 'ternaries' : All ternary subsystems.
        * 'mix' : All binary and ternary subsystems.
        * 'all' : All subsystems.

    """
    res = {}
    res['unaries'] = constituents
    res['binaries'] = [''.join(c) for c in it.combinations(constituents, 2)]
    res['ternaries'] = [''.join(c) for c in it.combinations(constituents, 3)]
    res['mix'] = res['binaries'] + res['ternaries']
    res['all'] = res['unaries'] + res['binaries'] + res['ternaries']
    return res


def index_binary_combination(i, j, n):
    """
    Compute index of binary combination in n-size system.

    Combination (i, j) with i > j.

    Example with n = 3:

    = = =====  
    i j index  
    = = =====  
    0 1 0      
    0 2 1      
    1 2 2      
    = = =====  

    """
    i, j = sorted([i, j])
    res = comb(n, 2) - 1
    res -= comb(n - i - 1, 2) + (n - j - 1)
    return res


def index_ternary_combination(i, j, k, n):
    """
    Compute index of ternary combination in n-size system.

    Combination (i, j, k) with i > j > k.

    Example with n = 4:

    = = = =====
    i j k index
    = = = =====
    0 1 2 0
    0 1 3 1
    0 2 3 2
    1 2 3 3
    = = = =====

    """
    i, j, k = sorted([i, j, k])
    res = comb(n, 3) - 1
    res -= comb(n - i - 1, 3) + comb(n - j - 1, 2) + (n - k - 1)
    return res


def make_permutations_samesize(solvent):
    """
    Make all permutations of constituents in solvent, preserving size.

    Parameters
    ----------
    solvent : str
        Constituents concatenated to string. A constituent is defined by one
        uppercase letter followed by any number of lowercase letters.

    Returns
    -------
    res : list of str
        All permutations, with same syntax as input.

    """
    elements = re.findall('[A-Z][a-z]*', solvent)
    res = [''.join(tup) for tup in it.permutations(elements, len(elements))]
    return res


# =============================================================================
# Integration, derivation

def integrate(f, u, F0, geometry):
    """
    Compute radial vector field from its divergence.

    The field is evaluated along one space coordinate, in either of 3 geometric
    configurations (planar, cylindrical, spherical).

    The `f` values are evaluated at midpoints of the `u` array.
    Ensures minimum error when composed with :func:`div`.

    Parameters
    ----------
    f : 1D array
        Divergence of the radial field, shape (u.size - 1,).
    u : 1D array
        Coordinates where the field is evaluated, shape (u.size,).
    F0 : float
        Integration constant.
    geometry : int
        Indicates the geometry (planar, cylindrical, spherical).

    Returns
    -------
    1D array
        Radial field, shape (u.size,).

    """
    du = np.diff(u)
    A = np.zeros(f.size + 1)
    res = None
    A[0] = F0
    if geometry == 'planar':
        for i in range(f.size):
            A[i + 1] = A[i] + f[i]*du[i]
        res = A
    elif geometry == 'cylindrical':
        um = (u[1:] + u[:-1])/2
        for i in range(f.size):
            A[i + 1] = A[i] + um[i]*f[i]*du[i]
        res = np.divide(A, u, out=A, where=(u!=0))
    elif geometry == 'spherical':
        um = (u[1:] + u[:-1])/2
        for i in range(f.size):
            A[i + 1] = A[i] + um[i]**2*f[i]*du[i]
        res = np.divide(A, u**2, out=A, where=(u!=0))
    return res


def div(F, u, geometry):
    """
    Compute the divergence of a radial vector field.

    The field is evaluated along one space coordinate, in either of 3 geometric
    configurations (planar, cylindrical, spherical).

    Ensures minimum error when composed with :func:`integrate`.

    Parameters
    ----------
    F : 1D array
        Radial field, shape (u.size,).
    u : 1D array
        Coordinates where the field is evaluated, shape (u.size,).
    geometry : int
        Indicates the geometry (planar, cylindrical, spherical).

    Returns
    -------
    1D array
        Divergence of the field, shape (u.size - 1,).

    """
    res = None
    if geometry == 'planar':
        res = np.diff(F)/np.diff(u)
    elif geometry == 'cylindrical':
        um = (u[1:] + u[:-1])/2
        res = 1/um*np.diff(u*F)/np.diff(u)
    elif geometry == 'spherical':
        um = (u[1:] + u[:-1])/2
        res = 1/um**2*np.diff(u**2*F)/np.diff(u)
    return res


# =============================================================================
# Custom errors
# =============================================================================


class AtomFractionError(Exception):
    """Exception raised when an atom fraction is invalid."""

    def __init__(self, n):
        msg = (f"Invalid atom fraction at step {n}. Try decreasing the time "
               "step\n(use dt_multiplier < 1 or nt_multiplier > 1).")
        super().__init__(msg)


class UserInputError(Exception):
    """Exception raised when a user input is invalid."""

    def __init__(self, msg):
        super().__init__(msg)


class ResultsError(Exception):
    """Exception raised when a user request to access results is invalid."""

    def __init__(self, msg):
        super().__init__(msg)
