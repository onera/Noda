# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Constants."""

import numpy as np

# =============================================================================
# Universal constants

R = 8.314
"""float: Ideal gas constant in J K-1 mol-1."""

EV = 1.602e-19
"""float: Electronvolt in J/eV."""

NA = 6.022e23    # mol-1
"""float: Avogadro number in mol-1."""

KB = 8.617e-5  # eV/K
"""float: Boltzmann constant in eV/K."""

# =============================================================================
# Thermodynamics

GV0_eV = np.log(2) - 1/2
GV0 = GV0_eV*EV*NA
"""float: Energy of vacancy end-member in J/mol."""

# =============================================================================
# Default parameters

factory_default_parameters = {
   'min_atom_fraction': 1e-9,
   'min_number_time_steps': 20,
   'num_out': 2,
   'grid_type': 'linear',
   'nz': 60,
   'q': 1.02,
   'Fourier_number': 0.25,
   'geometry': 'planar',
   'stencil': 'A',
   'partial_molar_volume': 1e-5,
   'vacancy_formation_energy': [2e5, 50],
    }
"""
dict: Factory default parameters.

'min_atom_fraction' : str
    Minimum atom fraction accepted in initial profiles and boundary conditions.
'min_number_time_steps' : int
    Minimum number of time steps.
'num_out' : int
    Number of saved time steps (see :ref:`time`).
'grid': str
    Type of space grid (see :ref:`space`).
'nz' : int
    Number of space steps (see :ref:`space`).
'q' : float
    Common ratio for geometric grid (see :ref:`space`).
'Fourier_number' : float
    Fourier number (see :ref:`time_step`).
'geometry' : str
    Domain geometry (see :ref:`space`).
'stencil' : str
    Space discretization stencil (see :ref:`stencil`).
'partial_molar_volume' : float
    Partial molar volume in m3/mol.
'vacancy_formation_energy' : list of floats
    Vacancy formation enthalpy and entropy, [HfV, SfV] in [J/mol, J/mol/K]
    GfV = HfV - T*SfV
"""
