# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define temperature program."""

import noda.utils as ut

class Temperature:
    """
    Store temperature parameters.
    
    Attributes
    ----------
    TK : float
        Temperature in Kelvin
    TC : float
        Temperature in Celsius
    
    """
    def __init__(self, params):
        """
        Class constructor.

        Parameters
        ----------
        params : dict
            Temperature-related input parameters.

        Raises
        ------
        :class:`utils.UserInputError`
            If both TC and TK are specified in input.

        """
        if 'TC' in params and 'TK' in params:
            msg = ("Cannot specify both TC and TK (temperature in Celsius and "
                   "Kelvin).")
            raise ut.UserInputError(msg) from None
        if 'TC' in params:
            self.TC = params['TC']
            self.TK = self.TC + 273.15
        elif 'TK' in params:
            self.TK = params['TK']
            self.TC = self.TK - 273.15
