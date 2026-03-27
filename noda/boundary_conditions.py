# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Handle boundary conditions."""

import numpy as np

import noda.utils as ut
import noda.composition_variables as cv


class BoundaryConditions:
    """
    Boundary conditions on one side of the domain (left or right).

    Can be either:

    * 'Dirichlet': prescribed composition, if 'atom_fractions' is provided in
      the input parameters.
    * 'Neumann': prescribed flux, if 'flux' is provided.

    Defaults to 0-flux if none is specified in the input parameters.
    
    The class is used through cvar_fun and J_fun, which are both function of
    time built from the input parameters, which may be a float or a str with
    t as the time variable (ex: (3*t + 2)**(1/2)).

    Attributes
    ----------
    thermo : :class:`thermodynamics.Thermodynamics`
        Thermodynamic properties handler.
    comps : list of str
        System components.
    inds : list of str
        Independent components.
    V_partial : dict
        Partial molar volumes. See :func:`data_io.get_volume_data`.
    min_atom_fraction : float
        Minimum atom fraction accepted.
    logger : :class:`log_utils.CustomLogger`
        Logger.
    side : str
        Side where boundary condition applies (left or right).
    type : str
        Type of boundary condition (Neumann or Dirichlet).
    
    Methods
    -------
    cvar_fun(t) :
        Function of time, returns a
        :class:`composition_variables.CompositionVariables` instance.
    J_fun(t)
        Function of time, returns a flux array.

    """
    def __init__(self, params, thermo, V_partial, min_atom_fraction, logger,
                 side):
        """Class constructor."""
        self.thermo = thermo
        self.comps = thermo.comps
        self.inds = self.comps[1:-1]
        self.V_partial = V_partial
        self.min_atom_fraction = min_atom_fraction
        self.logger = logger
        self.side = side
        self.type = self.make_BC_type(params)
        xparams = self.make_BC_dict(params, "atom_fraction")
        if self.type == "Dirichlet":
            self.cvar_fun = self.make_cvar_fun(xparams)
        else:
            self.cvar_fun = None
        info_list = []
        Jparams = self.make_BC_dict(params, "flux")
        for k in self.comps[1:]:
            if k not in Jparams:
                Jparams[k] = "0"
                info_list.append(k)
        if self.type == "Neumann":
            self.J_fun = self.make_J_fun(Jparams)
        else:
            self.J_fun = None
        if len(info_list) > 0:
            logger.info("Auto boundary conditions:")
            for k in info_list:
                text = f"* {side:5} BC for {k} set to 0 flux"
                logger.info(text)

    def make_BC_type(self, params):
        """Guess BC type from input dict and make sure input is consistent."""
        if len(params) > 1:
            msg = (f"{self.side} boundary condition: cannot specify more than "
                   "one variable.")
            raise ut.UserInputError(msg)

        if 'atom_fraction' in params:
            res = 'Dirichlet'
        elif ('flux' in params or params == {}):
            res = 'Neumann'
        else:
            var = list(params.keys())[0]
            msg = f"{self.side} boundary condition: invalid variable {var}."
            raise ut.UserInputError(msg)
        return res

    def make_BC_dict(self, params, var):
        """Format input dictionary."""
        dct = params.get(var, {})
        return {ut.format_element_symbol(k) : str(v) for k, v in dct.items()}

    def make_J_fun(self, Jparams):
        """
        Make function that returns a boundary flux array.

        The inner function returns a 1D array of fluxes, which include
        all System components.

        """
        self.check_BC_components(Jparams.keys(), self.comps[1:], 'flux')
        fun_dct = {k: lambda t, s=Jparams[k]: eval(s) for k in self.comps[1:]}
        def fun(t):
            J_dict = {k: fun_dct[k](t) for k in self.comps[1:]}
            J_arr = np.array(list(J_dict.values()))
            return J_arr
        return fun

    def make_cvar_fun(self, xparams):
        """
        Make function that returns a boundary composition variable.

        The inner function returns a
        :class:`composition_variables.CompositionVariables` instance. Two
        assumptions are made:

        * vacancies are at equilibrium,
        * the pore fraction is 0.

        """
        self.check_BC_components(xparams.keys(), self.inds,
                                      'atom_fraction')
        xparams = self.clip_xBC_values(xparams)
        fun_dct = {k: lambda t, s=xparams[k]: eval(s) for k in self.inds}

        def fun(t):
            x_dict = {k: fun_dct[k](t) for k in self.inds}
            x_arr = np.array([x_dict[k] for k in self.inds])
            yVa = self.thermo.yVa_fun(x_arr)
            cvar = cv.CompositionVariables(self.comps, x_dict, yVa,
                                           self.V_partial, fm=1)
            return cvar
        return fun

    def check_BC_components(self, found, expected, varname):
        """Make sure all expected components are present in the BC."""
        try:
            assert set(found) == set(expected)
        except AssertionError as exc:
            msg = (f"Missing or extra elements in '{varname}' boundary "
                   "condition.\n"
                   f"expected {expected}\n"
                   f"found {found}")
            raise ut.UserInputError(msg) from exc

    def clip_xBC_values(self, xparams):
        """Make sure atom fractions are within accepted bounds."""
        min_val = self.min_atom_fraction
        max_val = 1 - self.min_atom_fraction
        for k in xparams:
            v = float(xparams[k])
            if v < min_val:
                msg = (f"{self.side} boundary condition : input atom fraction "
                       f"{k} = {v} replaced by minimum allowed {min_val}.")
                self.logger.info(msg)
                xparams[k] = str(min_val)
            if v > max_val:
                msg = (f"{self.side} boundary condition : input atom fraction "
                       f"{k} = {v} replaced by maximum allowed {max_val}.")
                self.logger.info(msg)
                xparams[k] = str(max_val)
        return xparams
