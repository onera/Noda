# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define initial conditions."""

import numpy as np
import scipy.special as sp

import noda.utils as ut
import noda.composition_variables as cv

class InitialConditions:
    """
    Provide initial composition profiles.

    Make atom fraction, vacancy site fraction and pore volume fraction
    profiles. These are then used to initialize cvar, an instance of
    :class:`composition_variables.CompositionVariables`.

    Attributes
    ----------
    comps : list of str
        System components.
    inds : list of str
        Independent components.
    space : :class:`space.SpaceGrid`
        Space grid.
    work_dir : pathlib.Path
        Work directory.
    min_atom_fraction : float
        Minimum atom fraction accepted.
    logger : :class:`log_utils.CustomLogger`
        Logger.
    x : dict of 1D arrays
        Atom fractions.
    yVa : 1D array
        Vacancy site fractions.
    fm : 1D array
        Metal volume fraction profile.
    cvar : :class:`composition_variables.CompositionVariables`
        Composition variable, stores all composition data.

    """
    def __init__(self, params, V_partial, space,
                 work_dir, min_atom_fraction, thermo, logger):
        """
        Class constructor.

        Parameters
        ----------
        params : dict
            Input parameters related to initial conditions.
        V_partial : dict
            Partial molar volumes. See See :func:`data_io.get_volume_data`.
        space : :class:`space.SpaceGrid`
            Space grid.
        work_dir : pathlib.Path
            Work directory.
        min_atom_fraction : float
            Minimum atom fraction accepted.
        thermo : :class:`thermodynamics.Thermodynamics`
            Thermodynamic properties.
        logger : :class:`log_utils.CustomLogger`
            Logger.

        Raises
        ------
        :class:`utils.UserInputError`
            If custom pore fraction or vacancy fraction profile is provided and
            lattice is ideal (vacancy site fraction maintained at equilibrium).

        """
        self.comps = thermo.comps
        self.inds = self.comps[1:-1]
        self.space = space
        self.work_dir = work_dir
        self.min_atom_fraction = min_atom_fraction
        self.logger = logger

        # Make profiles
        self.x = self.make_x_profile(params["atom_fraction"])
        for key in ["vacancy_fraction", "pore_fraction"]:
            if key in params and thermo.ideal_lattice:
                msg = (f"Custom initial '{key}' profile not compatible with "
                       "ideal lattice.")
                raise ut.UserInputError(msg) from None
        self.yVa = self.make_yVa_profile(params, thermo.yVa_fun)
        self.fm = self.make_fm_profile(params)
        self.cvar = cv.CompositionVariables(self.comps, self.x, self.yVa,
                                            V_partial, self.fm)

    def make_x_profile(self, params):
        """Make initial atom fraction profile."""
        if 'file' in params:
            x = self.make_profile("atom_fraction", params)
            x = self.check_x_profile_from_file(x)
        else:
            params = self.prepare_x_params(params)
            x = self.make_profile("atom_fraction", params)
        x = {k: x[k] for k in self.inds}
        return x

    def check_x_profile_from_file(self, x):
        """
        Check initial atom fraction profile.

        * Make sure the components match the independent components declared
          in the input file.
        * Enforce bounds on initial atom fractions and print warning.

        """
        try:
            assert set(list(x)) == set(self.inds)
        except AssertionError as exc:
            msg = ("Missing or extra element in initial atom fraction profile."
                   f"User input contains data for {list(x)}.\n"
                   f"Declared independent components: {self.inds}.")
            raise ut.UserInputError(msg) from exc

        min_val = self.min_atom_fraction
        max_val = 1 - self.min_atom_fraction
        for k in x:
            prof = x[k]
            if any(prof < min_val) or any(prof > max_val):
                msg = (f"Some initial atom fractions were replaced"
                       f" by minimum allowed {min_val}"
                       f" or maximum allowed {max_val}.")
                self.logger.info(msg)
                prof[:] = np.clip(prof, min_val, max_val)
        return x

    def prepare_x_params(self, params):
        """
        Prepare initial atom fraction profile.

        * Make sure the constituents in the initial atom fraction profile match
          the list of independent constituents declared in the configuration.
        * Enforce bounds on initial atom fractions and print warning.

        """
        min_val = self.min_atom_fraction
        max_val = 1 - self.min_atom_fraction
        sides = [ side for side in ['left', 'right'] if side in params]
        msg = "Reading initial conditions for 'atom_fraction'. "
        for side in sides:
            x_side = params[side]
            for k in x_side:
                if x_side[k] < min_val:
                    msg = (f"{side} {k} = {x_side[k]} replaced "
                           f"by minimum allowed {min_val}.")
                    self.logger.info(msg)
                    x_side[k] = min_val
                if x_side[k] > max_val:
                    msg = (f"{side} {k} = {x_side[k]} replaced "
                           f"by maximum allowed {max_val}.")
                    self.logger.info(msg)
                    x_side[k] = max_val
            try:
                assert set(params[side]) == set(self.inds)
            except AssertionError as exc:
                side_comps = ", ".join(params[side])
                inds = ", ".join(self.inds)
                msg = ("Missing or extra element(s) in initial atom fraction "
                       "profile.\n"
                       f"Components on {side} side: {side_comps}\n"
                       f"Declared independent components: {inds}")
                raise ut.UserInputError(msg) from exc
        return params

    def make_yVa_profile(self, params, yVa_fun):
        """Make initial vacancy site fraction profile."""
        if "vacancy_fraction" in params:
            yVa = self.make_profile("vacancy_fraction",
                                    params["vacancy_fraction"])
            yVa = yVa['key']
        else:
            x_arr = np.array([self.x[k] for k in self.inds])
            yVa = yVa_fun(x_arr)
        return yVa

    def make_fm_profile(self, params):
        """Make initial metal volume fraction profile."""
        if "pore_fraction" in params:
            fp = self.make_profile("pore_fraction", params["pore_fraction"])
            fm = 1 - fp['key']
        else:
            fm = np.ones(self.space.nz_init - 1)
        return fm

    def make_profile(self, var, params):
        """
        Make initial profile of variable var.

        Parameters
        ----------
        var : str
            Variable name (atom_fraction, pore_fraction or vacancy_fraction).
        params : dict
            User input parameters.

        Raises
        ------
        :class:`utils.UserInputError`
            If params contains both 'file' and 'shape' parameters, or neither.

        Returns
        -------
        prof : dict of 1D arrays
            Initial profile of variable of interest.

        """
        msg = f"Reading initial conditions for '{var}'. "
        cond1 = "file" in params
        cond2 = "shape" in params
        if (cond1 and cond2) or (not cond1 and not cond2):
            msg += "Expecting either 'file' or 'shape' parameter."
            raise ut.UserInputError(msg) from None
        if "file" in params:
            prof = self.make_profile_from_file(params, msg)
        else:
            prof = self.make_profile_from_dict(var, params, msg)
        return prof

    def make_profile_from_file(self, params, msg):
        """
        Make initial profile from input file.

        The file is read using np.genfromtxt. It should have the independent
        component profiles arranged as columns, with the components as column
        names. The size of the columns must match that of zm (i.e., `nz - 1`).

        Parameters
        ----------
        params : dict
            User input parameters.
        msg : str
            Message indicating which variable is handled.

        Raises
        ------
        :class:`utils.UserInputError`
            If left or right parameter is specified in input dict.
        :class:`utils.UserInputError`
            If the profile length is not compatible with the space grid.

        Returns
        -------
        prof : dict of 1D arrays
            Initial profile of variable of interest.

        """
        if "left" in params or "right" in params:
            msg += ("'file' parameter provided. "
                    "Cannot specify 'left' or 'right' dictionary.")
            raise ut.UserInputError(msg) from None
        fpath = self.work_dir / params["file"]
        arr = np.genfromtxt(fpath, names=True, delimiter=',')
        zm = self.space.zm_init
        try:
            assert arr.shape[0] == zm.size
        except AssertionError as exc:
            msg = ("Initial atom fraction profile provided in "
                   f"{params['file']} has size {arr.shape[0]}. It is not "
                   f"compatible with initial grid of size {zm.size} "
                   f"(nz = {zm.size + 1} nodes i.e. {zm.size} volumes).")
            raise ut.UserInputError(msg) from exc
        prof = {k: arr[k] for k in arr.dtype.names}
        return prof

    def make_profile_from_dict(self, var, params, msg):
        """Make profile with prescribed shape (step, smooth step, flat)."""
        shape = params['shape']
        msg += f"Profile shape is '{shape}'. "
        if shape in ['step', 'smooth step']:
            prof = self.make_step_profile(var, params, msg, shape)
        elif shape == 'flat':
            prof = self.make_flat_profile(params, msg)
        else:
            msg += ("Invalid profile shape. Expecting 'flat', 'step' or "
                    "'smooth step'.")
            raise ut.UserInputError(msg) from None
        return prof

    def make_step_profile(self, var, params, msg, shape):
        """
        Make step or smooth step profile from dict.

        Check input validity and call :func:`make_step_profile` or
        :func:`make_smooth_step_profile`. The step profile is obtained
        with an heaviside function, the smooth step profile with an error
        function (see doc).
        The step position is specified using 'step_fraction' (fraction of the
        domain size, between 0 and 1) or 'step_position' (absolute position in
        m, must be within domain).
        The left and right end values are read from the 'left' and 'right'
        parameters, which must be dicts with the independent components as
        keys.

        Parameters
        ----------
        var : str
            Variable name.
        params : dict
            User input parameters.
        msg : str
            Message indicating which variable is handled.
        shape : str
            Profile shape, can be 'step' or 'smooth step'.

        Raises
        ------
        :class:`utils.UserInputError`
            If parameter 'left' or 'right' is missing.
        :class:`utils.UserInputError`
            If a value is smaller than 0 or larger than 1.
        :class:`utils.UserInputError`
            If parameters 'step_fraction' and 'step_position' are missing, or
            if both are given ;
        :class:`utils.UserInputError`
            If 'step_fraction' is smaller than 0 or larger than 1.
        :class:`utils.UserInputError`
            If 'step_position' is outside domain.

        Returns
        -------
        prof : dict of 1D arrays
            Initial profile of variable of interest.

        """
        zm = self.space.zm_init
        zmin, zmax = self.space.z_init[[0, -1]]
        try:
            left = params['left']
            right = params['right']
        except KeyError:
            msg += "Parameters 'left' and 'right' are required."
            raise ut.UserInputError(msg) from None
        if var in ["vacancy_fraction", "pore_fraction"]:
            for side in [left, right]:
                side = {'key' : side}
        for side_name, side in zip(['left', 'right'], [left, right]):
            for k, v in side.items():
                if not 0 <= v <= 1:
                    msg += (f"Value of {k} on {side_name} side must be "
                            "between 0 and 1.")
                    raise ut.UserInputError(msg) from None
        cond1 = "step_fraction" in params
        cond2 = "step_position" in params
        if (cond1 and cond2) or (not cond1 and not cond2):
            msg += ("Expecting either 'step_fraction' or "
                    "'step_position'.")
            raise ut.UserInputError(msg) from None
        if "step_fraction" in params:
            if (params["step_fraction"] < 0
                or params["step_fraction"] > 1):
                msg += "'step_fraction' should be between 0 and 1."
                raise ut.UserInputError(msg) from None
            zstep = zmin + params["step_fraction"] * (zmax - zmin)
        else:
            if (params["step_position"] < zmin
                or params["step_position"] > zmax):
                msg += ("'step_position' should be between zmin "
                        "({zmin}) and zmax ({zmax}).")
                raise ut.UserInputError(msg) from None
            zstep = params["step_position"]
        if shape == 'step':
            prof = make_step_profile(zm, zstep, left, right)
        else:
            prof = make_smooth_step_profile(zm, zstep, left, right)
        return prof

    def make_flat_profile(self, params, msg):
        """
        Make flat profile from dict.

        The value is read from the 'left' parameter, which must be a dict with
        the independent components as keys.

        Parameters
        ----------
        params : dict
            User input parameters.
        msg : str
            Message indicating which variable is handled.

        Raises
        ------
        :class:`utils.UserInputError`
            If parameter 'left' is missing.
        :class:`utils.UserInputError`
            If parameter 'right' is provided.

        Returns
        -------
        prof : dict of 1D arrays
            Initial profile of variable of interest.

        """
        zm = self.space.zm_init
        try:
            left = params['left']
        except KeyError:
            msg += "Parameter 'left' is required."
            raise ut.UserInputError(msg) from None
        prof = {k: np.ones(zm.shape)*v for k, v in left.items()}
        if 'right' in params:
            msg += "Cannot specify parameter 'right'."
            raise ut.UserInputError(msg) from None
        return prof


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
