# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define initial space grid."""

import numpy as np

import noda.utils as ut
from noda.log_utils import get_and_log


class SpaceGrid:
    """
    Store initial domain dimensions.

    Attributes
    ----------
    default_params : dict
        Parameters used when not specified in user input.
    geometry : str
        Domain geometry ('planar', 'cylindrical' or 'spherical').
    work_dir : pathlib.Path
        Work directory.
    logger : :class:`log_utils.CustomLogger`
        Logger.
    zmin : float
        Initial position of left-hand domain boundary (m).
    zmax : float
        Initial position of right-hand domain boundary (m).
    nz_init : int
        Number of steps in initial space grid.
    z_init : 1D array
        Initial node positions. Ranges from 0 to zmax, size nz. This is where
        fluxes are evaluated.
    dz_init : 1D array
        Initial space step, shape (nz - 1,).
    zm_init : 1D array
        Initial midpoint positions, shape (nz - 1,). This is where
        composition variables are evaluated.

    """
    def __init__(self, params, default_params, work_dir, logger):
        """Class constructor."""
        self.default_params = default_params
        self.work_dir = work_dir
        self.logger = logger
        self.geometry = get_and_log(params, 'geometry',
                                    default_params['geometry'],
                                    logger)
        if self.geometry not in ('planar', 'cylindrical', 'spherical'):
            msg = f"Invalid geometry parameter {self.geometry}."
            raise ut.UserInputError(msg) from None
        self.z_init = self.make_grid(params)
        self.zmin = self.z_init[0]
        self.zmax = self.z_init[-1]
        self.nz_init = self.z_init.size
        self.dz_init = np.diff(self.z_init)
        self.zm_init = (self.z_init[:-1] + self.z_init[1:])/2

    def make_grid(self, params):
        """
        Make initial space grid from input file or parameters.

        If 'file' is given, the grid is read from work_dir/filename using
        np.genfromtxt, and zmin, zmax and nz are inferred from the grid.
        If not, see :meth:`make_grid_from_dict`.

        Raises
        ------
        :class:`utils.UserInputError`
            If params contains both 'file' and space parameters.

        """
        if 'file' in params:
            for x in ['nz', 'zmin', 'zmax', 'grid_type']:
                if x in params:
                    msg = ("Found 'file' in space parameters. "
                           f"Cannot specify {x}.")
                    raise ut.UserInputError(msg) from None
            fpath = self.work_dir / params['file']
            z_init = np.genfromtxt(fpath)
        else:
            z_init = self.make_grid_from_dict(params)
        return z_init

    def make_grid_from_dict(self, params):
        """
        Make grid from parameters.

        The type of grid is specified via the optional parameter 'grid_type',
        which can be 'linear' or 'geometric' and defaults to 'linear' :

        * 'linear': linear grid from zmin to zmax with size nz.
        * 'geometric': geometric grid from zmin to zmax with size nz and common
          ratio q.

        In both cases :

        * 'zmax' must be included in input
        * 'zmin' is optional and defaults to 0
        * 'nz' is optional and defaults to value in default_params (can be set
          in 'user_data.toml').
        * 'q' is optional and defaults to value in default_params.

        Raises
        ------
        :class:`utils.UserInputError`
            If 'zmax' is not found in input dict.
        :class:`utils.UserInputError`
            If 'zmin' or 'zmax' are not compatible with the domain geometry.

        """
        grid_type = get_and_log(params, 'grid_type', 'linear', self.logger)
        zmin = params.get('zmin', 0)
        try:
            zmax = params['zmax']
        except AttributeError as exc:
            msg = "Space parameters : zmax is required."
            raise ut.UserInputError(msg) from exc
        nz = get_and_log(params, 'nz',
                         self.default_params['number_space_steps'],
                         self.logger)
        if grid_type == 'linear':
            z_init = np.linspace(zmin, zmax, num=nz)
        elif grid_type == 'geometric':
            q = get_and_log(params, 'q',
                            self.default_params['common_ratio'],
                            self.logger)
            z_init = geo_grid(zmin, zmax, nz, q)
        else:
            msg = f"Invalid grid_type '{grid_type}'."
            raise ut.UserInputError(msg) from None
        if self.geometry in ['cylindrical', 'spherical']:
            if zmin < 0:
                msg = (f"Found strictly negative zmin in {self.geometry} "
                       f"geometry (zmin = '{zmin}'). "
                       "zmin should be positive.")
                raise ut.UserInputError(msg) from None
            if zmax <= 0:
                msg = (f"Found negative zmax in {self.geometry} geometry "
                       f"(zmax = '{zmax}'). "
                       "zmax should be strictly positive.")
                raise ut.UserInputError(msg) from None
            if zmin >= zmax:
                msg = (f"Found zmin >= zmax in {self.geometry} geometry "
                       f"(zmin = '{zmin}', zmax = '{zmax}'). "
                       "zmin should be strictly smaller than zmax.")
                raise ut.UserInputError(msg) from None
        return z_init


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
