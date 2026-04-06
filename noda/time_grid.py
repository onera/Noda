# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define time grid."""

import numpy as np

import noda.utils as ut
from noda.log_utils import get_and_log


class TimeGrid:
    """
    Store time-related parameters.

    Attributes
    ----------
    th : float
        Simulation time in h.
    ts : float
        Simulation time in s.
    dt_multiplier : float
        Factor by which default time step is multiplied (see
        :meth:`make_time_steps`).
    dt : float
        Time step size in s.
    nt : int
        Number of time steps (positions in time sequence), including 0.
    num_out : int
        Number of time steps where simulation results are stored and saved on
        file.
    saved_steps : list of int
        Steps for which simulation results are stored and saved in log.
    saved_th : list of float
        Times in h for which simulation results are stored and saved in log.
    default_params : dict
        Parameters used when not specified in user input.
    logger : :class:`log_utils.CustomLogger`
        Logger.

    """

    def __init__(self, params, x, dz, DT_funx, default_parameters, logger):
        """
        Class constructor.

        The time step dt is calculated to ensure the stability of explicit
        schemes, with Fourier given in default_parameters. The DT value is
        calculated from the initial atom fraction profile. The time step
        can be made smaller using the user-specified dt_multiplier or
        nt_multiplier.

        Parameters
        ----------
        params : dict
            Time-related input parameters.
        x : dict of 1D array
            Initial atom fraction profile.
        dz : 1D array
            Initial space step.
        DT_funx : function
            Calculates tracer diffusion coefficients.
        default_parameters : dict
            Parameters used when not specified in user input.
        logger : :class:`log_utils.CustomLogger`
            Logger.

        Raises
        ------
        :class:`utils.UserInputError`
            If params contains both th and ts.
        :class:`utils.UserInputError`
            If params contains both nt_multipler and dt_multiplier.

        """
        if 'th' in params and 'ts' in params:
            msg = "Cannot specify both th and ts."
            raise ut.UserInputError(msg) from None
        if 'th' in params:
            self.th = params['th']
            self.ts = self.th * 3600
        elif 'ts' in params:
            self.ts = params['ts']
            self.th = self.ts / 3600
        if 'dt_multiplier' in params and 'nt_multiplier' in params:
            msg = "Cannot specify both dt_multiplier and nt_multiplier."
            raise ut.UserInputError(msg) from None
        if 'nt_multiplier' in params:
            nt_multiplier = params['nt_multiplier']
            self.dt_multiplier = 1/nt_multiplier
        else:
            self.dt_multiplier = params.get('dt_multiplier', 1)

        self.default_parameters = default_parameters
        self.logger = logger
        self.dt, self.nt = self.make_time_steps(x, dz, DT_funx)
        self.saved_steps, self.saved_th = self.make_saved_steps(params)
        self.num_out = self.saved_steps.size

    def make_time_steps(self, x, dz, DT_funx):
        """Calculate time steps."""
        x_arr = np.array(list(x.values()))
        DTmax = np.max(DT_funx(x_arr))
        Fo = self.default_parameters['Fourier_number']
        dt = Fo*self.dt_multiplier*dz.min()**2/DTmax
        nt = int(self.ts/dt) + 1
        # Adjust dt to recover ts (rounding nt introduces error)
        dt = self.ts/(nt - 1)

        # Impose lower bound on number of time steps
        if nt < self.default_parameters['min_number_time_steps']:
            nt = self.default_parameters['min_number_time_steps']
            dt = self.ts/(nt - 1)
        return dt, nt

    def make_saved_steps(self, params):
        """
        Add list of time steps where simulation results will be saved to file.

        The time steps are determined based on the num_out parameter.
        """
        if 'num_out' in params:
            if params['num_out'] < 2:
                msg = ("'num_out' must be greater than 2. "
                       f"Found 'num_out = {params['num_out']}'")
                raise ut.UserInputError(msg) from None
        num_out = get_and_log(params, 'num_out',
                              self.default_parameters['num_out'],
                              self.logger)
        if num_out == 'all':
            num_out = self.nt
        steps = np.linspace(0, self.nt - 1, num=num_out)
        saved_steps = np.around(steps).astype(int)
        saved_th = saved_steps*self.dt/3600

        if len(saved_steps) > 1e4:
            question = """
            The simulation will generate a large output file. Do you want to
            proceed ? [y/N]
            (Check "num_out" parameter in input file.)"""
            ans = input(question) or 'n'
            if ans != ('y' or 'Y'):
                raise ut.UserInputError('Canceled')
        return saved_steps, saved_th
