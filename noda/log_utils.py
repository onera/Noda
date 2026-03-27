# Copyright 2025-2026 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Custom logger and associated utilities."""

import logging
import textwrap as tw
import re
import json

import numpy as np

import noda.utils as ut

# Define new levels and add to logging internal dict
logging.DATA = 15
logging.addLevelName(logging.DATA, "DATA")

logging.RESULTS = 16
logging.addLevelName(logging.RESULTS, "RESULTS")

logging.INPUT = 17
logging.addLevelName(logging.INPUT, "INPUT")

logging.INFO_NOSTREAM = 18
logging.addLevelName(logging.INFO_NOSTREAM, "INFO")

# Set formatting options
FORMAT = logging.Formatter('%(levelname)-7s: %(message)s')
MESSAGE_WIDTH = 70


class CustomLogger(logging.Logger):
    """
    Custom logger.

    Features:

    * Uses a custom DATA level.
    * Implements a stream handler (level INFO) and a file handler (level DATA).
    * Wraps messages before logging.

    """

    def __init__(self, sdir, ref, log=True):
        """Class constructor."""
        super().__init__(name=ref)
        self.setLevel(logging.DATA)
        if self.hasHandlers():
            self.handlers.clear()
        if log is True:
            self.add_stream_handler()
            self.add_file_handler(f"{sdir}/{ref}.nod")

    def add_stream_handler(self):
        """Add stream handler with level INFO to logger object."""
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        stream_handler.setFormatter(FORMAT)
        self.addHandler(stream_handler)

    def add_file_handler(self, file_path):
        """Add file handler with level DATA to logger object."""
        file_handler = logging.FileHandler(file_path, "w+")
        file_handler.setLevel(logging.DATA)
        file_handler.setFormatter(FORMAT)
        self.addHandler(file_handler)

    def data(self, msg, *args, **kwargs):
        """Log message with DATA level."""
        self.log(logging.DATA, msg, *args, **kwargs)

    def input(self, msg, *args, **kwargs):
        """Log message with INPUT level."""
        self.log(logging.INPUT, msg, *args, **kwargs)

    def results(self, msg, *args, **kwargs):
        """Log message with RESULTS level."""
        msg = {str(n): {var: d[var].tolist() for var in d}
               for n, d in msg.items()}
        self.log(logging.RESULTS, msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        """Log message with INFO level (overrides parent class method)."""
        stream = kwargs.pop('stream', True)
        if stream is True:
            self.log_wrapper(logging.INFO, msg, *args, **kwargs)
        else:
            self.log_wrapper(logging.INFO_NOSTREAM, msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        """Log message with WARNING level (overrides parent class method)."""
        self.log_wrapper(logging.WARNING, msg, *args, **kwargs)

    def log_wrapper(self, level, msg, *args, **kwargs):
        """Wrap message and log lines as multiple messages."""
        wrapped_lines = tw.wrap(msg.lstrip(), MESSAGE_WIDTH)
        for line in wrapped_lines:
            self.log(level, line, *args, **kwargs)


def parse_log(fpath):
    """
    Get simulation results from log file.

    Parameters
    ----------
    fpath: pathlib.Path instance
        Path to log file.

    Returns
    -------
    config : dict
        Input simulation parameters.
    results : dict
        Simulation results (see :meth:`simu.NewSimulation.run` and
        :func:`solvers.solver`).

    """
    with open(fpath, 'r', encoding='utf-8') as file:
        raw = file.read()
    config_str = re.findall(r"INPUT  : (.*)", raw)[0]
    config = json.loads(config_str.replace("'", '"'))
    try:
        res_str = re.findall(r"RESULTS: (.*)", raw)[0]
    except IndexError as exc:
        msg = "No results in output file. Run simulation first."
        raise ut.UserInputError(msg) from exc
    res_di = json.loads(res_str.replace("'", '"'))
    results = {int(n): {var: np.array(d[var]) for var in d}
               for n, d in res_di.items()}
    return config, results


def get_and_log(dct, key, default, logger, key_alias=None, stream=True):
    """
    Get dict value and send log if the key is absent from the dict.

    Parameters
    ----------
    dct : dictionary
        Dictionary to be probed.
    key : str
        Requested key.
    default : any
        Value returned if key is not in dict.
    key_alias : str, optional
        Name to be logged instead of key if key is absent. The default is None.

    Returns
    -------
    res : any
        Value obtained from dict or default.

    """
    try:
        res = dct[key]
    except KeyError:
        res = default
        if key_alias is not None:
            key = key_alias
        text = f"'{key}' absent from input file, using '{default}' as default."
        logger.info(text, stream=stream)
    return res
