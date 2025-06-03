# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define custom logger and logging wrapper function."""

import logging
import textwrap as tw

# Define new levels and add to logging internal dict
logging.DATA = 15
logging.addLevelName(logging.DATA, "DATA")

logging.RESULTS = 16
logging.addLevelName(logging.RESULTS, "RESULTS")

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

    def __init__(self, ref):
        """Class constructor."""
        super().__init__(name=ref)
        self.setLevel(logging.DATA)
        if self.hasHandlers():
            self.handlers.clear()
        self.add_stream_handler()
        self.add_file_handler(f'{ref}.nod')

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

    def data(self, message, *args, **kwargs):
        """Log message with DATA level."""
        self.log(logging.DATA, message, *args, **kwargs)

    def results(self, message, *args, **kwargs):
        """Log message with RESULTS level."""
        message = {str(n): {var: d[var].tolist() for var in d}
                   for n, d in message.items()}
        self.log(logging.RESULTS, message, *args, **kwargs)

    def info(self, message, *args, **kwargs):
        """Log message with INFO level (overrides parent class method)."""
        self.log_wrapper(logging.INFO, message, *args, **kwargs)

    def warning(self, message, *args, **kwargs):
        """Log message with WARNING level (overrides parent class method)."""
        self.log_wrapper(logging.WARNING, message, *args, **kwargs)

    def log_wrapper(self, level, message, *args, **kwargs):
        """Wrap message and log lines as multiple messages."""
        wrapped_lines = tw.wrap(message.lstrip(), MESSAGE_WIDTH)
        for line in wrapped_lines:
            self.log(level, line, *args, **kwargs)


# Logging wrappers


def get_and_log(di, key, default, logger, key_alias=None):
    """
    Get dict value and send log if the key is absent from the dict.

    Parameters
    ----------
    di : dictionary
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
        res = di[key]
    except KeyError:
        res = default
        if key_alias is not None:
            key = key_alias
        text = f"'{key}' absent from input file, using '{default}' as default."
        logger.info(text)
    return res
