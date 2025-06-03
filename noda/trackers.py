# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Decorators used to monitor functions."""

import time
import functools


def timer(func):
    """
    Timer function to be used as decorator.

    Prints run time of decorated function.
    Printing can be turned off when running the decorated function by passing
    the extra kwarg print_runtime to the function.
    The functools wrapper is used to retain the name and docstring of the
    wrapped function (useful for autodoc with sphinx).
    """
    @functools.wraps(func)
    def wrapper(*args, print_runtime=True, **kwargs):
        start_time = time.time()
        res = func(*args, **kwargs)
        dur = time.time() - start_time
        if print_runtime:
            print(f'"{func.__name__}" run time = {dur:.3f} s')
        return res
    return wrapper
