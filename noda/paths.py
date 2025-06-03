# Copyright 2025 Onera
# This file is part of the Noda package
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate paths to useful directories."""

import os
from pathlib import Path

package_dir = Path(__file__).parent
pkg_data_dir = package_dir / 'data'


def get_data_dir(work_dir, logger):
    """
    Determine the appropriate data directory.

    Use either of these, based on the work directory:

    * if the job being run is part of the package test suite (validation tests
      and tutorials), use the data folder in the installation directory.
    * else (i.e. the job is a user simulation), if the 'NODA_HOME' has been
      defined and a 'data' folder is found in the user directory, use it; if
      not, fall back to the data folder in the installation directory.
    """
    if is_in_tests_dir(work_dir):
        data_dir = pkg_data_dir
        msg = f"Using data from package installation directory ({data_dir})."
        logger.info(msg)
    else:
        if 'NODA_HOME' in os.environ:
            data_dir = Path(os.environ['NODA_HOME']) / 'data'
            if data_dir.exists():
                msg = f"Using user-provided data in {data_dir}."
                logger.info(msg)
            else:
                msg = (f"The directory indicated by the environment variable "
                       f"'NODA_HOME', {data_dir}, was not found. Using "
                       "data from package installation directory instead "
                       f"({pkg_data_dir}).")
                logger.warning(msg)
                data_dir = pkg_data_dir

        else:
            data_dir = pkg_data_dir
            msg = ("The environment variable 'NODA_HOME' was not found. "
                   "Cannot use user-provided data. Using data from package "
                   "installation directory instead ({pkg_data_dir}).")
            logger.warning(msg)

    return data_dir

def is_in_tests_dir(job_dir):
    """Check if the job directory is part of the package test suite."""
    cond1 = job_dir.parents[0].name == 'jobs'
    cond2 = job_dir.parents[1].name == 'tests'
    cond3 = (job_dir.parents[1] / 'conftest.py').exists()
    return bool(cond1 and cond2 and cond3)
