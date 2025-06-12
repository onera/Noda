Getting started
===============

.. contents:: :local:

File structure
--------------

The package root directory (or installation directory) is structured as follows::

    NODA
       ├───doc
       │   ├───_build
       │   ├───...
       │   └───source
       ├───noda
       │   │   alloy_system.py
       │   │   ...
       │   └───data
       │           fcc-AB-mob-bin_ideal.ods
       │           fcc-AB-thermo-bin_ideal.ods
       │           fcc-NiCrSi-mob-du2001.ods
       │           fcc-NiCrSi-thermo-schuster2000.ods
       │           user_data.toml
       └───tests
           │   conftest.py
           └───jobs
               ├───ideal_couple
               │       ideal_couple-input.txt
               │       ideal_couple.png
               │       run.py
               │       run_test.py
               ├───...

* The ``doc`` folder contains the source code of the documentation, and the html
  files once these are built.
* The ``noda`` folder contains the source code of the package and a set of data
  files used to run test jobs.
* The ``tests`` folder contains input files used to run test jobs.

The user directory, whose location is set by the ``NODA_HOME`` environment
variable (see :ref:`environment_variable`) is used to store user data --- it is
also typically used to store job folders, although these can be located anywhere
in the file system::

    NODA_user
       ├───data
       │       my_thermodynamics_database.ods
       │       my_mobility_database.xlsx
       │       user_data.toml
       └───jobs
            └───fancy_simulation
                    fancy1-input.txt
                    fancy2-input.txt
                    make_fancy_graphs.py

Running a simulation requires an input file, which must be a txt file named with
a reference followed by ``-input.txt`` (for example, ``fancy1-input.txt``).
Simulations can be run from a Python script (for example,
``make_fancy_graphs.py``) or interactively. The contents of the input file and
common use cases are described in :ref:`basic_use`.

When a simulation is run, if the ``NODA_HOME`` environment variable is provided
and a ``data`` folder is found in the user directory, Noda will look for data
there. If not, Noda will look for data in the installation directory ``noda/data``
folder.

.. note::

   It is recommended to not modify the files present in the installation directory
   ``noda/data`` folder. To get started, users may instead copy these files
   into the user directory, and modify the files and add new files there.

Tests and documentation
-----------------------

The ``tests`` folder in the installation directory serves several purposes. Test
jobs are used for validation during code development. It is also recommended to
run the test suite after installing to make sure the installation proceeded
properly (see :ref:`running_tests`). When test jobs are run, simulation results
are compared with reference results, including results from analytical
solutions in the cases where these exist. Running simulations requires data ---
this is why Noda ships with some example data. Jobs in the ``tests`` folder are
always run using the data in the ``noda/data`` folder. Test jobs are also used as
examples in the documentation: see :ref:`basic_use`. The data and job files
are a good starting point to prepare your own data and jobs.

As shown above, in addition to an input file and a running script, test job
folders contain a script named ``run_test.py``: this is used as part of the
automatic testing procedure, which relies on the Pytest framework
(https://docs.pytest.org). Some test job folders also contain image files used
in the documentation.

User data
---------

The file "user_data.toml" is used to store several types of data
required to run simulations. It uses the TOML format (see https://toml.io) and
contains four tables:

* ``molar_volume`` gathers partial molar volumes of pure metals, vacancies and
  pores.
* ``vacancy_formation_energy`` gathers vacancy formation energies in pure metals.
* ``thermodynamics`` stores names of thermodynamic data files.
* ``mobility`` does the same for mobility databases.

Each table contains one or more subtables, which define databases users can
choose from when setting up a simulation (:ref:`setting_up`).

For example, the "user_data.toml" file included in the installation directory
(``noda/data``) contains the following:

.. literalinclude:: /../../noda/data/user_data.toml
   :caption:

Two partial molar volume databases are present, named ``standard`` and
``Vegard``. When creating a simulation, Noda looks for the partial molar
volume of all system constituents in the selected database. If a constituent is
not present, it looks for the ``default`` key. If no ``default`` key is
provided in the selected database, it falls back to a system-wide default value
\ [#f1]_. Here, if using the ``standard`` database, all constituents will be
assigned the value given by the ``default`` key of this database, 1e-5 m3/mol;
if using the ``Vegard`` database, Al, Cr, Ni and Si will be given the indicated
partial molar volumes, and any other constituent will be given the system-wide
default.

.. note::

   Partial molar volumes of pure metals are numerical values. Vacancies and
   pores are special in that their partial molar volume can be a float or the
   string ``local``, in which case they take the local average molar volume ---
   see :ref:`background`.

The ``vacancy_formation_energy`` table works the same way as the
``molar_volume`` table, with a slightly different syntax: the crystal structure
of the pure metals is specified, and the energy values are given as
[enthalpy, entropy] lists (in [eV, eV/K]).

Subtables in the ``thermodynamics`` and ``mobility`` tables contain only one
key/value pair, used to indicate the name of the file containing the actual
data. The thermodynamics and mobility database files should be located in the
user data folder. Their content is described next.

.. _thermokin_database_files:

Thermodynamic and mobility database files
-----------------------------------------

Thermodynamic and mobility database files are spreadsheets in either ods or
xslx format.

The Gibbs free energy of the metal phase is described with the Calphad method,
using a Redlich-Kister polynomial for the excess term (see the :ref:`thermo`
Section). Thermodynamic database files contain two sheets:

* `Pure elements` stores coefficients that describe the temperature
  dependence of the Gibbs free energy of pure elements, in the form
  :math:`G - H_\mathrm{SER}`, according to [Dinsdale_1991]_:

  .. math::
   
     G = a + b\,T + c\,T\,\ln{T} + d_2\,T^2 + d_3\,T^3 + \frac{e}{T} + G_\mathrm{mag}

  (for the magnetic part :math:`G_\mathrm{mag}`, see [Dinsdale_1991]_).
* `Interactions` stores coefficients that describe the temperature
  dependence of the binary and ternary interaction parameters, in the
  form:
  
  .. math::
   
     _{}^{0}\Lambda_{ij} = A + B \cdot T,\\
     _{}^{1}\Lambda_{ij} = C + D \cdot T,\\
     \Lambda_{ijk} = A + B \cdot T.

.. note::

   Noda only considers ternary interactions of order 0. Entries for coefficients
   `C` and `D` of a ternary endmember will be ignored.

The logarithm of tracer diffusion coefficients is described with a
Redlich-Kister polynomial (see :ref:`mobility`). Mobility database files contain
one sheet, that stores coefficients that describe the temperature dependence of
unary terms, and binary and ternary interaction parameters:

.. math::
   
   \phi_\mathrm{solute}^\mathrm{solvent} = A + B \cdot T.

In this context, `solute` is the diffusing species, while `solvent` is the
endmember in which it diffuses.

The parameters that populate the thermodynamics and mobility database files
are typically found in journal articles. Users are recommended to use the
files provided in the ``tests`` folder of the installation directory and the
corresponding journal articles as examples to build their own database files.

.. rubric:: Footnotes

.. [#f1] A "factory" default value is defined in the package installation and
   can be overridden, see :ref:`default_parameters`.

.. rubric:: References

.. [Dinsdale_1991] A.T. Dinsdale, Calphad 15 (1991) 317-425,
   `DOI: 10.1016/0364-5916(91)90030-N <https://doi.org/10.1016/0364-5916(91)90030-N>`_
