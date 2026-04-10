Getting started
===============

.. contents:: :local:

File structure
--------------

The package installation directory is structured as follows::

    NODA
       ├───doc
       │   ├───_build
       │   ├───...
       │   └───source
       ├───noda
       │   │   alloy_system.py
       │   │   ...
       │   └───data
       │           fcc-AB-mob-ideal.ods
       │           fcc-AB-thermo-ideal.ods
       │           fcc-NiCrSi-mob-Du2001.ods
       │           fcc-NiCrSi-thermo-Schuster2000.ods
       │           user_data.toml
       └───tests
           │   conftest.py
           └───jobs
               ├───ideal_couple
               │       ideal_couple.toml
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
       │       my_mobility_database.ods
       │       user_data.toml
       └───jobs
            └───fancy_simulation
                    simulation1.toml
                    simulation2.toml
                    make_fancy_graphs.py

Simulations are created from a Python script (for example,
``make_fancy_graphs.py``) or interactively in a Python interpreter. The
configuration (system properties, simulation conditions) is provided
through an input file with the `TOML <https://toml.io>`__ format (for example,
``simulation1.toml``) or directly through a Python dictionary. The
contents of the configuration and common use cases are described in
:doc:`basic_use`.

When a simulation is run, if the ``NODA_HOME`` environment variable is provided
and a ``data`` folder is found in the user directory, Noda will look for data
there. If not, Noda will look for data in the installation directory ``noda/data``
folder.

.. note::

   Users should not modify the files present in the installation directory
   ``noda/data`` folder. To get started, users may instead copy these files
   to the user directory, and modify them there.

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
examples in the documentation: see :doc:`basic_use`. The data and job files
are a good starting point to prepare your own data and jobs.

As shown above, in addition to an input file and a running script, test job
folders contain a script named ``run_test.py``: this is used as part of the
automatic testing procedure, which relies on the Pytest framework
(https://docs.pytest.org). Some test job folders also contain image files used
in the documentation.

.. _user_data:

User data
---------

The file "user_data.toml" is a database register. It uses the
`TOML <https://toml.io>`__ format and contains five tables:

* ``[thermodynamics]`` stores thermodynamic database file names.
* ``[mobility]`` does the same for mobility databases.
* ``[partial_molar_volume]`` stores databases of partial molar volumes of pure
  elements, vacancies and pores.
* ``[vacancy_formation_energy]`` stores databases of vacancy formation energies
  in pure elements.
* ``[default_parameters]`` stores default values for various parameters.

For example, the "user_data.toml" file included in the installation directory
(``noda/data``) contains the following:

.. literalinclude:: /../../noda/data/user_data.toml
   :caption:

The ``[thermodynamics]`` and ``[mobility]`` tables associate database names with
the names of the files containing the actual data. These files
should be located in the user data folder. When setting up a simulation, users
may provide one of these database names or directly indicate a file path (see
:ref:`setting_up`). The content of the database files is described next
(:ref:`thermokin_database_files`).

The ``[molar_volume]`` and ``[vacancy_formation_energy]`` tables associate
database names with partial molar volume and vacancy formation energy data.
When setting up a simulation, users may provide one of these database names, or
directly provide data (see :ref:`setting_up`). The program looks
for the partial molar volume of all components in the provided database. If a
component is not present, it falls back to a system-wide default value. Here,
if using the ``Vegard`` database, Al, Cr, Ni and Si will be
assigned the indicated partial molar volumes, and any other constituent will be
assigned the system-wide default. The ``vacancy_formation_energy`` table works
the same way as the ``molar_volume`` table, with a slightly different syntax:
the crystal structure of the pure metals is specified, and the energy values
are given as [H, S] lists (in [J/mol, J/mol/K]). The vacancy
formation energy in element `A` is calculated as
:math:`G_\mathrm{f}^\mathrm{A} = H_\mathrm{f}^\mathrm{A} - TS_\mathrm{f}^\mathrm{A}`
. The values are then used to calculate chemical potentials, see :ref:`thermo`.

.. note::

   Partial molar volumes of pure metals are numerical values. Vacancies and
   pores are special in that their partial molar volume can be a float or the
   string ``local``, in which case they take the local average molar volume ---
   see :ref:`background`.

Factory default values for a number of parameters are provided in the package
files (in the :mod:`constants` module). These can be overridden by specifying
values in the ``[default_parameters]`` table, see :ref:`default_parameters`.

.. _thermokin_database_files:

Thermodynamic and mobility database files
-----------------------------------------

Thermodynamic and mobility database files are spreadsheets in csv, ods or
xlsx format.

The Gibbs free energy of the metal phase is described with the Calphad method,
using a Redlich-Kister polynomial for the excess term (see the :ref:`thermo`
Section). Thermodynamic database files in ods and xlsx format contain two sheets:

* `Elements` stores coefficients that describe the temperature
  dependence of the Gibbs free energy of pure elements, in the form
  :math:`G - H_\mathrm{SER}`, according to [Dinsdale_1991]_:

  .. math::
   
     G = a + b\,T + c\,T\,\ln{T} + d_2\,T^2 + d_3\,T^3 + \frac{e}{T} + G_\mathrm{mag}

  (for the magnetic part :math:`G_\mathrm{mag}`, see [Dinsdale_1991]_).
* `Interactions` stores coefficients that describe the temperature
  dependence of the binary and ternary interaction parameters, in the
  form:
  
  .. math::
   
     _{}^{0}\Lambda_{ij} = a + b \cdot T,\\
     _{}^{1}\Lambda_{ij} = c + d \cdot T,\\
     \Lambda_{ijk} = a + b \cdot T.

.. note::

   Noda only considers ternary interactions of order 0. Entries for coefficients
   `c` and `d` of a ternary endmember will be ignored.

Thermodynamic database files in csv format store the parameters under
`Elements` and `Interactions` headings.

The logarithm of tracer diffusion coefficients is described with a
Redlich-Kister polynomial (see :ref:`mobility`). Mobility database files contain
coefficients that describe the temperature dependence of
unary terms, and binary and ternary interaction parameters:

.. math::
   
   \phi_\mathrm{solute}^\mathrm{solvent} = a + b \cdot T.

In this context, `solute` is the diffusing species, while `solvent` is the
endmember in which it diffuses, or the subsystem that generates interactions.
For instance, in a Ni-Cr-Al system, there are three solutes (Ni, Cr, Al),
three endmembers (Ni, Cr, Al), three binary subsystems (NiCr, NiAl, CrAl) and
one ternary subsystem (NiCrAl).

The parameters that populate the thermodynamics and mobility database files
are typically found in journal articles. Users are recommended to use the
files provided in the ``noda/data`` folder of the installation directory as
examples to build their own database files.

.. rubric:: References

.. [Dinsdale_1991] A.T. Dinsdale, Calphad 15 (1991) 317-425,
   `DOI: 10.1016/0364-5916(91)90030-N <https://doi.org/10.1016/0364-5916(91)90030-N>`_
