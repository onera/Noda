.. _basic_use:

Basic use
=========

.. contents:: :local:

.. |br| raw:: html

     <br>

.. _setting_up:

Setting up a simulation
-----------------------

A typical input file contains the following:

.. literalinclude:: /../../tests/jobs/couple_NiCrSi/couple_NiCrSi-input.txt
   :caption:

Blank lines and code located after the comment character ``#`` are ignored when
the file is parsed. Parameters are to be entered with the following syntax: 
``key = value``, with any number of spaces around the ``=`` sign. Some keys
have an accepted alias (shorter name), as shown below.

Element symbols are case-insensitive --- for instance, ``fe`` is a valid symbol
for iron.

Some of the parameters are optional. When an optional parameter is not present
in the input file, a default value is assigned. "Factory" default values are
defined in the package installation but can be overridden by user-defined
values (see :ref:`default_parameters`).

Thermodynamic and mobility databases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``thermo_db``: name of the thermodynamic database to be used in the calculation.
* ``mob_db``: name of the mobility database to be used in the calculation.

These names must be associated with database files in ``user_data.toml``.

System constituents
^^^^^^^^^^^^^^^^^^^

* ``independent``: name independent constituents, with the following syntax:|br|
  ``element1, element2, ...``. Alias: ``inds``.
* ``dependent``: name of the dependent constituent. This is typically the alloy
  base element.|br|
  Alias: ``dep``.

System conditions
^^^^^^^^^^^^^^^^^

* ``phase``: name of the phase.
* ``TC``: temperature in Celsius.

.. _simulation_conditions:

Simulation conditions
^^^^^^^^^^^^^^^^^^^^^

* ``th``: simulation time in h.
* ``num_out``: number of saved time steps, optional (factory default: 2). A
  simulation often requires a large number of time steps; only the number of
  steps specified by ``num_out`` are saved to file, and accessible after the
  simulation is done.

.. note::

   The initial and last step are always saved, and these two steps are included 
   in ``num_out``. Therefore ``num_out`` must be greater than or equal to 2. For 
   instance, if the simulation time is 2 h and one wants to access the results 
   that correspond to 1 h of simulation, one may use |br|
   ``num_out = 3``. This will save results at 0, 1 and 2 h of simulation
   (``saved_times = [0, 1, 2]``).

.. _domain_geometry:

Domain geometry
^^^^^^^^^^^^^^^

Noda solves the diffusion problem along one space coordinate, in either of three
geometric configurations:

    * planar : a slab of given thickness in one direction
      and translational symmetry in the two other directions, described in
      Cartesian coordinates.
    * cylindrical : a cylinder with translational and rotational
      symmetry with respect to its axis, described in cylindrical coordinates.
    * spherical : a sphere with rotational symmetry with respect
      to its center, described in spherical coordinates.

The ``geometry`` parameter takes either argument ('planar', 'cylindrical' or
'spherical'). It is optional (factory default: 'planar').

.. note::

   The rotational symmetry in the cylindrical and spherical geometries has
   implications regarding the choice of the domain boundary positions and
   boundary conditions (see :ref:`initial_space_grid`).

.. _initial_space_grid:

Initial space grid
^^^^^^^^^^^^^^^^^^

* ``zmin``: position of left-hand domain boundary in m, optional (system
  default: 0).
* ``zmax``: position of right-hand domain boundary in m.
* ``nz``: number of space steps, optional (factory default: 51).
* ``grid``: optional (factory default: 'linear'). Can be either a standard grid
  type or a file name:

    * 'linear': linear grid from ``zmin`` to ``zmax`` with size ``nz``.
    * 'geo': geometric grid from ``zmin`` to ``zmax`` with size ``nz`` and
      common ratio ``q``:|br|:math:`\Delta z_{i+1} = q \Delta z_{i}`. The
      common ratio is indicated as an additional parameter, for example
      ``q = 1.02``. It is optional (factory default: 1.02).
    * filename: read grid positions from file in the current job folder. The
      file can be in any format readable by the Numpy ``genfromtxt`` method.
      ``zmin``, ``zmax`` and ``nz`` are inferred from the grid file, and must
      not be included in the input file.

.. note::

   As shown in :ref:`implementation_diffusion`, grid positions correspond to the
   edges of the volumes (nodes, noted ``z``). Fluxes are evaluated on nodes.
   Composition and related variables (concentrations, atom fractions, chemical
   potentials, molar volumes, ...) are evaluated in the center of the volumes,
   noted ``zm``. For instance, with ``grid = linear``, ``nz = 4`` and ``zmax = 3``
   will produce |br| ``z = [0, 1, 2, 3]`` and ``zm = [0.5, 1.5, 2.5]``.

.. note::

   At the moment, Noda only includes an explicit (forward Euler) time
   integration scheme, which is conditionnally stable. The minimum time step is
   based on the smallest space step (see :ref:`time_step`). Using a geometric
   grid with a large common ratio will produce locally small space steps and
   therefore require a large number of time steps.

.. note::

   In the cylindrical and spherical geometries:

   * The space coordinate is named  :math:`z`, instead of the usual :math:`\rho` 
     or :math:`r`, to favor compatibility across the three geometries.
   * The simulation domain is represented by positive radial coordinates, and
     the left-hand domain border is closest to the centre of symmetry.
     Therefore the following is enforced:
     :math:`0 \leq z_\mathrm{min} < z_\mathrm{max}`.
   * If ``zmin = 0``, the left-hand domain boundary is at the axis (or center)
     of symmetry. In this case, it is recommended to use a 0-flux boundary
     condition at this boundary.
   * It is possible to represent a hollow cylinder (or sphere) by giving a
     strictly positive value to ``zmin``.

.. _initial_atom_fraction:

Initial atom fraction profile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``x_profile``: can be either a standard profile type or a file name:
      
    * 'step' followed by a float argument: heaviside profile with step at
      ``zstep``. |br|
      If the argument is > 0.1, it is interpreted as a fraction of the domain
      length |br|
      :math:`\rightarrow` ``zstep = zmin + arg * (zmax - zmin)``. |br|
      If the argument is < 0.1, it is interpreted as ``zstep`` in m. |br|
      The atom fractions on the left and right sides of the step are read from
      the ``x_left`` and ``x_right`` parameters (both required).
    * 'smooth step' followed by a float: same with an error function instead of
      heaviside.
    * 'flat': flat profile. The atom fractions are read from the ``x_left``
      parameter (required).
    * filename: read values from file in the current job folder. The file can be
      in any format readable by the Numpy ``genfromtxt`` method. It must have
      the independent constituent profiles arranged by columns, with the
      constituent names on the first line. The size of the columns must match
      that of ``zm`` (i.e., ``nz - 1``). ``x_left`` and ``x_right`` must not be
      used.

* ``x_left``: name and atom fractions of each independent constituent on the
  left-hand side, with the following syntax: ``element1: value, element2: value, ...``.
* ``x_right``: same on the right-hand side.

Boundary conditions
^^^^^^^^^^^^^^^^^^^

Noda supports Neumann (fixed flux) and Dirichlet (fixed composition)
boundary conditions. For each of the left and right boundary, the user can
specify a composition (atom fractions of the independent constituents) or fluxes
(of all atom constituents). If no condition is given, a zero flux condition is
applied.

* ``JBC_left``: flux of atom constituents on left boundary, with the following
  syntax: |br|
  ``element1: value, element2: value, ...`` |br|
  The value is in SI units (mol/m2/s). It can also be an expression of time
  (noted ``t``) with basic Python operators, which will be evaluated with
  ``eval``, ex: ``(3*t + 2)**(1/2)``. If ``JBC_left`` is used, values must be
  given for all atom constituents.
* ``JBC_right``: same on right boundary.
* ``xBC_left``: atom fractions of independent constituents on left boundary,
  with the following syntax: ``element1: value, element2: value, ...`` |br|
  It can also be an expression of time (noted ``t``) with basic Python
  operators. If ``xBC_left`` is used, values must be given for all independent
  constituents.
* ``xBC_right``: same on right boundary.

.. note::

   Internally, Dirichlet boundary conditions are implemented by setting fluxes
   of the independent constituants in the lattice reference frame. By default,
   a flux of the dependent constituent is set so that the sum of the volumetric
   fluxes at the boundary is 0:
   
   .. math::
   
      J_n^\mathrm{lat} = -\frac{1}{V_n}\sum_{i=1}^{n - 1} J_i^\mathrm{lat}V_i,
   
   where the :math:`V_i` are the partial molar volumes. The volume of the
   simulation domain is conserved. The quantity of atom matter (in mol) may
   change due to differences in the partial molar volumes.

Optional parameters
^^^^^^^^^^^^^^^^^^^

The user may also specify:

* A partial molar volume database, with key ``volume_db``. The factory default is
  ``standard``.
* A vacancy formation energy database, with key ``vacancy_db``. The factory
  default is ``standard``.

Like the other optional parameters, the factory default values can be overridden
by user-specified default values (see :ref:`default_parameters`).

Running a diffusion simulation
------------------------------

A new simulation is created using the :class:`simu.NewSimulation` class::

   from noda import simu
   foo = simu.NewSimulation('couple_NiCrSi')
   
This supposes that an input file named ``couple_NiCrSi-input.txt`` is present
in the working directory. When the simulation is created, some information is
logged:

* Part of it (with the `INFO` level) is printed to screen and saved to the log
  file (here ``couple_NiCrSi.nod``): this indicates what default choices are
  made for all parameters not present in the input file, and gives general
  information as to how thermodynamics and mobility functions are generated.
* Another part (with the `DATA` level) is only saved to the log file: this
  contains data that will be used when loading the simulation.

The simulation is then run with the :meth:`simu.Simulation.run` method::

   foo.run()

Once the calculation is done, Noda logs the run time (screen and log file) and
the results (to the log file only).

Accessing and plotting simulation results
-----------------------------------------

Simulation results may be accessed either directly after a new simulation is
run (see above), or by loading the log file from a simulation run in a previous
session. In the latter case, the simulation is created using the
:class:`simu.ReadSimulation` class::

   foo = simu.ReadSimulation('couple_NiCrSi')
   
This supposes that a log file named ``couple_NiCrSi.nod`` is present in the
working directory.

The simulation results are stored at a number of time steps. The steps and the
associated simulation times (in h) can be accessed using the ``saved_steps``
and ``saved_times`` attributes of the simulation object:

>>> foo.saved_steps
array([  0,  47,  94, 141, 188, 235, 282, 329, 376, 423, 470])
>>> foo.saved_times
array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10.])

.. note::
   The unit of ``saved_times`` is hours. Here, there are 11 saved steps
   including the initial step (:math:`t=0`), and the simulated time is 10 h.
   However, the time step is not exactly a rational number multiplied by 3600.
   It follows that the saved times are not necessarily integer numbers of hours.

Other useful attributes are:

>>> foo.th       # time (h)
10.0
>>> foo.ts       # time (s)
36000.0
>>> foo.dt       # time step (s)
76.59574468085107
>>> foo.nt       # number of time steps
471
>>> foo.num_out  # number of saved time steps
11

Accessing results
^^^^^^^^^^^^^^^^^

The results are accessed using the :meth:`results.SimulationResults.result`
method. The simulation time one would like to access is specified with either
of two keyword arguments:

* ``step_index``, which is the index of the required step in ``saved_steps``.
  ``step_index`` must be an integer between 0 and the number of saved steps.
* ``th``, a time in h. In this case the simulation time accessed will be the
  closest to the specified value.

The results are also accessible as a dictionary of results objects. The
dictionary is called ``results`` and uses time steps as keys.

To wrap up, these are three ways to access the results after 4 h of simulation::

   bar = foo.result(step_index=4)
   bar = foo.result(th=4)
   bar = foo.results[188]

Both ``step_index`` and ``th`` are optional arguments of
:meth:`results.SimulationResults.result`. If no argument is given, the method
uses the default ``step_index = -1``, which is an alias for the last time step,
i.e., ``foo.result()`` returns the results at the last time step. Similarly,
``foo.results[-1]`` returns the results at the last time step.

The results objects (instances of :class:`results.UnitResult`) have attributes
which store the simulation variables as Numpy arrays (or dictionaries of Numpy
arrays). Commonly used attributes are:

* z : Node positions (m).
* zm : Midpoint positions (m).
* c : Concentrations (mol/m3).
* y : Metal site fractions.
* x : Metal atom fractions.
* Vm : Average molar volume of metal (m3/mol).
* mu : Chemical potentials (J/mol).
* Jlat : Fluxes in lattice-fixed frame (mol m-2 s-1).
* v : Velocity field of lattice relative to laboratory frame (m/s).
* deformation : Relative length variation.
    
Composition variables such as ``x`` or ``c`` are stored as dictionaries of 1D
arrays, with the relevant constituents as keys. The relevant constituents are
vacancies and all atom constituents, except for atom fractions ``x`` which only
apply to atom constituents:

>>> bar.x.keys()
dict_keys(['Cr', 'Si', 'Ni'])
>>> bar.c.keys()
dict_keys(['Va', 'Cr', 'Si', 'Ni'])
>>> bar.Jlat.keys()
dict_keys(['Va', 'Cr', 'Si', 'Ni'])

Composition variables are also accessible as 2D arrays of shape ``(nc, nz - 1)``
where ``nc`` is the number of constituents. These 2D arrays are named with the
suffix '_arr', for example ``x_arr``.

Plotting results
^^^^^^^^^^^^^^^^

Profiles of the simulated variables can be plotted using the
:meth:`results.UnitResult.plot` method. The variable is specified with the
optional argument ``varname``, which defaults to ``x``:

* ``bar.plot(varname='Jlat')`` will plot fluxes in the lattice reference frame.
* ``bar.plot()`` will plot atom fractions.

It is also possible to call plot directly from the simulation object; in this
case, the time step is specified with the optional arguments ``th`` or
``step_index``, which defaults to the last time step:

* ``foo.plot(varname='v', th=3)`` will plot the lattice velocity after 3 h.
* ``foo.plot(th=3)`` will plot atom fractions after 3 h.
* ``foo.plot(varname='mu')`` will plot chemical potentials at the last time step.
* ``foo.plot()`` will plot atom fractions at the last time step.

.. figure:: /../../tests/jobs/couple_NiCrSi/couple_NiCrSi.png
    :width: 400px
    :align: center

|

The :meth:`results.UnitResult.plot` method returns matplotlib figure and axis
objects. This allows modifying the graph settings after it is generated, for
example::

   fig, ax = foo.plot()
   ax.set_title('Custom title')

The :meth:`results.SimulationResults.interactive_plot` method allows accessing
time steps dynamically on a plot using a slider. Again the variable to be
plotted is specified with the optional argument ``varname``, which defaults to
``x``, and a shortcut is accessible from the simulation object:

* ``foo.interactive_plot(varname='mu')`` will plot chemical potentials.
* ``foo.interactive_plot()`` will plot atom fractions.

.. note::
   Interactive plots require an interactive graphics backend. If you are using
   Spyder, you will need to set the graphics backend to 'Automatic' rather than
   'Inline' (see Tools/Preferences/IPython console/Graphics/Graphics backend).

Accessing thermodynamic and diffusion properties
------------------------------------------------

Noda contains methods to compute the thermodynamic and diffusion properties of
a system, at given compositions. This can be done from a simulation object (an
instance of the :class:`simu.NewSimulation` class, see above), or by creating
a system object with the :class:`alloy_system.AlloySystem` class::

   from Noda import alloy_system
   foo = alloy_system.AlloySystem('NiCrSi')

The input file for :class:`alloy_system.AlloySystem` instances only requires
the following parameters:

.. literalinclude:: /../../tests/jobs/properties_NiCrSi/NiCrSi-input.txt
   :caption:

The following example illustrates how properties are computed:

.. literalinclude:: /../../tests/jobs/properties_NiCrSi/run.py
   :caption:

The methods shown here are also accessible from :class:`simu.NewSimulation`
instances.

Thermodynamic and diffusion properties are obtained as ndarrays and may be
plotted using standard matplotlib commands:

.. literalinclude:: /../../tests/jobs/properties_NiCrSi/plot.py
   :caption:

.. figure:: /../../tests/jobs/properties_NiCrSi/chemical_potentials_NiCrSi.png
    :width: 400px
    :align: center

|

The following example shows how to plot tracer diffusion coefficients and
intrinsic diffusion coefficients in a binary system. It reproduces Figure 3 in
[Gheno_2022]_:

.. literalinclude:: /../../tests/jobs/properties_NiCr/plot.py
   :caption:

.. figure:: /../../tests/jobs/properties_NiCr/diffusivities_NiCr.png
    :width: 400px
    :align: center

|

Working with model systems
--------------------------

Users can easily explore the effects of thermodynamic and diffusion properties
on the shape of concentration profiles by creating their own database files.
The examples below shows how the known solution to the diffusion problem
is recovered when the system properties follow a particular set of constraints.

Let us consider a single-phase binary system AB, where the average molar volume
is composition-independent. In 1D, the diffusion problem may be written in the
following form:

.. math::

   \frac{1}{V_\mathrm{m}}\frac{\partial x_B}{\partial t}
   = -\frac{\partial \tilde{J}_B}{\partial z}
   
   \tilde{J}_B = -\frac{\tilde{D}}{V_\mathrm{m}}\frac{\partial x_B}{\partial z}

This problem has analytical solutions when the interdiffusion coefficient is
constant --- in this case:

.. math::

   \frac{\partial x_B}{\partial t} = \tilde{D}\frac{\partial^2 x_B}{\partial z^2}
   
If A and B are both substitutional constituents of a disordered phase,
the interdiffusion coefficient is given by Darken's formula:

.. math::

   \tilde{D} = (x_A D_B^* + x_B D_A^*)\phi

where :math:`\phi` is the thermodynamic factor:

.. math::

   \phi = \frac{x_B}{RT}\frac{\partial \mu_B}{\partial x_B}

The simplest way for :math:`\tilde{D}` to be composition-independent is that the
system follows the following set of conditions:

* AB is an ideal solution: :math:`\mu_B = \mu_B^0 + RT\ln{x_B}`, which yields
  :math:`\frac{\partial \mu_B}{\partial x_B} = \frac{RT}{x_B}` and
  :math:`\phi = 1`;
* A and B have equal diffusivities: :math:`D_A^* = D_B^*`;
* The diffusivities are composition-independent: :math:`D_B^* \neq f(x_B)`.

This yields :math:`\tilde{D} = D_B^*`.

These conditions are fullfilled for the model AB system provided in the test
base, with databases named "thermo_bin_ideal" and "mob_bin_ideal".

A number of diffusion problems involving a binary system with
composition-independent molar volume and diffusivity have analytical solutions.
Some typical initial distributions and boundary conditions are illustrated
below.

Diffusion couple (planar geometry)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of an infinitely-thick planar diffusion couple:

.. math::

   &t = 0,\ z < z_\mathrm{step},\ x_B = x_B^\mathrm{left}\\
   &t = 0,\ z > z_\mathrm{step},\ x_B = x_B^\mathrm{right}\\
   &t > 0,\ \left.\frac{\partial x_B}{\partial z}\right|_{z\rightarrow -\infty}
   = \left.\frac{\partial x_B}{\partial z}\right|_{z\rightarrow +\infty} = 0,

the analytical solution to the diffusion problem is ([Crank_1975]_, p. 14):

.. math::

   \frac{x_B(z, t) - x_B^\mathrm{right}}{x_B^\mathrm{left} - x_B^\mathrm{right}}
   = \frac{1}{2} \mathrm{erfc}\left(\frac{z - z_\mathrm{step}}{2\sqrt{\tilde{D}t}}\right)

This is compared with the Noda simulation in the example named "couple_AB":

.. literalinclude:: /../../tests/jobs/couple_AB/couple_AB-input.txt
   :caption:

.. literalinclude:: /../../tests/jobs/couple_AB/run.py
   :caption:

.. figure:: /../../tests/jobs/couple_AB/couple_AB.png
    :width: 400px
    :align: center

|

Constant surface concentration (planar geometry)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case of a semi-infinite planar solid initially at a uniform concentration
:math:`x_B^\mathrm{bulk}`, with its left-hand surface maintained at a constant
concentration :math:`x_B^\mathrm{surf}`:

.. math::

   &t = 0,\ x_B = x_B^\mathrm{bulk}\\
   &t > 0,\ z = 0,\ x_B = x_B^\mathrm{surf},

the solution reads ([Crank_1975]_, p. 32):

.. math::

   \frac{x_B(z, t) - x_B^\mathrm{bulk}}{x_B^\mathrm{surf} - x_B^\mathrm{bulk}}
   = \mathrm{erfc}\left(\frac{z}{2\sqrt{\tilde{D}t}}\right)

This is compared with the Noda simulation in the example named
"source_AB", which illustrates the use of a geometric grid, well suited to
the Dirichlet boundary condition:

.. literalinclude:: /../../tests/jobs/source_AB/source_AB-input.txt
   :caption:

.. literalinclude:: /../../tests/jobs/source_AB/run.py
   :caption:

.. figure:: /../../tests/jobs/source_AB/source_AB.png
    :width: 400px
    :align: center

|

Constant surface concentration (spherical geometry)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a sphere of radius :math:`R` is initially at a uniform concentration
:math:`x_B^\mathrm{bulk}`, and its surface is maintained at a constant
concentration :math:`x_B^\mathrm{surf}`, the concentration at distance
:math:`r` from the center is ([Crank_1975]_, p. 91):

.. math::

   \frac{x_B(r, t) - x_B^\mathrm{bulk}}{x_B^\mathrm{surf} - x_B^\mathrm{bulk}}
   = 1 + \frac{2R}{\pi r} \sum_{n=1}^{\infty} \frac{(-1)^n}{n}
   \sin{\left(\frac{n\pi r}{R}\right)} \exp\left(\frac{-Dn^2\pi^2t}{R^2}\right)

This is compared with the Noda simulation after different diffusion times in
the example named "sphere_AB":

.. literalinclude:: /../../tests/jobs/sphere_AB/sphere_AB-input.txt
   :caption:

.. literalinclude:: /../../tests/jobs/sphere_AB/run.py
   :caption:

.. figure:: /../../tests/jobs/sphere_AB/sphere_AB.png
    :width: 400px
    :align: center

|

.. rubric:: References

.. [Crank_1975] J. Crank, The mathematics of diffusion, 2nd ed., Oxford University Press,
   Oxford, 1975

