Advanced use
============

.. contents:: :local:

.. _stencil:

Stencil
-------

As shown in the :ref:`implementation_diffusion` Section, fluxes are discretized
using a resistance :math:`R_i` which represents a local average of
:math:`\Delta z/L`:

.. math::
   
   J^{\text{lat}}_i = - \frac{\mu_i - \mu_{i - 1}}{R_i}.

Three discretization schemes are provided with the parameter ``stencil``
(factory default: ``A``):

* ``A``: Quotient of arithmetic averages of :math:`\Delta z` and :math:`L`:
  
  .. math::
  
     R_i = \frac{\Delta z_i + \Delta z_{i - 1}}{L_i + L_{i - 1}}.

* ``H``: Arithmetic average of :math:`\Delta z/L`:

  .. math::
  
     R_i = \frac{1}{2}\left(\frac{\Delta z_i}{L_i} + \frac{\Delta z_{i - 1}}{L_{i - 1}}\right).

* ``G``: arithmetic average of :math:`\Delta z`, geometric average of :math:`L`:

  .. math::

     R_i = \frac{1}{2}\frac{\Delta z_i + \Delta z_{i - 1}}{\sqrt{L_iL_{i - 1}}}.

.. _time_step:

Time step
---------

The default time step is calculated as:

.. math::

   \Delta t_\mathrm{default} = \mathrm{Fo}
                               \frac{\Delta z^2_\mathrm{min}}{D^*_\mathrm{max}}

where :math:`\mathrm{Fo}` is the Fourier number (factory default: 0.4),
:math:`\Delta z_\mathrm{min}` is the minimum value of the initial
:math:`\Delta z` array and :math:`D^*_\mathrm{max}` is the maximum :math:`D^*`
value among all atom species over the initial concentration profile (see
:meth:`simu.Simulation.add_time` for implementation details). The actual time
time step is obtained by multiplying :math:`\Delta t_\mathrm{default}` by an
optional input parameter ``dt_multiplier``. This is usually used to reduce the
time step when the solver is unstable (``dt_multiplier`` < 1). Alternatively,
the time step may be modified with input parameter ``nt_multiplier`` (alias
``nt_mult``), which is used as the inverse of ``dt_multiplier``. One would
decrease the time step (i.e., increase the number of time steps) with
``nt_mult`` > 1.

Unstabilities mostly arise in two cases:

* When an initial or boundary atom fraction is close to 0. In this case,
  ``nt_multiplier`` values of 2-3 are usually sufficient to recover stability.

* When running simulations with finite sink strengths
  (see :ref:`non_ideal`).

.. _non_ideal:

Non-ideal lattice
-----------------

By default, simulations are run with the commonly used assumption of an ideal
lattice, where vacancies are maintained at equilibrium through the action of 
dislocations, which produces Kirkendall shift and no pores. Instead, users may
set a finite sink strength for dislocations, which will reduce the Kirkendall
shift and let vacancy supersaturation develop. A sink strength associated with
porosity may also be specified, which will produce Kirkendall porosity.

Sink strengths :math:`k_\mathrm{d}` and :math:`k_\mathrm{p}` as defined in
:ref:`lattice_velocity` can be specified in either of two ways:

* Using input parameters ``k_dislo`` and ``k_pores``, respectively. These may be
  either:

  * A numerical value.
  * A file name, which allows position-dependent values to be read from a file in
    the current job folder. The file can be in any format readable by the Numpy
    ``genfromtxt`` method and should contain `nz - 1` values.

* Or using input parameters ``rho_dislo`` and ``rho_pores``, which represent sink
  densities (in m\ :sup:`-2`). In this case :math:`k_i\ (i = d\ \text{or}\ p)`
  is obtained as :math:`k_i = \rho_i D_0^*`, where :math:`\rho_i` is the
  user-specified density and :math:`D_0^*` is the vacancy diffusion coefficient:
  
  .. math::
  
     D_0^* = \frac{1}{y_0}\sum_{k=1}^n{y_k\,D_k^*}

  :math:`D_0^*` is composition-dependent, therefore :math:`k_i` will be
  composition-dependent. Like ``k_dislo`` and ``k_pores``, ``rho_dislo`` and
  ``rho_pores`` can be scalar values or file names.

The lattice is considered non-ideal when either ``k_dislo``, ``k_pores``,
``rho_dislo`` or ``rho_pores`` is specified. If only one :math:`k_i` is given a
value (via ``k_xxxxx`` or ``rho_xxxxx``), the other takes the default value 0.

With a non-ideal lattice, the user can specify initial vacancy site fraction
and pore volume fraction profiles, using input parameters ``yVa_profile`` and
``fp_profile``, respectively. The syntax is the same as that applying to initial
atom fractions (see :ref:`initial_atom_fraction`), using keys ``yVa_left``,
``yVa_right``, ``fp_left``, ``fp_right``.

The :meth:`results.UnitResult.plot_quartet` method produces a 2 x 2 grid of
subplots with profiles of the following variables:

    * atom fraction;
    * relative difference between simulated and equilibrium vacancy site
      fraction;
    * flux in the lattice-fixed frame;
    * volume fraction of pores.

It is convenient to analyse the results of non-ideal lattice simulations.

.. note::

   Non-ideal lattice simulations require much smaller time steps than ideal lattice
   simulations: typically, ``nt_multiplier`` should be in the order of
   :math:`1/y_0^\mathrm{eq}`, which can be in the order of :math:`10^4 - 10^6`.
   These simulations take several hours to complete.

.. warning::
   Non-ideal lattice simulations may produce results that are difficult to
   interpret and should be run with caution. This type of simulation has not
   been thoroughly tested. Users are encouraged to be well acquainted with
   ideal lattice simulations and to read the :ref:`background` Section before
   working with non-ideal cases.

Run options
-----------

The :meth:`simu.Simulation.run` method has two optional arguments:

* ``show_completion`` (bool): whether to show times and completion rates,
  defaults to ``False``. Can be used to preview the duration of long simulations.
* ``verbose`` (int): verbosity level. Passed to the diffusion solver, controls
  the level of information logged (in particular, when volumes are created or
  deleted, see :func:`solvers.remesh`). Defaults to 1 (logs remeshing info).

.. _default_parameters:

Default parameters
------------------

When an optional parameter is not present in the input file, it is assigned a
default value. "Factory" default values are defined as part of the package
installation (see :data:`constants.factory_default_parameters`). These can be
overridden by creating a table named ``default_parameters`` in the user data file
(“user_data.toml”) and indicating the new default values as key-value pairs,
for example::

   # user_data.toml

   [default_parameters]

   number_space_steps = 100
   volume_db = 'Vegard'

In this case, when ``nz`` is absent from the input file, it will take the value
100. Likewise, ``volume_db`` will default to ``Vegard``.

In addition to optional parameters, :data:`constants.factory_default_parameters`
includes numerical parameters that cannot be set on a per-job basis, but can be
set system-wide via the ``default_parameters`` table in “user_data.toml”:

   * ``min_atom_fraction``: Noda solves the diffusion problem using transport
     coefficients which are proportional to concentrations (see Eq. 
     :eq:`Onsager_coefficients` in :ref:`mobility`). As a consequence, Noda
     cannot handle concentrations strictly equal to 0. Atom fractions given by
     the user to specify initial profiles and boundary conditions are therefore
     clipped between a lower bound, ``min_atom_fraction``, and a higher bound,
     ``1 - min_atom_fraction``.  The factory default is 1e-9.
   
   * ``min_number_time_steps``: this is a lower bound for the number of time
     steps. The factory default is 20.