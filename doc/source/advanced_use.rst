Advanced use
============

.. contents:: :local:

.. _L_mean_kind:

L_mean_kind
-----------

As shown in the :ref:`implementation_diffusion` Section, Onsager coefficients
(L) are known in the volumes (mid points), while fluxes are evaluated
between volumes (at node points), and therefore require L to be evaluated
at node points as well. The latter can be obtained via three kinds of means 
(factory default: ``arithmetic``):

* ``arithmetic``: :math:`\bar{L}_i = (L_i + L_{i - 1})/2`.
* ``harmonic``: :math:`\bar{L}_i = \frac{2}{1/L_i + 1/L_{i - 1}}`.
* ``geometric``: :math:`\bar{L}_i = \sqrt{L_iL_{i - 1}}`.

This parameter can be supplied in the ``options`` table of the configuration
file/dictionary (see examples in :ref:`couple_AB` and :ref:`source_AB`).

.. _time_step:

Time step
---------

The default time step is calculated as:

.. math::

   \Delta t_\mathrm{default} = \mathrm{Fo}
                               \frac{\Delta z^2_\mathrm{min}}{D^*_\mathrm{max}}

where :math:`\mathrm{Fo}` is the Fourier number (factory default: 0.25),
:math:`\Delta z_\mathrm{min}` is the minimum value of the initial
:math:`\Delta z` array and :math:`D^*_\mathrm{max}` is the maximum :math:`D^*`
value among all atom species over the initial concentration profile (see
:meth:`time_grid.TimeGrid.make_time_steps` for implementation details). The
actual time time step is obtained by multiplying :math:`\Delta t_\mathrm{default}`
by an optional input parameter ``dt_multiplier``. This is usually used to reduce
the time step when the solver is unstable (``dt_multiplier`` < 1). Alternatively,
the time step may be modified with input parameter ``nt_multiplier``, which is
used as the inverse of ``dt_multiplier``. One would decrease the time step
(i.e., increase the number of time steps) with ``nt_multiplier`` > 1.

Instabilities mostly arise in two cases:

* When an initial or boundary atom fraction is close to 0. In this case,
  ``nt_multiplier`` values of 2-3 are usually sufficient to recover stability.

* When running simulations with finite sink strengths
  (see :ref:`non_ideal`).

Parameters ``nt_multiplier`` or ``dt_multiplier`` can be supplied in the
``time`` table of the configuration file/dictionary.

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

Parameters ``k_dislo``, ``k_pores``, ``rho_dislo`` or ``rho_pores`` can be
supplied in the ``options`` table of the configuration file/dictionary.

With a non-ideal lattice, the user can specify initial vacancy site fraction
and pore volume fraction profiles, using tables ``vacancy_fraction`` and
``pore_fraction``, respectively, in the ``initial_conditions`` table of the
configuration file/dictionary. The syntax is the same as that applying to initial
atom fractions (see :ref:`initial_conditions`), for instance::

  config = {...
            'initial_conditions':
                {'atom_fractions' : {...},
                'vacancy_fraction':
                  {'shape': 'step',
                  'step_fraction': 0.5,
                  'left': 1e-6,
                  'right': 1e-2
                  }
               },
            ...
            }

The :meth:`results.UnitResult.plot_quartet` method produces a 2 x 2 grid of
subplots with profiles of the following variables:

    * atom fraction;
    * relative difference between simulated and equilibrium vacancy site
      fraction;
    * flux in the lattice-fixed frame;
    * volume fraction of pores.

It is convenient to analyze the results of non-ideal lattice simulations.

.. note::

   Non-ideal lattice simulations require much smaller time steps than ideal lattice
   simulations: typically, ``nt_multiplier`` should be in the order of
   :math:`1/y_0^\mathrm{eq}`, which can be in the order of :math:`10^4 - 10^6`.
   These simulations may take several hours to complete.

.. warning::
   Non-ideal lattice simulations may produce results that are difficult to
   interpret and should be run with caution. Users are encouraged to be well
   acquainted with ideal lattice simulations and to read the :ref:`background`
   Section before working with non-ideal cases.

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
installation (in :data:`constants.factory_default_parameters`). These can be
overridden by creating a table named ``default_parameters`` in the
“user_data.toml” file and indicating the new default values as key-value pairs,
for example::

   # user_data.toml

   [default_parameters]

   nz = 100
   num_out = 10

In this case, when the number of space steps ``nz`` is absent from the input
file, it will take the value 100. Likewise, the number of saved time steps
``num_out`` will default to 10.

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
  * ``partial_molar_volume`` : default value (m3/mol) assigned to a species
    when it is absent from a partial molar volume database, or when no partial
    molar volume database has been specified in the simulation configuration.
    The factory default is 1e-5.
  * ``vacancy_formation_energy`` : same with the vacancy formation energy in
    pure elements. The factory default is [2e5, 50] (GfV = HfV - T*SfV with
    [HfV, SfV] in [J/mol, J/mol/K]).
