.. _background:

Background
==========

.. contents:: :local:

Models
------

This section describes the models used in Noda. Additionnal information is
available in [Gheno_2022]_.

.. note::

   Noda has been developped in part to study the porosity stemming from the
   Kirkendall effect. This requires explicitely solving a particular
   (non-conserved) continuity equation for vacancies, with *finite sink strengths*,
   and considering a pore phase separate from the metal. In many cases, one is
   not interested in the Kirkendall porosity, and one wants to solve the diffusion
   equation for the atom species, assuming vacancies remain at equilibrium, and
   no pore forms. This is the default behavior. Vacancies may depart from
   equilibrium (and porosity develop) only if sink strengths are specified in the
   input file, see :ref:`non_ideal`.

.. _composition_variables:

Composition variables
^^^^^^^^^^^^^^^^^^^^^

This version of Noda is restricted to systems with one metal phase and (possibly)
one pore phase. Let :math:`W` be a representative volume comprising some metal
and some pores:

.. math::
   :label: W_sum
   
   W = W_\mathrm{m} + W_\mathrm{p}.
   
Within :math:`W`, metal and pores are not spatially resolved; we consider the
average properties of the volume.

The metal is a cristalline, disordered solid solution, with a single sublattice
(interstitial species may be represented as substitutional species with a
partial molar volume equal to zero). It contains :math:`N` mol of lattice sites,
occupied by atoms and vacancies. The pore contains no site. The overall system
composition is described using the variables gathered in
:numref:`Table %s <compvar>`, where index 0 represents vacancies and indices
1 to :math:`n` represent atoms.

.. _compvar:

.. table:: System-wide composition variables

   ============== ================================= ============================================ ================= =============================================
   Variable       Definition                        Species                                      Unit              Closure
   ============== ================================= ============================================ ================= =============================================
   Amount         :math:`N_{k}`                     :math:`k \in \left\lbrack 0,n \right\rbrack` mol               :math:`\sum_{k = 0}^{n}N_{k} = N`
   Site fraction  :math:`y_k = \frac{N_{k}}{N}`     :math:`k \in \left\lbrack 0,n \right\rbrack` \-                :math:`\sum_{k = 0}^{n}y_k = 1`
   Atom fraction  :math:`x_k = \frac{y_k}{1 - y_0}` :math:`k \in \left\lbrack 1,n \right\rbrack` \-                :math:`\sum_{k = 1}^{n}x_k = 1`
   Concentration  :math:`c_k = \frac{y_k}{\Omega}`  :math:`k \in \left\lbrack 0,n \right\rbrack` mol/m\ :sup:`3`   :math:`\sum_{k = 0}^{n}c_k = \frac{1}{\Omega}`
   ============== ================================= ============================================ ================= =============================================

Concentrations in the metal phase are defined as:

.. math::
   :label: conc_m
   
   c_k^\mathrm{m} = \frac{y_k}{\Omega_\mathrm{m}},

where :math:`\Omega_\mathrm{m}` is the molar volume of the metal. The latter is
defined as:

.. math::
   :label: Omega_my
   
   \Omega_\mathrm{m} = \sum_{k = 0}^{n}{y_k\Omega_k\ },

where :math:`\Omega_k` is the partial molar volume of species :math:`k` in the
metal phase, assumed to be constant. The partial molar volume of vacancies,
:math:`\Omega_0`, can be assigned either an independent value, or the local
metal molar volume (see [Svoboda_2006]_). In the latter case,

.. math::
   :label: Omega_mx
   
   \Omega_\mathrm{m} = \sum_{k = 1}^{n}{x_k\Omega_k}.

The volume occupied by the metal is

.. math::
   :label: W_m

    W_\mathrm{m} = N \Omega_\mathrm{m}.

The global molar volume :math:`\Omega` is defined such that:

.. math::
   :label: Omega

   W = N \Omega.
   
Volume fractions are defined as

.. math::
   :label: volume_fractions

   f_{i} = \frac{W_{i}}{W},\ \ \ i = \mathrm{m\ or\ p}.
   
The global and metal molar volumes are related through:

.. math::
   :label: Omega_Omega
   
   \Omega = \frac{\Omega_\mathrm{m}}{f_\mathrm{m}}.

Concentrations are related through:

.. math::
   :label: conc_conc
   
   c_k = c_k^\mathrm{m} f_\mathrm{m}.

.. _thermo:

Thermodynamics of the metal phase
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermodynamic properties of the metal phase are modeled according to the
Calphad method. The Gibbs free energy per mole of lattice site is written
[Saunders_1998]_ [Lukas_2007]_:

.. math::
   :label: G_M
   
   G_\mathrm{M} = \sum_{i = 0}^{n}{y_i G_{i}} + RT\sum_{i = 0}^{n}{y_i \ln y_i} +_{}^{\mathrm{xs}}G_\mathrm{M},

where :math:`G_{i}` is the Gibbs energy of endmember :math:`i`, :math:`R`
and :math:`T` have their usual meaning and :math:`_{}^{\mathrm{xs}}G_\mathrm{M}`
is the excess term. The latter is modeled using a Redlich-Kister
polynomial, with binary interactions of order 1 and, in systems with 3 or more
elements, ternary interactions of order 0:

.. math::
   :label: G_xs
   
   _{}^{\mathrm{xs}}G_\mathrm{M} = \sum_{i = 0}^{n - 1}{\sum_{j = i + 1}^{n}{y_i y_j\sum_{\nu = 0}^{1}{_{}^{\nu}\Lambda_{ij}}\left(y_i - y_j\right)^{\nu}}}
                                   + \sum_{i = 0}^{n - 2}{\sum_{j = i + 1}^{n - 1}{\sum_{k = j + 1}^{n}{y_i y_j y_k\Lambda_{ijk}}}}.

The binary and ternary interaction terms, :math:`_{}^{\nu}\Lambda_{ij}`
and :math:`\Lambda_{ijk}`, are expressed in the form :math:`A + B \cdot T`.

Chemical potentials are defined as:

.. math::
   :label: mu
   
   \mu_k = \frac{\partial G}{\partial N_{k}},

where :math:`G` is the total Gibbs free energy, :math:`G = N \cdot G_\mathrm{M}`.

While the Gibbs energies of metals can readily be calculated from parameters
in the literature, the Gibbs energy of the vacancy has no trivial definition.
A common practice is to choose a value arbitrarily to ensure the existence and
uniqueness of an equilibrium state between the alloy and a vacancy phase
(a pore). Here we follow the recommendation of [Franke_2014]_ and use
:math:`G_{\mathrm{Va}} = RT(\ln 2\ –\ 1/2)`. The value of
:math:`G_{\mathrm{Va}}` has no significant impact on the simulation
results.

Likewise, atom-vacancy interaction parameters are not generally known.
Following [Abe_2018]_, we assume 0-order binary interactions (i.e.,
:math:`_{}^{1}\Lambda_{0i} = 0\ \forall\ i` and
:math:`\Lambda_{0ij} = 0\ \forall\ (i,j)`), and
the interaction parameters are determined based on the vacancy formation
energy in pure metals

.. math::
   :label: vacancy_interaction_parameters
   
   \Lambda_{0k} = G_\mathrm{f, Va}^{k} - G_{\mathrm{Va}}.

The equilibrium vacancy fraction is obtained by solving :math:`\mu_0 = 0`.
This yields an expression for :math:`y_0^{\mathrm{eq}}`, which can
be written:

.. math::
   :label: y0_eq
   
   y_0^{\mathrm{eq}} = \exp{\left( - \frac{G_\mathrm{f, Va}^{\mathrm{alloy}}}{RT} \right)},

where :math:`G_\mathrm{f, Va}^{\mathrm{alloy}}` is the vacancy formation energy in the
alloy. The analytical form of :math:`G_\mathrm{f, Va}^{\mathrm{alloy}}` depends on the
model chosen for :math:`_{}^{\mathrm{xs}}G_\mathrm{M}`. The chemical potential of
vacancies can be written in the form:

.. math::
   :label: mu_0
   
   \mu_0 = RT \cdot \ln\left( \frac{y_0}{y_0^{\mathrm{eq}}} \right).

.. _mobility:

Mobility
^^^^^^^^

Flux densities are described using two reference frames: (i) the lattice
reference frame, based on a local coordinate system attached to inert
markers in the lattice, and (ii) the laboratory reference frame, based
on a global coordinate system that coincides with the lattice reference
frame at :math:`t = 0` but does not deform thereafter. The lattice moves
with respect to the laboratory reference frame with a velocity field
:math:`v`. The fluxes of species :math:`k` in the lattice and
laboratory reference frames, :math:`J_k^{\mathrm{lat}}` and
:math:`J_k^{\mathrm{lab}}`, are related by:

.. math::
   :label: change_frame
   
   J_k^{\mathrm{lab}} = J_k^{\mathrm{lat}} + c_k \cdot v\ ,

where we recognize diffusive (:math:`J_k^{\mathrm{lat}}`) and advective
(:math:`c_k \cdot v`) contributions. :math:`J_k^{\mathrm{lat}}` is written
[Philibert_1991]_:

.. math::
   :label: Jlat_full
   
   J_k^{\mathrm{lat}} = - \sum_{i = 0}^{n}{L_{ki}\ \nabla\mu_i\ },

where the :math:`L_{ki}` are the Onsager coefficients, which form a
symmetric matrix. Lattice sites are conserved by substitutional diffusion:

.. math::
   :label: Jlat_closure
   
   \sum_{k = 0}^{n}J_k^{\mathrm{lat}} = 0.

This allows the vacancy flux to be expressed as a function of :math:`n`
independent atom fluxes, and all :math:`L_{k0}` to be expressed as a function
of the atom transport coefficients, which are related to measurable mobilities.
After some manipulation one obtains:

.. math::
   :label: Jlat_novac
   
   J_k^{\mathrm{lat}} = - \sum_{i = 1}^{n}{L_{ki}\ \nabla{\tilde{\mu}}_i}\ ,

where :math:`{\tilde{\mu}}_i` is the diffusion potential,
defined as :math:`{\tilde{\mu}}_i = \mu_i - \mu_0`.

The :math:`L_{ki}` reflect diffusional transport in the {metal + pore} system,
and we assume these can be expressed as a weighted average of transport
coefficients in the metal and in the pore:

.. math::
   :label: Lki_full

   L_{ki} = f_\mathrm{m}L_{ki}^\mathrm{m} + f_\mathrm{p}L_{ki}^\mathrm{p}.

We further assume :math:`L_{ki}^\mathrm{p} = 0\ \ \forall(k,i)`, i.e.,
no diffusion in the pores. This yields:

.. math::
   :label: Lki

   L_{ki} = f_\mathrm{m}L_{ki}^\mathrm{m}

where :math:`L_{ki}^\mathrm{m}` are the transport coefficients in the metal.

Following Ågren [Agren_1982]_, we consider the vacancy-exchange
mechanism, in which diffusion occurs by atoms jumping into neighboring
vacant sites, and vacant sites are distributed randomly. The transport
coefficients are expressed as:

.. math::
   :label: Onsager_coefficients
   
   \left\{
   \begin{array}{l}
   L_{kk}^\mathrm{m} = c_k^\mathrm{m}y_0 M_{k0}\\
   L_{ki}^\mathrm{m} = 0\quad\quad\quad\quad \mathrm{for}\ k \neq i,
   \end{array}
   \right.

where :math:`M_{k0}` is a kinetic parameter representing the rate of exchange
between a :math:`k` atom and a neighboring vacancy. The mobility of species
:math:`k` is [Andersson_1992]_:

.. math::
   :label: mobility
   
   M_k = y_0 M_{k0}.

The mobility is related to the tracer diffusion coefficient :math:`D_k^*`
by Einstein relation:

.. math::
   :label: Einstein
   
   M_k = \frac{D_k^*}{RT}.

Mobilities evaluated from tracer diffusion coefficients are considered
equilibrium mobilities, i.e.,

.. math::
   :label: Einstein_equilibrium
   
   \frac{D_k^*}{RT} = M_k^{\mathrm{eq}} = y_0^{\mathrm{eq}}M_{k0}.

Consequently, the mobility can be written with the general form:

.. math::
   :label: mobility2
   
   M_k = \frac{y_0}{y_0^{\mathrm{eq}}}M_k^{\mathrm{eq}}.

Using Eqs. :eq:`Jlat_novac`-:eq:`mobility2`, the fluxes are finally written:

.. math::
   :label: Jlat_reduced
   
   J_k^{\mathrm{lat}} = - L_{kk}\ \nabla {\tilde{\mu}}_{k},

with

.. math::
   :label: Lkk
   
   L_{kk} &= f_\mathrm{m}c_k^\mathrm{m}\frac{y_0}{y_0^{\mathrm{eq}}}\frac{D_k^*}{RT}\\
          &= c_k\frac{y_0}{y_0^{\mathrm{eq}}}\frac{D_k^*}{RT}.

The temperature dependence of tracer diffusion coefficients is described with
the Arrhenius relation:

.. math::
   :label: Arrhenius
   
   D_k^* = D_k^{0}\exp{\left( - \frac{Q_{k}}{RT} \right)\ }

where :math:`D_k^{0}` is the preexponential factor and :math:`Q_{k}` the
activation energy. Following a common practice in the mobility literature, Eq.
:eq:`Arrhenius` is expressed:

.. math::
   :label: lnDT
   
   \ln D_k^* = \frac{\phi_k}{RT},

with :math:`\phi_k = RT\ln D_k^{0} - Q_{k}`. The composition dependence is 
then given by expanding :math:`\phi_k` with a Redlich-Kister polynomial, with
binary and ternary interactions of order 0:

.. math::
   :label: lnDT_RK
   
   \phi_k = \sum_{i = 1}^{n}{x_i \phi_k^i}
            + \sum_{i = 1}^{n - 1}{\sum_{j = i + 1}^{n}{x_i x_j\Lambda_{ij}}}
            + \sum_{i = 1}^{n - 2}{\sum_{j = i + 1}^{n - 1}{\sum_{k = j + 1}^{n}{x_i x_j x_k\Lambda_{ijk}}}}.

The binary and ternary interaction terms, :math:`\Lambda_{ij}` and
:math:`\Lambda_{ijk}`, are expressed in the form :math:`A + B \cdot T`.

Diffusion
^^^^^^^^^

Mass balance
""""""""""""

Following [Svoboda_2006]_, the mass balance in the lattice reference frame is:

.. math::
   :label: N_dot

   \left\{
   \begin{array}{l}
   {\dot{N}}_{k} = - W\ \nabla\cdot J_k^{\mathrm{lat}},\quad\quad\quad\quad k \in \left\lbrack 1,n \right\rbrack\\
   {\dot{N}}_{0} = - W\ \nabla\cdot J_{0}^{\mathrm{lat}} + \dot{N}.
   \end{array}
   \right.   

Using :numref:`Table %s <compvar>`, this is expressed in terms of site fractions:

.. math::
   :label: y_dot

   \left\{
   \begin{array}{l}
   \displaystyle{\dot{y}}_{k} = - \Omega\ \nabla\cdot J_k^{\mathrm{lat}} - y_k\frac{\dot{N}}{N},\quad\quad\quad\quad k \in \left\lbrack 1,n \right\rbrack\\
   \displaystyle{\dot{y}}_{0} = - \Omega\ \nabla\cdot J_{0}^{\mathrm{lat}} - y_0\frac{\dot{N}}{N} + \frac{\dot{N}}{N},
   \end{array}
   \right.

and then concentrations:

.. math::
   :label: c_dot

   \left\{
   \begin{array}{l}
   \displaystyle{\dot{c}}_{k} = - \nabla\cdot J_k^{\mathrm{lat}} - c_k\frac{\dot{W}}{W},\quad\quad\quad\quad k \in \left\lbrack 1,n \right\rbrack\\
   \displaystyle{\dot{c}}_{0} = - \nabla\cdot J_{0}^{\mathrm{lat}} - c_{0}\frac{\dot{W}}{W} + \frac{1}{\Omega}\frac{\dot{N}}{N}.
   \end{array}
   \right.
   
In deriving Eq. :eq:`c_dot`, we have also used the following relation:

.. math::
   :label: W_dot_W

   \frac{\dot{W}}{W} = \frac{\dot{N}}{N} + \frac{\dot{\Omega}}{\Omega},

which derives from Eq. :eq:`Omega` and reflects the fact that volume variations
are due to variations in the number of lattice sites and in the molar volume.
The rate of volume variation can be expressed as the divergence of the lattice
velocity:

.. math::
   :label: W_div_v

   \frac{\dot{W}}{W} = \nabla\cdot v.

The time derivatives in Eqs. :eq:`c_dot` are material (or total) derivatives:

.. math::
   :label: total_derivatives

   {\dot{c}}_{k} = \frac{\partial c_k}{\partial t} + v\ \nabla c_k.

Introducing Eq. :eq:`total_derivatives` in Eqs. :eq:`c_dot` and making use of Eq.
:eq:`W_div_v` yields the system of continuity equations in the laboratory
reference frame:

.. math::
   :label: continuity
   
   \left\{
   \begin{array}{l}
   \displaystyle\frac{\partial c_k}{\partial t} = - \nabla\cdot\left( J_k^{\mathrm{lat}} + c_kv \right),\quad\quad\quad\quad k \in \left\lbrack 1,n \right\rbrack\\
   \displaystyle\frac{\partial c_{0}}{\partial t} = - \nabla\cdot\left( J_{0}^{\mathrm{lat}} + c_{0}v \right) + \frac{1}{\Omega}\frac{\dot{N}}{N}
   \end{array}
   \right.

Equation :eq:`continuity` reflects the fact that atom concentrations vary
because of diffusion and advection, while the vacancy concentration may
additionally change due to the sink term :math:`\dot{N}/N`.

.. _lattice_velocity:

Lattice velocity
""""""""""""""""

The divergence of the velocity field is evaluated from the rate of volume
variation through Eqs. :eq:`W_div_v` and :eq:`W_dot_W`. The term
:math:`\dot{\Omega}/\Omega` in Eq. :eq:`W_dot_W` cannot be directly evaluated.
Instead, we note that combining Eqs. :eq:`W_sum` and :eq:`W_m`, the rate of
volume variations can be decomposed into contributions from the metal and the
pore, which are more easily accessed:

.. math::
   :label: div_v

   \nabla\cdot v = f_\mathrm{m}\left( \frac{\dot{N}}{N} + \frac{{\dot{\Omega}}_\mathrm{m}}{\Omega_\mathrm{m}} \right) + f_\mathrm{p}\frac{{\dot{W}}_\mathrm{p}}{W_\mathrm{p}}. 

Evaluating :math:`\nabla\cdot v` requires a model for the sink term,
:math:`\alpha = \dot{N}/N`. We use a linearized form and assume :math:`\alpha`
is proportional to the vacancy supersaturation,
:math:`y_0 - y_0^{\mathrm{eq}}`. We further postulate that
:math:`\alpha` can be expressed as the sum of two contributions, due to
lattice sinks and pores:

.. math::
   :label: sink_term

   \left\{
   \begin{array}{l}
   \dfrac{\dot{N}}{N} = \alpha = \alpha_\mathrm{d} + \alpha_\mathrm{p},\\
   \alpha_\mathrm{d} = - k_\mathrm{d}\left( y_0 - y_0^{\mathrm{eq}} \right),\\
   \alpha_\mathrm{p} = - k_\mathrm{p}\left( y_0 - y_0^{\mathrm{eq}} \right).
   \end{array}
   \right.

In Eq. :eq:`sink_term`, :math:`k_\mathrm{d}` and :math:`k_\mathrm{p}` (unit
:math:`s^{-1}`) are the sink strengths and reflect both the sink
densities and the frequency at which they operate; these parameters may
depend on alloy composition and microstructure. The index :math:`\mathrm{d}`
refers to dislocations: :math:`\alpha_\mathrm{d}`
produces the Kirkendall shift. The index :math:`\mathrm{p}` refers to pores.

The rate of pore volume variation is written:

.. math::
   :label: W_p_dot

   {\dot{W}}_\mathrm{p} = - N\alpha_\mathrm{p}\Omega_\mathrm{p}.

where :math:`\Omega_\mathrm{p}` is the partial molar volume of vacancies in the
pore. Like :math:`\Omega_0`, :math:`\Omega_\mathrm{p}` can be assigned either an
independent value or the local metal molar volume. The relative rate of pore
volume variation is:

.. math::
   :label: W_p_dot_w_p

   \frac{{\dot{W}}_\mathrm{p}}{W_\mathrm{p}} = - \frac{f_\mathrm{m}}{f_\mathrm{p}}\frac{\Omega_\mathrm{p}}{\Omega_\mathrm{m}}\alpha_\mathrm{p}.

The analytical expression of
:math:`\frac{{\dot{\Omega}}_\mathrm{m}}{\Omega_\mathrm{m}}` depends on the choice
made for :math:`\Omega_0`:

.. math::
   :label: Omega_m_dot_Omega_m

   \frac{{\dot{\Omega}}_\mathrm{m}}{\Omega_\mathrm{m}} = \left\{
   \begin{array}{ll}
    \displaystyle\frac{1}{f_\mathrm{m}(1 - y_0)}\sum_{k = 1}^{n}{\left(\Omega_\mathrm{m} - \Omega_k \right)\nabla\cdot J_k^{\mathrm{lat}}} & \mathrm{if}\ \Omega_0=\Omega_\mathrm{m},\\
    \displaystyle\frac{1}{f_\mathrm{m}}\sum_{k = 1}^{n}{(\Omega_0 - \Omega_k)\mathrm{div}\,J_k^\mathrm{lat}} + \alpha \frac{\Omega_0 - \Omega_\mathrm{m}}{\Omega_\mathrm{m}} & \mathrm{otherwise.}
   \end{array}
   \right.

Finally, the divergence of the velocity field is:

.. math::
   :label: div_v_detailed

   \nabla\cdot v = \left\{
   \begin{array}{ll}
    \displaystyle f_\mathrm{m}\left(\alpha_\mathrm{d} + \alpha_\mathrm{p}\frac{\Omega_\mathrm{m} - \Omega_\mathrm{p}}{\Omega_\mathrm{m}}\right) + \frac{1}{1 - y_0}\sum_{k = 1}^{n}{\left( \Omega_\mathrm{m} - \Omega_k \right)\nabla\cdot J_k^{\mathrm{lat}}} & \mathrm{if}\ \Omega_0=\Omega_\mathrm{m},\\
    \displaystyle f_\mathrm{m}\left(\alpha_\mathrm{d} \frac{\Omega_0}{\Omega_\mathrm{m}} + \alpha_\mathrm{p} \frac{\Omega_0 - \Omega_\mathrm{p}}{\Omega_\mathrm{m}}\right) + \sum_{k = 1}^{n}{(\Omega_0 - \Omega_k)\mathrm{div}\,J_k^\mathrm{lat}} & \mathrm{otherwise.}
   \end{array}
   \right.

Limiting case: ideal lattice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ideal lattice refers to a situation where vacancies are maintained at
equilibrium, everywhere and at all times, through the action of dislocations (no
pores are formed). This may be represented by :math:`k_\mathrm{d} = + \infty`.
If :math:`y_0^{\mathrm{eq}}` is composition-independent, and assuming that
initially :math:`y_0 = y_0^{\mathrm{eq}}`, the corresponding sink rate, noted
:math:`\alpha_{\infty}`, can be obtained analytically by solving 
:math:`{\dot{y}}_{0} = 0`. Using Eq. :eq:`y_dot`, this yields:

.. math::
   :label: alpha_ideal

   \alpha_{\infty} = \frac{\Omega}{1 - y_0^{\mathrm{eq}}}\nabla\cdot J_{0}^{\mathrm{lat}}.

If :math:`y_0^{\mathrm{eq}}` depends on the local composition, however,
:math:`\alpha_{\infty}` cannot be found by solving :math:`{\dot{y}}_{0} = 0`,
since :math:`y_0^{\mathrm{eq}}` changes over time due to composition variations.
In this case, :math:`\alpha_{\infty}` is obtained by iteration, as shown in
:ref:`implementation_ideal_lattice`.

Implementation
--------------

Thermodynamics
^^^^^^^^^^^^^^

Gibbs energies of pure metals and atom-atom interaction parameters are to be
obtained from assessed descriptions of the relevant systems published in
the literature. These are typically available as functions of atom fractions,
not site fractions, i.e., instead of :eq:`G_M`, the data are based on:

.. math::
   :label: G_M_x

   G_\mathrm{M}\left( \ldots x_{i}\ldots \right) = \sum_{i = 1}^{n}{x_{i}G_{i}} + RT\sum_{i = 1}^{n}{x_{i}\ln x_{i}} +_{}^{\mathrm{xs}}G_\mathrm{M}\left( \ldots x_{i}\ldots \right).

In principle, a conversion is required from the {atoms} description,
:math:`G_\mathrm{M}(\ldots x_{i}\ldots)`, to the {atoms + vacancy}
description of Eq. :eq:`G_M`, :math:`G_\mathrm{M}(\ldots y_i\ldots)`.
We assume :math:`y_0 \ll 1`, and
therefore :math:`y_k \cong x_k`. Consequently, we use Eq. :eq:`G_M_x` to
obtain interaction parameters, and assume these would be the same as
those obtained with Eq. :eq:`G_M`, i.e., we use these parameters in
Eq. :eq:`G_M` to compute equilibrium vacancy fractions and chemical potentials.

Similarly, calculating the equilibrium vacancy fraction using Eqs. :eq:`y0_eq`
requires the system composition to be described in terms of
site fractions. Most often, however, one needs to calculate
:math:`y_0^{\mathrm{eq}}` in a metal whose composition is known as mole
fractions of the metal species, not site fractions. Obtaining
:math:`y_k` from :math:`x_k` requires :math:`y_0`, which is the
unknown quantity. However, the assumption :math:`y_0 \ll 1` implies
:math:`y_k \cong x_k,\ \ k \in \left\lbrack 1,n \right\rbrack`. It
follows that :math:`y_0^{\mathrm{eq}}` can be directly computed from
atom fractions, without the need for an iterative scheme to solve Eq.
:eq:`y0_eq`.

Mobility
^^^^^^^^

The parameters describing the composition dependence of :math:`\phi_k`
in Eq. :eq:`lnDT_RK` are to be obtained from the literature. These are
typically available in terms of atom fractions:
:math:`\phi_k = \phi_k(\ldots x_{i}\ldots)`. They are directly used to
calculate tracer diffusion coefficients with Eqs. :eq:`lnDT` and :eq:`lnDT_RK`.
Non-equilibrium vacancy fractions are taken into account when calculating
:math:`L_{\mathrm{kk}}` via Eq. :eq:`Lkk`.

.. _implementation_diffusion:

Diffusion
^^^^^^^^^

The time evolution of the system composition is described by Eq.
:eq:`continuity`. This is solved using an explicit (forward Euler) finite
difference scheme, in one dimension of space. The total length is
divided into segments of size :math:`\Delta z_{i}` separated by grid
points at positions :math:`z_{i}`. The composition variables
(:math:`x_k,\ y_k,\ c_k,\ \Omega_{m},\ \Omega,\ f_{m},\mu_k`)
are associated with positions
:math:`z_i^\mathrm{m}=z_i + \Delta z_i/2`, and
represent the average system composition in the :math:`\Delta z_{i}`
segments. The fluxes (and the velocity field) are evaluated on grid
points :math:`z_{i}`; :math:`J_{i}^{\mathrm{lat}}` represents a flux
between segments :math:`\Delta z_{i - 1}` and
:math:`\Delta z_{i}`. The discretization is illustrated in
:numref:`Figure %s <grid>`.

.. _grid:

.. figure:: figures/grid.png
   :width: 300px
   :align: center

   Space discretization.

Fluxes are discretized using local averages of the volume size and Onsager
coefficients. This is done via the concept of the local diffusion resistance.
The continuous form:

.. math::
   :label: Jlat_1D
   
   J^{\mathrm{lat}} = - L\ \nabla \mu = - L\ \frac{\partial \mu}{\partial z}

is discretized into:

.. math::
   :label: Jlat_discrete
   
   J^{\mathrm{lat}}_i = - \frac{\mu_i - \mu_{i - 1}}{R_i},

where :math:`R_i` is the resistance between volumes :math:`i - 1` and :math:`i`.
(In Eq. :eq:`Jlat_discrete`, indices refer to positions on the grid, not chemical
species.) Noda provides several ways to compute :math:`R_i` from neighboring
:math:`\Delta z` and :math:`L` values, see :ref:`stencil`.

The initial conditions comprise profiles of the atom fractions :math:`x_k`,
vacancy site fraction :math:`y_0` and metal volume fraction :math:`f_{m}`.
These are used to compute :math:`y_k`, :math:`\Omega_{m}`, and then
:math:`\Omega` and :math:`c_k` through the composition relationships in
:ref:`composition_variables`.

The divergence of the velocity field is determined via Eq. :eq:`div_v_detailed`.
The velocity is obtained by integrating :math:`\nabla v` along the space
dimension, where the integration constant :math:`v(z_{\min})` depends on the
boundary condition. The divergence operator and its reciprocal depend on the
geometry, see definitions in :numref:`Table %s <divergence_geometries>` and
implementation in :func:`utils.div` and :func:`utils.integrate`.

.. _divergence_geometries:

.. table:: Divergence operator and its reciprocal in the three supported geometries. :math:`F` is a vector field, :math:`f` is its divergence.

   ============= ====================================================================== =====================================================
   Geometry                Divergence                                                             Reciprocal
   ============= ====================================================================== =====================================================
   Planar        :math:`\nabla F = \frac{\partial F}{\partial z}`                       :math:`F = \int f\mathrm{d}u`
   Cylindrical   :math:`\nabla F = \frac{1}{\rho}\frac{\partial \rho F}{\partial \rho}` :math:`F = \frac{1}{\rho} \int \rho f \mathrm{d}\rho`
   Spherical     :math:`\nabla F = \frac{1}{r^2}\frac{\partial r^2 F}{\partial r}`      :math:`F = \frac{1}{r^2} \int r^2 f \mathrm{d}r`
   ============= ====================================================================== =====================================================

The velocity field is then used to compute the time derivatives in Eq.
:eq:`continuity`. The :math:`n + 1` concentrations are related through the
closure relationship:

.. math::
   :label: closure_concentration
   
   \sum_{k = 0}^{n}c_k = \frac{1}{\Omega},

where :math:`\Omega` is the molar volume (see :ref:`composition_variables`).
However, :math:`\Omega` may vary over time, due to the variation of the metal molar
volume and volume fraction (see :ref:`composition_variables`). Obtaining one
of the concentrations through the :math:`n` others and :math:`\Omega` would
therefore require that the time evolution of :math:`\Omega` be determined
independently. Alternatively, we solve Eq. :eq:`continuity` for all
:math:`n + 1` concentrations, and then use Eq. :eq:`closure_concentration` to
obtain :math:`\Omega`. The :math:`c_k` and :math:`\Omega` values are then
used to compute :math:`y_k`, then :math:`\Omega_{m}` and :math:`f_{m}`.

.. _implementation_ideal_lattice:

Limiting case: ideal lattice
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ideal lattice configuration, represented
here by :math:`k_\mathrm{d} = + \infty`, is of particular interest for comparison
with the literature. This is readily implemented if
:math:`y_0^{\mathrm{eq}}` is constant, using Eq. :eq:`alpha_ideal` instead of Eq.
:eq:`sink_term` to compute :math:`\alpha`. In Noda,
however, :math:`y_0^{\mathrm{eq}}` is derived from the Gibbs free
energy of the metal, and is in general composition-dependent (see :ref:`thermo`).
An alternative method must therefore be implemented to compute
:math:`\alpha`. The solution retained here is to first assume
that :math:`y_0^{\mathrm{eq}}` is constant, compute
:math:`\alpha` via Eq. :eq:`alpha_ideal`, and solve the continuity equation
on this basis, to obtain a virtual :math:`c_k^{n + 1}`. Then, a
virtual :math:`y_0^{\mathrm{eq}}` is calculated from this
:math:`c_k^{n + 1}`, and a corrected :math:`\alpha` is calculated
such that
:math:`y_0^{n} + \Delta t \cdot {\dot{y}}_{0} = y_0^{\mathrm{eq}}`.
The continuity equation is then solved using the corrected
:math:`\alpha`. The new :math:`y_0^{\mathrm{eq}}` will necessarily
differ from the one calculated in the virtual step; however, the method
produces values of :math:`y_0` typically within
:math:`10^{- 13} - 10^{- 14}` of :math:`y_0^{\mathrm{eq}}`, which
is a satisfying approximation.

|
|
|
|

.. rubric:: References

.. [Agren_1982] J. Ågren, Journal of Physics and Chemistry of Solids 43 (1982) 421–430,
   `DOI: 10.1016/0022-3697(82)90152-4 <https://doi.org/10.1016/0022-3697(82)90152-4>`_

.. [Philibert_1991] J. Philibert, Atom movements - Diffusion and mass transport in solids,
   Les Editions de Physique, 1991.

.. [Andersson_1992] J.-O. Andersson, J. Ågren, Journal of Applied Physics 72 (1992) 1350–1355,
   `DOI: 10.1063/1.351745 <https://doi.org/10.1063/1.351745>`_

.. [Saunders_1998] N. Saunders, A.P. Miodownik, CALPHAD (Calculation of Phase Diagrams): A Comprehensive Guide,
   Pergamon Press, 1998.

.. [Svoboda_2006] J. Svoboda, F.D. Fischer, P. Fratzl, Acta Materialia 54 (11) (2006) 3043–3053,
   `DOI: 10.1016/j.actamat.2006.02.041 <https://doi.org/10.1016/j.actamat.2006.02.041>`_

.. [Lukas_2007] H. Lukas, S.G. Fries, B. Sundman, Computational Thermodynamics - The Calphad Method,
   Cambridge University Press, 2007.

.. [Franke_2014] P. Franke, J. Phase Equilib. Diffus. 35 (2014) 780–787,
   `DOI: 10.1007/s11669-014-0348-0 <https://doi.org/10.1007/s11669-014-0348-0>`_
   
.. [Abe_2018] T. Abe, K. Hashimoto, M. Shimono, Materials Transactions 59 (2018) 580–584,
   `DOI: 10.2320/matertrans.M2017328 <https://doi.org/10.2320/matertrans.M2017328>`_

.. [Gheno_2022] T. Gheno, V. Szczepan, C. Salsi, C. Desgranges, D. Monceau,
   Computational Materials Science 215 (2022) 111785,
   `DOI: 10.1016/j.commatsci.2022.111785 <https://doi.org/10.1016/j.commatsci.2022.111785>`_
