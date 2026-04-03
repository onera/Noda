Working with model systems
==========================

.. contents:: :local:

Users can easily explore the effects of thermodynamic and diffusion properties
on the shape of concentration profiles by creating their own database files.
The examples below shows how the known solution to the diffusion problem
is recovered when the system properties follow a particular set of constraints.

Let us consider a single-phase binary system AB, where the partial molar
volumes of A abd B are equal (therefore the average molar volume is
composition-independent). In 1D, the diffusion problem may be written in the
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

These conditions are fulfilled for the model AB system provided in the test
base, with databases named "AB_thermo_ideal" and "AB_mob_ideal".

A number of diffusion problems involving a binary system with
composition-independent molar volume and diffusivity have analytical solutions.
Some typical initial distributions and boundary conditions are illustrated
below.

Diffusion couple (planar geometry)
----------------------------------

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

.. literalinclude:: /../../tests/jobs/couple_AB/couple_AB.toml
   :caption:

.. literalinclude:: /../../tests/jobs/couple_AB/run.py
   :caption:

.. figure:: /../../tests/jobs/couple_AB/couple_AB.png
    :width: 400px
    :align: center

|

.. _source_AB:

Constant surface concentration (planar geometry)
------------------------------------------------

In the case of a semi-infinite solid initially at a uniform concentration
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

.. literalinclude:: /../../tests/jobs/source_AB/source_AB.toml
   :caption:

.. literalinclude:: /../../tests/jobs/source_AB/run.py
   :caption:

.. figure:: /../../tests/jobs/source_AB/source_AB.png
    :width: 400px
    :align: center

|

Constant surface concentration (spherical geometry)
---------------------------------------------------

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

.. literalinclude:: /../../tests/jobs/sphere_AB/sphere_AB.toml
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

