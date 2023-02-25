.. _methods_depletion:

=========
Depletion
=========

When materials in a system are subject to irradiation over a long period of
time, nuclides within the material will transmute due to nuclear reactions as
well as spontaneous radioactive decay. The time-dependent process by which
nuclides transmute under irradiation is known as *depletion* or *burnup*. To
accurately analyze nuclear systems, it is often necessary to predict how the
composition of materials will change since this change results in a
corresponding change in the solution of the transport equation. The equation
that governs the transmutation and decay of nuclides inside of an irradiated
environment can be written as

.. math::

    \begin{aligned} \frac{dN_i(t)}{dt} = &\sum\limits_j
    \underbrace{\left [ \underbrace{f_{j \rightarrow i} \int_0^\infty dE \;
    \sigma_j (E, t) \phi(E,t)}_\text{transmutation} +
    \underbrace{\lambda_{j\rightarrow i}}_\text{decay} \right ]
    N_j(t)}_{\text{Production of nuclide }i\text{ from nuclide }j} \\
    &- \underbrace{\left [\underbrace{\int_0^\infty dE \; \sigma_i
    (E,t) \phi(E,t)}_\text{transmutation} +
    \underbrace{\sum\limits_j \lambda_{i\rightarrow j}}_\text{decay} \right ]
    N_i(t)}_{\text{Loss of nuclide }i} \end{aligned}

where :math:`N_i` is the density of nuclide :math:`i` at time :math:`t`,
:math:`\sigma_i` is the transmutation cross section for nuclide :math:`i` at
energy :math:`E`, :math:`f_{j \rightarrow i}` is the fraction of transmutation
reactions in nuclide :math:`j` that produce nuclide :math:`i`, and
:math:`\lambda_{j \rightarrow i}` is the decay constant for decay modes in
nuclide :math:`j` that produce nuclide :math:`i`. Note that we have not included
the spatial dependence of the flux or cross sections. As one can see, the
equation simply states that the rate of change of :math:`N_i` is equal to the
production rate minus the loss rate. Because the equation for nuclide :math:`i`
depends on the nuclide density for possibly many other nuclides, we have a
system of first-order differential equations. To form a proper initial value
problem, we also need the nuclide densities at time 0:

.. math::

    N_i(0) = N_{i,0}.

These equations can be written more compactly in matrix notation as

.. math::
    :label: depletion-matrix

    \frac{d\mathbf{n}}{dt} = \mathbf{A}(\mathbf{n},t)\mathbf{n}, \quad \mathbf{n}(0) =
    \mathbf{n}_0

where :math:`\mathbf{n} \in \mathbb{R}^n` is the nuclide density vector,
:math:`\mathbf{A}(\mathbf{n},t) \in \mathbb{R}^{n\times n}` is the burnup matrix
containing the decay and transmutation coefficients, and :math:`\mathbf{n}_0` is
the initial density vector. Note that the burnup matrix depends on
:math:`\mathbf{n}` because the solution to the transport equation depends on the
nuclide densities.

.. _methods_depletion_integration:

---------------------
Numerical Integration
---------------------

A variety of numerical methods exist for solving Eq. :eq:`depletion-matrix`. The
simplest such method, known as the "predictor" method, is to divide the overall
time interval of interest :math:`[0,t]` into smaller timesteps over which it is
assumed that the burnup matrix is constant. Let :math:`t \in [t_i, t_i + h]` be
one such timestep. Over the timestep, the solution to Eq. :eq:`depletion-matrix`
can be written analytically using the matrix exponential

.. math::

    \mathbf{A}_i = \mathbf{A}(\mathbf{n}_i, t_i) \\

    \mathbf{n}_{i+1} = e^{\mathbf{A}_i h} \mathbf{n}_i

where :math:`\mathbf{n}_i \equiv \mathbf{n}(t_i)`. The exponential of a matrix
:math:`\mathbf{X}` is defined by the power series expansion

.. math::

    e^{\mathbf{X}} = \sum\limits_{k=0}^\infty \frac{1}{k!} \left ( \mathbf{X}
    \right )^k

where :math:`\mathbf{X}^0 = \mathbf{I}`. A series of so-called
predictor-corrector methods that use multiple stages offer improved accuracy
over the predictor method. The simplest of these methods, the CE/CM algorithm,
is defined as

.. math::

    \mathbf{n}_{i+1/2} = e^{\frac{h}{2}\mathbf{A}(\mathbf{n}_i, t_i)} \mathbf{n}_i \\
    \mathbf{n}_{i+1} = e^{h \mathbf{A}(\mathbf{n}_{i+1/2},t_{i+1/2})} \mathbf{n}_i

Here, the value of :math:`\mathbf{n}` at the midpoint is estimated using
:math:`\mathbf{A}` evaluated at the beginning of the timestep. Then,
:math:`\mathbf{A}` is evaluated using the densities at the midpoint and used to
integrate over the entire timestep.

Our aim here is not to exhaustively describe all integration methods but rather
to give a few examples that elucidate the main considerations one must take into
account when choosing a method. Generally, there is a tradeoff between the
accuracy of the method and its computational expense. In the case of
transport-coupled depletion, the expense is driven almost entirely by the time
to compute a transport solution, i.e., to evaluate :math:`\mathbf{A}` for a
given :math:`\mathbf{n}`. Thus, the cost of a method scales with the number of
:math:`\mathbf{A}` evaluations that are performed per timestep. On the other
hand, methods that require more evaluations generally achieve higher accuracy.
The predictor method only requires one evaluation and its error converges as
:math:`\mathcal{O}(h)`. The CE/CM method requires two evaluations and is thus
twice as expensive as the predictor method, but achieves an error of
:math:`\mathcal{O}(h^2)`. An exhaustive description of time integration methods
and their merits can be found in the `thesis of Colin Josey
<http://dspace.mit.edu/handle/1721.1/7582>`_.

OpenMC does not rely on a single time integration method but rather has several
classes that implement different algorithms. For example, the
:class:`openmc.deplete.PredictorIntegrator` class implements the predictor
method, and the :class:`openmc.deplete.CECMIntegrator` class implements the
CE/CM method. A full list of the integrator classes available can be found in
the documentation for the :mod:`openmc.deplete` module.

------------------
Matrix Exponential
------------------

As we saw in the :ref:`previous section <methods_depletion_integration>`,
numerically integrating Eq. :eq:`depletion-matrix` requires evaluating one or
more matrix exponentials. OpenMC uses the Chebyshev rational approximation
method (CRAM), which was introduced in a series of papers by Pusa (`1
<https://doi.org/10.13182/NSE09-14>`_, `2
<https://doi.org/10.13182/NSE10-81>`_), to evaluate matrix exponentials. In
particular, OpenMC utilizes an `incomplete partial fraction
<https://doi.org/10.13182/NSE15-26>`_ (IPF) form of CRAM that provides a good
balance of numerical stability and efficiency. In this representation the matrix
exponential is approximated as

.. math::

    e^{\mathbf{A}t} \approx \alpha_0 \prod\limits_{\ell=1}^{k/2} \left (
    \mathbf{I} + 2 \text{Re} \left ( \widetilde{\alpha}_\ell \left (\mathbf{A}t
    - \theta_\ell \mathbf{I} \right )^{-1} \right ) \right )

where :math:`k` is the order of the approximation and :math:`\alpha_0`,
:math:`\widetilde{\alpha}_\ell`, and :math:`\theta_\ell` are coefficients that
have been tabulated for orders up to :math:`k=48`. Rather than computing the
full approximation and then multiplying it by a vector, the following algorithm
is used to incrementally apply the terms within the product (note that the
original description of the algorithm presented by `Pusa
<https://doi.org/10.13182/NSE15-26>`_ contains a typo):

1. :math:`\mathbf{n} \gets \mathbf{n_0}`
2. For :math:`\ell = 1, 2, \dots, k/2`

   - :math:`\mathbf{n} \gets \mathbf{n} + 2\text{Re}(\widetilde{\alpha}_\ell
     (\mathbf{A}t - \theta_\ell)^{-1})\mathbf{n}`

3. :math:`\mathbf{n} \gets \alpha_0 \mathbf{n}`

The :math:`k`\ th order approximation for CRAM requires solving :math:`k/2`
sparse linear systems. OpenMC relies on functionality from
:mod:`scipy.sparse.linalg` for solving the linear systems.

-------------------
Data Considerations
-------------------

In principle, solving Eq. :eq:`depletion-matrix` using CRAM is fairly simple:
just construct the burnup matrix at various times and solve a set of sparse
linear systems. However, constructing the burnup matrix itself involves not
only solving the transport equation to estimate transmutation reaction rates
(in the case of transport-coupled depletion) or to obtain microscopic cross
sections (in the case of transport-independent depletion), but also a series of
choices about what data to include. In OpenMC, the burnup matrix is constructed
based on data inside of a *depletion chain* file, which includes fundamental
data gathered from ENDF incident neutron, decay, and fission product yield
sublibraries. For each nuclide, this file includes:

- What transmutation reactions are possible, their Q values, and their products;
- If a nuclide is not stable, what decay modes are possible, their branching
  ratios, and their products; and
- If a nuclide is fissionable, the fission products yields at any number of
  incident neutron energies.

Transmutation Reactions
-----------------------

In transport-coupled depletion, OpenMC will setup tallies in a problem based on
what transmutation reactions are available in a depletion chain file, so any
arbitrary number of transmutation reactions can be tracked. In
transport-independent depletion, OpenMC will calculate reaction rates for every
reaction that is present in both the available cross sections and the depletion
chain file. The pregenerated chain files that are available on
https://openmc.org include the following transmutation reactions: fission, (n,\
:math:`\gamma`\ ), (n,2n), (n,3n), (n,4n), (n,p), and (n,\ :math:`\alpha`\ ).

Capture Branching Ratios
------------------------

Some (n,\ :math:`\gamma`\ ) reactions may result in a product being in either the
ground or a metastable state. The most well-known example is capture in Am241,
which can produce either Am242 or Am242m. Because the metastable state of Am242m
has a significantly longer half-life than the ground state, it is important to
accurately model the branching of the capture reaction in Am241. This is
complicated by the fact that the branching ratio may depend on the incident
neutron energy causing capture.

OpenMC's transport solver does not currently allow energy-dependent capture
branching ratios. However, the depletion chain file does allow a transmutation
reaction to be listed multiple times with different branching ratios resulting
in different products. Spectrum-averaged capture branching ratios have been
computed in LWR and SFR spectra and are available at
https://openmc.org/depletion-chains.

Fission Product Yields
----------------------

Fission product yields (FPY) are also energy-dependent in general. ENDF fission
product yield sublibraries typically include yields tabulated at 2 or 3
energies. It is an open question as to what the best way to handle this energy
dependence is. OpenMC includes three methods for treating the energy dependence
of FPY:

1. Use FPY data corresponding to a specified energy. This is used by default in
   both transport-coupled and transport-independent depletion.
2. Tally fission rates above and below a specified cutoff energy. Assume that
   all fissions below the cutoff energy correspond to thermal FPY data and all
   fission above the cutoff energy correspond to fast FPY data. Only applicable
   to transport-coupled depletion.
3. Compute the average energy at which fission events occur and use an effective
   FPY by linearly interpolating between FPY provided at neighboring energies.
   Only applicable to transport-coupled depletion.

The method for transport-coupled depletion can be selected through the
``fission_yield_mode`` argument to the :class:`openmc.deplete.CoupledOperator`
constructor.

Power Normalization
-------------------

In transport-coupled depletion, the reaction rates provided OpenMC are given in
units of reactions per source particle. For depletion, it is necessary to
compute an absolute reaction rate in reactions per second. To do so, the
reaction rates are normalized based on a specified power. A complete
description of how this normalization can be performed is described in
:ref:`usersguide_tally_normalization`. Here, we simply note that the main
depletion class, :class:`openmc.deplete.CoupledOperator`, allows the user to
choose one of two methods for estimating the heating rate, including:

1. Using fixed Q values from a depletion chain file (useful for comparisons to
   other codes that use fixed Q values), or
2. Using the ``heating`` or ``heating-local`` scores to obtain an nuclide- and
   energy-dependent estimate of the true heating rate.

The method for normalization can be chosen through the ``normalization_mode``
argument to the :class:`openmc.deplete.CoupledOperator` class.

----------------
MSR Capabilities
----------------

A nice feature of Molten Salt Reactors (MSRs), and in general of
*any*  liquid fuel reactor, is the capability to continuously reprocess the fuel
to remove fission products or to feed fresh fuel into the system.

Removal rate
------------

A removal rate as defined here is the rate at which nuclides are
continuously removed from a material. Mathematically, it
can be thought as an additional proportional term :math:`\mathbf{T}`
to the depletion equation, which can be written as:

.. math::

  \begin{aligned}\frac{dN_i(t)}{dt} = &\underbrace{\sum\limits_j f_{j\rightarrow i}
  \int_0^\infty dE  \; \sigma_j (E,t) \phi(E,t) N_j(t)  - \int_0^\infty dE \; \sigma_i(E,t)
  \phi(E,t) N_i(t)}_\textbf{R} \\
  &+ \underbrace{\sum_j \left [ \lambda_{j\rightarrow i} N_j(t) - \lambda_{i\rightarrow j} N_i(t) \right ]}_\textbf{D} \\
  &- \underbrace{\epsilon_i \lambda_i N_i(t)}_\textbf{T} \end{aligned}

where the transmutation term :math:`\mathbf{R}`, the decay term :math:`\mathbf{D}`
and the new removal term :math:`\mathbf{T}` have been grouped together so that
:math:`\mathbf{A} = \mathbf{R}+\mathbf{D}-\mathbf{T}`.
The removal efficiency :math:`\epsilon_i` and the removal rate coefficient
:math:`\lambda_i` define the continuous removal of the nuclide :math:`i`, which
behaviour is similar to radioactive decay.
:math:`\lambda_i` can also be defined as the reciprocal of a cycle time
:math:`T_{cyc}`, intended as the time needed to process the whole inventory.

For simplicity, :math:`\epsilon_i` and :math:`\lambda_i`
can be combined together in one single user-defined parameter that is defined
again as :math:`\lambda_i`.
Thus, setting a removal rate coefficient of :math:`1 s^{-1}` at :math:`100\%`
efficiency, would be the same as setting :math:`10 s^{-1}` at :math:`10\%`.

Note that this formulation assumes first order removal and
homogeneous distribution of nuclide :math:`i` throughout the material.

A more rigorous description of removal rate and its implementation can be found
in the paper by `Hombourger
<https://www.sciencedirect.com/science/article/pii/S0306454920302024?via%3Dihub>`_.

The resulting burnup matrix can be solved with the same integration algorithms,
as we've seen before.

.. note::

    If no ``destination_material`` argument is passed to the
    :meth:`~openmc.deplete.MsrContinuous.set_removal_rate()` method of the
    :class:`~openmc.deplete.MsrContinuous` class, nuclides that are removed will
    not be tracked afterwards.

Feeding rate
------------
To keep track of removed nuclides or to feed nuclides from one depletable material
to another, the respective depletion equations have to be coupled. This can be
achieved by defining one multidimensional square matrix with dimensions equal to
the number of depletable materials that transfer nuclides.
The diagonal positions are filled with the usual depletion matrices
:math:`\mathbf{A_{ii}}`, where the index :math:`i` indicates the depletable
material id, and the off-diagonal positions are filled with the removal matrices
:math:`\mathbf{T_{ij}}`, positioned so that that the indices :math:`i` and
:math:`j` indicate the nuclides receiving and loosing materials, respectively.
The nuclide vectors are assembled together in one single vector and the resulting
system is solved with the same integration algorithms seen before.

For sake of example, let's consider the  case of two depletable materials and one
feed defined from material 1 to material 2. The final system will look like:

.. math::

  \begin{aligned}\frac{d}{dt}\begin{pmatrix}\vec{N_1}\\ \vec{N_2}\end{pmatrix} &=
  \begin{pmatrix}\mathbf{A_{11}} & \mathbf{0}\\ \mathbf{T_{21}} & \mathbf{A_{22 }}
  \end{pmatrix} \begin{pmatrix}\vec{N_1}\\ \vec{N_2}\end{pmatrix} \end{aligned}

where:

:math:`\mathbf{A_{11}} = \mathbf{R_{11}}+\mathbf{D_{11}}-\mathbf{T_{21}}`, and

:math:`\mathbf{A_{22}} = \mathbf{R_{22}}+\mathbf{D_{22}}`.

Note that mass conservation during feeding is assured by transferring the number
of atoms directly.
