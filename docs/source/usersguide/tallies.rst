.. _usersguide_tallies:

==================
Specifying Tallies
==================

.. currentmodule:: openmc

In order to obtain estimates of physical quantities in your simulation, you need
to create one or more tallies using the :class:`openmc.Tally` class. As
explained in detail in the :ref:`theory manual <methods_tallies>`, tallies
provide estimates of a scoring function times the flux integrated over some
region of phase space, as in:

.. math::

    X = \underbrace{\int d\mathbf{r} \int d\mathbf{\Omega} \int
    dE}_{\text{filters}} \underbrace{f(\mathbf{r}, \mathbf{\Omega},
    E)}_{\text{scores}} \psi (\mathbf{r}, \mathbf{\Omega}, E)

Thus, to specify a tally, we need to specify what regions of phase space should
be included when deciding whether to score an event as well as what the scoring
function (:math:`f` in the above equation) should be used. The regions of phase
space are generally called *filters* and the scoring functions are simply
called *scores*.

The only cases when filters do not correspond directly with the regions of
phase space are when expansion functions are applied in the integrand, such as
for Legendre expansions of the scattering kernel.

-------
Filters
-------

To specify the regions of phase space, one must create a
:class:`openmc.Filter`. Since :class:`openmc.Filter` is an abstract class, you
actually need to instantiate one of its sub-classes (for a full listing, see
:ref:`pythonapi_tallies`). For example, to indicate that events that occur in a
given cell should score to the tally, we would create a
:class:`openmc.CellFilter`::

   cell_filter = openmc.CellFilter([fuel.id, moderator.id, reflector.id])

Another commonly used filter is :class:`openmc.EnergyFilter`, which specifies
multiple energy bins over which events should be scored. Thus, if we wanted to
tally events where the incident particle has an energy in the ranges [0 eV, 4
eV] and [4 eV, 1 MeV], we would do the following::

  energy_filter = openmc.EnergyFilter([0.0, 4.0, 1.0e6])

Energies are specified in eV and need to be monotonically increasing.

.. caution:: An energy bin between zero and the lowest energy specified is not
             included by default as it is in MCNP.

Once you have created a filter, it should be assigned to a :class:`openmc.Tally`
instance through the :attr:`Tally.filters` attribute::

  tally.filters.append(cell_filter)
  tally.filters.append(energy_filter)

  # This is equivalent
  tally.filters = [cell_filter, energy_filter]

.. note:: You are actually not required to assign any filters to a tally. If you
          create a tally with no filters, all events will score to the
          tally. This can be useful if you want to know, for example, a reaction
          rate over your entire model.

.. _usersguide_lifetime_moments:

Time Moments Since Birth
------------------------

The :class:`openmc.LifetimeMomentFilter` multiplies the scores of a tally by
the powers :math:`\tau^n`, :math:`n = 0, 1, \ldots, N`, of the time
:math:`\tau` that has elapsed since the scoring particle was born, where
:math:`N` is the order of the filter. Bin :math:`n` of the filter therefore
estimates the :math:`n`-th moment of the score with respect to the time since
birth, in units of the score multiplied by :math:`\text{s}^n`. Bin 0 is the
score itself. All bins are obtained with the same estimator from the same
events, so their statistical errors are correlated, which generally makes
ratios of bins more precise than ratios of independent tallies.

The time since birth, :math:`\tau`, is measured in seconds from the moment the
particle was started from its source site. For a source particle in a
fixed-source calculation this is the time elapsed since its emission, and the
fission neutrons transported in a fixed-source calculation are likewise
started from their fission sites. In an eigenvalue calculation, a neutron is
started from the fission site that produced it, so for a delayed neutron
:math:`\tau` counts from its emission rather than from the fission that created
its precursor.

.. note:: The time since birth is not carried over to secondary or split
   particles that continue a history. The additional neutrons of
   :math:`(n,xn)` reactions and the particles created by weight-window
   splitting restart from :math:`\tau = 0` at the collision or split that
   created them, instead of continuing from the birth of the neutron of the
   source site. This is a limitation of the implementation rather than part of
   the definition of the moments: when such reactions or weight windows are
   present, the moments of order :math:`n \ge 1` are biased low. The birth
   position and birth cell used by :class:`openmc.MeshBornFilter` and
   :class:`openmc.CellBornFilter`, and the neutron lifetime of the iterated
   fission probability method, have the same limitation.

Because the filter weight changes along a track, a tally with this filter
cannot use a track-length estimator. A tally read from ``tallies.xml`` gets the
estimator that it would get without this filter, except that a collision
estimator replaces the default track-length one. For example, a tally with
only this filter and a ``flux`` score uses a collision estimator, and one with
this filter, an :class:`openmc.EnergyoutFilter`, and a ``nu-fission`` score
uses an analog estimator. Requesting a track-length estimator for such a tally
in ``tallies.xml`` results in an error. A tally created at run time through
:mod:`openmc.lib` that has a track-length estimator, whether by default or set
explicitly, is instead switched, with a warning, when the simulation is
initialized, to the estimator that the same tally gets from ``tallies.xml``.
Since collision and analog estimators only score at collisions, a tally with
this filter scores nothing in void regions. In surface current tallies on a mesh
(:class:`openmc.MeshSurfaceFilter`), :math:`\tau` is the time since birth at
the end of the flight that crossed the mesh surface rather than at the
crossing, as for :class:`openmc.TimeFilter`.

A typical use is to obtain fission-to-fission time moments in an eigenvalue
calculation by combining the filter with a fission production score. The
following tallies give, for each pair of mesh elements, the number of fission
neutrons produced in one element by neutrons born in another one, together
with the first and second moments of the time from the birth of the incident
neutron to the fission event::

  born = openmc.MeshBornFilter(mesh)
  where = openmc.MeshFilter(mesh)
  moments = openmc.LifetimeMomentFilter(2)

  tally = openmc.Tally(name='fission time moments')
  tally.filters = [born, where, moments]
  tally.scores = ['nu-fission', 'prompt-nu-fission']

  delayed = openmc.Tally(name='delayed fission time moments')
  delayed.filters = [born, where,
                     openmc.DelayedGroupFilter([1, 2, 3, 4, 5, 6]), moments]
  delayed.scores = ['delayed-nu-fission']

For each element pair, the ratio of bin 1 to bin 0 is the mean time from the
birth of a neutron to a fission that it causes, weighted by the number of
fission neutrons produced, and bin 2 provides the second moment needed, for
example, for the variance of that time. The time ends at the fission event:
for the delayed neutrons produced, the decay time of their precursors is not
included.

.. _usersguide_scores:

------
Scores
------

To specify the scoring functions, a list of strings needs to be given to the
:attr:`Tally.scores` attribute. You can score the flux ('flux'), or a reaction
rate ('total', 'fission', etc.). For example, to tally the elastic scattering
rate and the fission neutron production, you'd assign::

  tally.scores = ['elastic', 'nu-fission']

With no further specification, you will get the total elastic scattering rate
and the total fission neutron production. If you want reaction rates for a
particular nuclide or set of nuclides, you can set the :attr:`Tally.nuclides`
attribute to a list of strings indicating which nuclides. The nuclide names
should follow the same :ref:`naming convention <usersguide_naming>` as that used
for material specification. If we wanted the reaction rates only for U235 and
U238 (separately), we'd set::

  tally.nuclides = ['U235', 'U238']

To get a separate reaction rate for every nuclide in a material, ask the
material for its nuclides::

  tally.nuclides = material.get_nuclides()

Note that these are nuclide names; an element symbol such as 'Al' is not a valid
entry, even for photon tallies, whose data is tabulated per element.

The following tables show all valid scores:

.. table:: **Flux scores: units are particle-cm per source particle.**

    +----------------------+---------------------------------------------------+
    |Score                 | Description                                       |
    +======================+===================================================+
    |flux                  |Total flux.                                        |
    +----------------------+---------------------------------------------------+

.. table:: **Reaction scores: units are reactions per source particle.**

    +------------------------+-------------------------------------------------+
    |Score                   |Description                                      |
    +========================+=================================================+
    |absorption              |Total absorption rate. For incident neutrons,    |
    |                        |this accounts for all reactions that do not      |
    |                        |produce secondary neutrons as well as fission.   |
    |                        |For incident photons, this includes              |
    |                        |photoelectric and pair production.               |
    +------------------------+-------------------------------------------------+
    |elastic                 |Elastic scattering reaction rate.                |
    +------------------------+-------------------------------------------------+
    |fission                 |Total fission reaction rate.                     |
    +------------------------+-------------------------------------------------+
    |scatter                 |Total scattering rate.                           |
    +------------------------+-------------------------------------------------+
    |total                   |Total reaction rate.                             |
    +------------------------+-------------------------------------------------+
    |(n,2nd)                 |(n,2nd) reaction rate.                           |
    +------------------------+-------------------------------------------------+
    |(n,2n)                  |(n,2n) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,3n)                  |(n,3n) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,na)                  |(n,n\ :math:`\alpha`\ ) reaction rate.           |
    +------------------------+-------------------------------------------------+
    |(n,n3a)                 |(n,n3\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,2na)                 |(n,2n\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,3na)                 |(n,3n\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,np)                  |(n,np) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,n2a)                 |(n,n2\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,2n2a)                |(n,2n2\ :math:`\alpha`\ ) reaction rate.         |
    +------------------------+-------------------------------------------------+
    |(n,nd)                  |(n,nd) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,nt)                  |(n,nt) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,n3He)                |(n,n\ :sup:`3`\ He) reaction rate.               |
    +------------------------+-------------------------------------------------+
    |(n,nd2a)                |(n,nd2\ :math:`\alpha`\ ) reaction rate.         |
    +------------------------+-------------------------------------------------+
    |(n,nt2a)                |(n,nt2\ :math:`\alpha`\ ) reaction rate.         |
    +------------------------+-------------------------------------------------+
    |(n,4n)                  |(n,4n) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,2np)                 |(n,2np) reaction rate.                           |
    +------------------------+-------------------------------------------------+
    |(n,3np)                 |(n,3np) reaction rate.                           |
    +------------------------+-------------------------------------------------+
    |(n,n2p)                 |(n,n2p) reaction rate.                           |
    +------------------------+-------------------------------------------------+
    |(n,npa)                 |(n,np\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,n*X*)                |Level inelastic scattering reaction rate. The    |
    |                        |*X* indicates which inelastic level, e.g.,       |
    |                        |(n,n3) is third-level inelastic scattering.      |
    +------------------------+-------------------------------------------------+
    |(n,nc)                  |Continuum level inelastic scattering             |
    |                        |reaction rate.                                   |
    +------------------------+-------------------------------------------------+
    |(n,gamma)               |Radiative capture reaction rate.                 |
    +------------------------+-------------------------------------------------+
    |(n,p)                   |(n,p) reaction rate.                             |
    +------------------------+-------------------------------------------------+
    |(n,d)                   |(n,d) reaction rate.                             |
    +------------------------+-------------------------------------------------+
    |(n,t)                   |(n,t) reaction rate.                             |
    +------------------------+-------------------------------------------------+
    |(n,3He)                 |(n,\ :sup:`3`\ He) reaction rate.                |
    +------------------------+-------------------------------------------------+
    |(n,a)                   |(n,\ :math:`\alpha`\ ) reaction rate.            |
    +------------------------+-------------------------------------------------+
    |(n,2a)                  |(n,2\ :math:`\alpha`\ ) reaction rate.           |
    +------------------------+-------------------------------------------------+
    |(n,3a)                  |(n,3\ :math:`\alpha`\ ) reaction rate.           |
    +------------------------+-------------------------------------------------+
    |(n,2p)                  |(n,2p) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,pa)                  |(n,p\ :math:`\alpha`\ ) reaction rate.           |
    +------------------------+-------------------------------------------------+
    |(n,t2a)                 |(n,t2\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,d2a)                 |(n,d2\ :math:`\alpha`\ ) reaction rate.          |
    +------------------------+-------------------------------------------------+
    |(n,pd)                  |(n,pd) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,pt)                  |(n,pt) reaction rate.                            |
    +------------------------+-------------------------------------------------+
    |(n,da)                  |(n,d\ :math:`\alpha`\ ) reaction rate.           |
    +------------------------+-------------------------------------------------+
    |photon-total            |Total photo-atomic reaction rate.                |
    +------------------------+-------------------------------------------------+
    |coherent-scatter        |Coherent (Rayleigh) scattering reaction rate.    |
    +------------------------+-------------------------------------------------+
    |incoherent-scatter      |Incoherent (Compton) scattering reaction rate.   |
    +------------------------+-------------------------------------------------+
    |photoelectric           |Photoelectric absorption reaction rate.          |
    +------------------------+-------------------------------------------------+
    |photoelectric-*S*       |Subshell photoelectric absorption rate for the   |
    |                        |*S* shell. For example, "photoelectric-N3" is the|
    |                        |rate for the N3 subshell.                        |
    +------------------------+-------------------------------------------------+
    |pair-production         |Pair production reaction rate (total).           |
    +------------------------+-------------------------------------------------+
    |pair-production-electron|Pair production reaction rate in the electron    |
    |                        |field.                                           |
    +------------------------+-------------------------------------------------+
    |pair-production-nuclear |Pair production reaction rate in the nuclear     |
    |                        |field.                                           |
    +------------------------+-------------------------------------------------+
    |*Arbitrary integer*     |An arbitrary integer is interpreted to mean the  |
    |                        |reaction rate for a reaction with a given ENDF   |
    |                        |MT number.                                       |
    +------------------------+-------------------------------------------------+

.. table:: **Particle production scores: units are particles produced per
           source particles.**

    +----------------------+---------------------------------------------------+
    |Score                 | Description                                       |
    +======================+===================================================+
    |delayed-nu-fission    |Total production of delayed neutrons due to        |
    |                      |fission.                                           |
    +----------------------+---------------------------------------------------+
    |prompt-nu-fission     |Total production of prompt neutrons due to         |
    |                      |fission.                                           |
    +----------------------+---------------------------------------------------+
    |nu-fission            |Total production of neutrons due to fission.       |
    +----------------------+---------------------------------------------------+
    |nu-scatter            |This score is similar in functionality to the      |
    |                      |``scatter`` score except the total production of   |
    |                      |neutrons due to scattering is scored vice simply   |
    |                      |the scattering rate. This accounts for             |
    |                      |multiplicity from (n,2n), (n,3n), and (n,4n)       |
    |                      |reactions.                                         |
    +----------------------+---------------------------------------------------+
    |H1-production         |Total production of H1.                            |
    +----------------------+---------------------------------------------------+
    |H2-production         |Total production of H2 (deuterium).                |
    +----------------------+---------------------------------------------------+
    |H3-production         |Total production of H3 (tritium).                  |
    +----------------------+---------------------------------------------------+
    |He3-production        |Total production of He3.                           |
    +----------------------+---------------------------------------------------+
    |He4-production        |Total production of He4 (alpha particles).         |
    +----------------------+---------------------------------------------------+

.. table:: **Miscellaneous scores: units are indicated for each.**

    +----------------------+---------------------------------------------------+
    |Score                 | Description                                       |
    +======================+===================================================+
    |current               |It may not be used in conjunction with any other   |
    |                      |score except flux.                                 |
    |                      |                                                   |
    |                      |When used in combination with a meshsurface filter:|
    |                      |Partial currents on the boundaries of each cell in |
    |                      |a mesh.                                            |
    |                      |                                                   |
    |                      |When used in combination with a surface filter:    |
    |                      |Net currents on any surface previously defined in  |
    |                      |the geometry. It may be used along with any other  |
    |                      |filter, except meshsurface filters.                |
    |                      |Surfaces can alternatively be defined with cell    |
    |                      |from and cell filters thereby resulting in tallying|
    |                      |partial currents.                                  |
    |                      |                                                   |
    |                      |Units are particles per source particle.           |
    +----------------------+---------------------------------------------------+
    |events                |Number of scoring events. Units are events per     |
    |                      |source particle.                                   |
    +----------------------+---------------------------------------------------+
    |inverse-velocity      |The flux-weighted inverse velocity where the       |
    |                      |velocity is in units of centimeters per second.    |
    +----------------------+---------------------------------------------------+
    |heating               |Total nuclear heating in units of eV per source    |
    |                      |particle. For neutrons, this corresponds to MT=301 |
    |                      |produced by NJOY's HEATR module while for photons, |
    |                      |this is tallied from direct photon energy          |
    |                      |deposition. See :ref:`methods_heating`.            |
    +----------------------+---------------------------------------------------+
    |heating-local         |Total nuclear heating in units of eV per source    |
    |                      |particle assuming energy from secondary photons is |
    |                      |deposited locally. Note that this score should only|
    |                      |be used for incident neutrons. See                 |
    |                      |:ref:`methods_heating`.                            |
    +----------------------+---------------------------------------------------+
    |kappa-fission         |The recoverable energy production rate due to      |
    |                      |fission. The recoverable energy is defined as the  |
    |                      |fission product kinetic energy, prompt and delayed |
    |                      |neutron kinetic energies, prompt and delayed       |
    |                      |:math:`\gamma`-ray total energies, and the total   |
    |                      |energy released by the delayed :math:`\beta`       |
    |                      |particles. The neutrino energy does not contribute |
    |                      |to this response. The prompt and delayed           |
    |                      |:math:`\gamma`-rays are assumed to deposit their   |
    |                      |energy locally. Units are eV per source particle.  |
    +----------------------+---------------------------------------------------+
    |fission-q-prompt      |The prompt fission energy production rate. This    |
    |                      |energy comes in the form of fission fragment       |
    |                      |nuclei, prompt neutrons, and prompt                |
    |                      |:math:`\gamma`-rays. This value depends on the     |
    |                      |incident energy and it requires that the nuclear   |
    |                      |data library contains the optional fission energy  |
    |                      |release data. Energy is assumed to be deposited    |
    |                      |locally. Units are eV per source particle.         |
    +----------------------+---------------------------------------------------+
    |fission-q-recoverable |The recoverable fission energy production rate.    |
    |                      |This energy comes in the form of fission fragment  |
    |                      |nuclei, prompt and delayed neutrons, prompt and    |
    |                      |delayed :math:`\gamma`-rays, and delayed           |
    |                      |:math:`\beta`-rays. This tally differs from the    |
    |                      |kappa-fission tally in that it is dependent on     |
    |                      |incident neutron energy and it requires that the   |
    |                      |nuclear data library contains the optional fission |
    |                      |energy release data. Energy is assumed to be       |
    |                      |deposited locally. Units are eV per source         |
    |                      |paticle.                                           |
    +----------------------+---------------------------------------------------+
    |decay-rate            |The delayed-nu-fission-weighted decay rate where   |
    |                      |the decay rate is in units of inverse seconds.     |
    +----------------------+---------------------------------------------------+
    |damage-energy         |Damage energy production in units of eV per source |
    |                      |particle. This corresponds to MT=444 produced by   |
    |                      |NJOY's HEATR module.                               |
    +----------------------+---------------------------------------------------+
    |pulse-height          |The energy deposited by an entire photon's history |
    |                      |(including its progeny). Units are eV per source   |
    |                      |particle. Note that this score can only be combined|
    |                      |with a cell filter and an energy filter.           |
    +----------------------+---------------------------------------------------+
    |ifp-time-numerator    |Adjoint-weighted lifetime of neutron produced by   |
    |                      |fission in units of seconds per source particle.   |
    |                      |This score is used to compute kinetics parameters  |
    |                      |using the iterated fission probability (IFP)       |
    |                      |method.                                            |
    +----------------------+---------------------------------------------------+
    |ifp-beta-numerator    |Adjoint-weighted number of delayed fission events  |
    |                      |in units of number of delayed fission event per    |
    |                      |source particle. This score is used to compute     |
    |                      |kinetics parameters using the iterated fission     |
    |                      |probability (IFP) method.                          |
    +----------------------+---------------------------------------------------+
    |ifp-denominator       |Weights corresponding to the number of fission     |
    |                      |events in units of number of fission event per     |
    |                      |source particle. This score is used to compute     |
    |                      |kinetics parameters using the iterated fission     |
    |                      |probability (IFP) method.                          |
    +----------------------+---------------------------------------------------+

.. _usersguide_virtual_material:

-----------------
Virtual Materials
-----------------

Reaction-rate scores normally include the atom densities of the materials in the
model. Setting :attr:`Tally.multiply_density` to ``False`` instead produces
microscopic, nuclide-wise results. After a simulation, densities from a *virtual
material* can be applied to these results with
:meth:`Tally.apply_virtual_material`. This is useful for determining the
response of a material that is not actually present in the transport model.

For example, the absorbed dose in silicon can be calculated in a region that
does not contain silicon as follows. In this example, ``dose_cell`` identifies
the region of interest that has its volume set in cm³.

::

    silicon = openmc.Material(name='virtual silicon')
    silicon.add_element('Si', 1.0)
    silicon.set_density('g/cm3', 2.329)

    tally = openmc.Tally(name='silicon heating')
    tally.filters = [openmc.CellFilter(dose_cell)]
    tally.scores = ['heating']
    tally.nuclides = silicon.get_nuclides()
    tally.multiply_density = False

    model.tallies = [tally]
    model.run(apply_tally_results=True)

    # Apply the silicon atom densities without collapsing the nuclide axis
    tally.apply_virtual_material(silicon)

    # Sum over nuclides and convert deposited energy to dose
    eV_per_source = tally.mean.sum()
    J_per_source = eV_per_source * openmc.data.JOULE_PER_EV
    J_per_cm3_source = J_per_source / dose_cell.volume
    kg_per_cm3 = silicon.get_mass_density() * 1.0e-3
    Gy_per_source = J_per_cm3_source / kg_per_cm3

Before the virtual material is applied, the nuclide-wise heating results have
units of eV-b-cm/(atom-source). Multiplication by atom densities in atom/b-cm
gives eV/source while retaining the nuclide dimension. The explicit sum above
combines the nuclide contributions. Finally, dividing the deposited energy in
joules by the virtual silicon mass in kilograms gives absorbed dose in
Gy/source.

The tally must contain individual nuclide bins; a ``'total'`` nuclide bin cannot
be used because a material density cannot be assigned to its constituent
nuclides. Tally nuclides that do not occur in the virtual material are treated
as having zero atom density.

.. _usersguide_tally_normalization:

------------------------------
Normalization of Tally Results
------------------------------

As described in :ref:`usersguide_scores`, all tally scores are normalized per
source particle simulated. However, for analysis of a given system, we usually
want tally scores in a more natural unit. For example, neutron flux is often
reported in units of particles/cm\ :sup:`2`\ -s. For a fixed source simulation,
it is usually straightforward to convert units if the source rate is known. For
example, if the system being modeled includes a source that is emitting 10\
:sup:`4` neutrons per second, the tally results just need to be multipled by 10\
:sup:`4`. This can either be done manually or using the
:attr:`openmc.SourceBase.strength` attribute.

For a :math:`k`\ -eigenvalue calculation, normalizing tally results is not as
simple because the source rate is not actually known. Instead, we typically know
the system power, :math:`P`, which represents how much energy is deposited per
unit time. Most of this energy originates from fission, but a small percentage
also results from other reactions (e.g., photons emitted from :math:`(n,\gamma)`
reactions). The most rigorous method to normalize tally results is to run a
coupled neutron-photon calculation and tally the ``heating`` score over the
entire system. This score provides the heating rate in units of [eV/source],
which we'll denote :math:`H`. Then, calculate the heating rate in J/source as

.. math::

    H' = 1.602\times10^{-19} \left [ \frac{\text{J}}{\text{eV}} \right ] \cdot H
    \left [\frac{\text{eV}}{\text{source}} \right ].

Dividing the power by the observed heating rate then gives us a normalization
factor that can be applied to other tallies:

.. math::

    f = \frac{P}{H'} = \frac{[\text{J}/\text{s}]}{[\text{J}/\text{source}]} =
    \left [ \frac{\text{source}}{\text{s}} \right ].

Multiplying by the normalization factor and dividing by volume, we can then get
the flux in typical units:

.. math::

    \phi' = \frac{f\phi}{V} =
    \frac{[\text{source}/\text{s}][\text{particle-cm}/\text{source}]}
    {[\text{cm}^3]} = \left [\frac{\text{particle}}{\text{cm}^2\cdot\text{s}}
    \right ]

There are several slight variations on this procedure:

- Run a neutron-only calculation and estimate the total heating using the
  ``heating-local`` score (this requires that your nuclear data has local
  heating data available, such as in the official data library at
  https://openmc.org. See :ref:`methods_heating` for more information.)
- Run a neutron-only calculation and use the ``kappa-fission`` or
  ``fission-q-recoverable`` scores along with an estimate of the extra heating
  due to neutron capture reactions.
- Calculate the overall fission rate and then used a fixed Q value to estimate
  the heating rate.

Note that the only difference between these and the above procedures is in how
:math:`H'` is estimated.
