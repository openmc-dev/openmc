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

You can also list 'all' as a nuclide which will give you a separate reaction
rate for every nuclide in the model.

The following tables show all valid scores:

.. table:: **Flux scores: units are particle-cm per source particle.**

    +----------------------+---------------------------------------------------+
    |Score                 | Description                                       |
    +======================+===================================================+
    |flux                  |Total flux.                                        |
    +----------------------+---------------------------------------------------+

.. table:: **Reaction scores: units are reactions per source particle.**

    +----------------------+---------------------------------------------------+
    |Score                 | Description                                       |
    +======================+===================================================+
    |absorption            |Total absorption rate. This accounts for all       |
    |                      |reactions which do not produce secondary neutrons  |
    |                      |as well as fission.                                |
    +----------------------+---------------------------------------------------+
    |elastic               |Elastic scattering reaction rate.                  |
    +----------------------+---------------------------------------------------+
    |fission               |Total fission reaction rate.                       |
    +----------------------+---------------------------------------------------+
    |scatter               |Total scattering rate.                             |
    +----------------------+---------------------------------------------------+
    |total                 |Total reaction rate.                               |
    +----------------------+---------------------------------------------------+
    |(n,2nd)               |(n,2nd) reaction rate.                             |
    +----------------------+---------------------------------------------------+
    |(n,2n)                |(n,2n) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,3n)                |(n,3n) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,na)                |(n,n\ :math:`\alpha`\ ) reaction rate.             |
    +----------------------+---------------------------------------------------+
    |(n,n3a)               |(n,n3\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,2na)               |(n,2n\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,3na)               |(n,3n\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,np)                |(n,np) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,n2a)               |(n,n2\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,2n2a)              |(n,2n2\ :math:`\alpha`\ ) reaction rate.           |
    +----------------------+---------------------------------------------------+
    |(n,nd)                |(n,nd) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,nt)                |(n,nt) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,nHe-3)             |(n,n\ :sup:`3`\ He) reaction rate.                 |
    +----------------------+---------------------------------------------------+
    |(n,nd2a)              |(n,nd2\ :math:`\alpha`\ ) reaction rate.           |
    +----------------------+---------------------------------------------------+
    |(n,nt2a)              |(n,nt2\ :math:`\alpha`\ ) reaction rate.           |
    +----------------------+---------------------------------------------------+
    |(n,4n)                |(n,4n) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,2np)               |(n,2np) reaction rate.                             |
    +----------------------+---------------------------------------------------+
    |(n,3np)               |(n,3np) reaction rate.                             |
    +----------------------+---------------------------------------------------+
    |(n,n2p)               |(n,n2p) reaction rate.                             |
    +----------------------+---------------------------------------------------+
    |(n,n*X*)              |Level inelastic scattering reaction rate. The *X*  |
    |                      |indicates what which inelastic level, e.g., (n,n3) |
    |                      |is third-level inelastic scattering.               |
    +----------------------+---------------------------------------------------+
    |(n,nc)                |Continuum level inelastic scattering reaction rate.|
    +----------------------+---------------------------------------------------+
    |(n,gamma)             |Radiative capture reaction rate.                   |
    +----------------------+---------------------------------------------------+
    |(n,p)                 |(n,p) reaction rate.                               |
    +----------------------+---------------------------------------------------+
    |(n,d)                 |(n,d) reaction rate.                               |
    +----------------------+---------------------------------------------------+
    |(n,t)                 |(n,t) reaction rate.                               |
    +----------------------+---------------------------------------------------+
    |(n,3He)               |(n,\ :sup:`3`\ He) reaction rate.                  |
    +----------------------+---------------------------------------------------+
    |(n,a)                 |(n,\ :math:`\alpha`\ ) reaction rate.              |
    +----------------------+---------------------------------------------------+
    |(n,2a)                |(n,2\ :math:`\alpha`\ ) reaction rate.             |
    +----------------------+---------------------------------------------------+
    |(n,3a)                |(n,3\ :math:`\alpha`\ ) reaction rate.             |
    +----------------------+---------------------------------------------------+
    |(n,2p)                |(n,2p) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,pa)                |(n,p\ :math:`\alpha`\ ) reaction rate.             |
    +----------------------+---------------------------------------------------+
    |(n,t2a)               |(n,t2\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,d2a)               |(n,d2\ :math:`\alpha`\ ) reaction rate.            |
    +----------------------+---------------------------------------------------+
    |(n,pd)                |(n,pd) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,pt)                |(n,pt) reaction rate.                              |
    +----------------------+---------------------------------------------------+
    |(n,da)                |(n,d\ :math:`\alpha`\ ) reaction rate.             |
    +----------------------+---------------------------------------------------+
    |*Arbitrary integer*   |An arbitrary integer is interpreted to mean the    |
    |                      |reaction rate for a reaction with a given ENDF MT  |
    |                      |number.                                            |
    +----------------------+---------------------------------------------------+

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
    |current               |Used in combination with a meshsurface filter:     |
    |                      |Partial currents on the boundaries of each cell in |
    |                      |a mesh. It may not be used in conjunction with any |
    |                      |other score. Only energy and mesh filters may be   |
    |                      |used.                                              |
    |                      |Used in combination with a surface filter:         |
    |                      |Net currents on any surface previously defined in  |
    |                      |the geometry. It may be used along with any other  |
    |                      |filter, except meshsurface filters.                |
    |                      |Surfaces can alternatively be defined with cell    |
    |                      |from and cell filters thereby resulting in tallying|
    |                      |partial currents.                                  |
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
    |                      |this is tallied from either direct photon energy   |
    |                      |deposition (analog estimator) or pre-generated     |
    |                      |photon heating number. See :ref:`methods_heating`  |
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
