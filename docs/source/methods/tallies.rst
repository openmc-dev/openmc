.. _methods_tallies:

=======
Tallies
=======

------------------
Filters and Scores
------------------

The tally capability in OpenMC takes a similar philosophy as that employed in
the MC21_ Monte Carlo code to give maximum flexibility in specifying tallies
while still maintaining scalability. Any tally in a Monte Carlo simulation can
be written in the following form:

.. math::
    :label: tally-integral

    X = \underbrace{\int d\mathbf{r} \int d\mathbf{\Omega} \int
    dE}_{\text{filters}} \underbrace{f(\mathbf{r}, \mathbf{\Omega},
    E)}_{\text{scores}} \psi (\mathbf{r}, \mathbf{\Omega}, E)


A user can specify one or more filters which identify which regions of phase
space should score to a given tally (the limits of integration as shown in
equation :eq:`tally-integral`) as well as the scoring function (:math:`f` in
equation :eq:`tally-integral`). For example, if the desired tally was the
:math:`(n,\gamma)` reaction rate in a fuel pin, the filter would specify the
cell which contains the fuel pin and the scoring function would be the radiative
capture macroscopic cross section. The following quantities can be scored in
OpenMC: flux, total reaction rate, scattering reaction rate, neutron production
from scattering, higher scattering moments, (n,xn) reaction rates, absorption
reaction rate, fission reaction rate, neutron production rate from fission, and
surface currents. The following variables can be used as filters: universe,
material, cell, birth cell, surface, mesh, pre-collision energy, and
post-collision energy.

With filters for pre- and post-collision energy and scoring functions for
scattering and fission production, it is possible to use OpenMC to generate
cross sections with user-defined group structures. These multigroup cross
sections can subsequently be used in deterministic solvers such as coarse-mesh
finite difference (CMFD) diffusion.

------------------------------
Using Maps for Filter-Matching
------------------------------

Some Monte Carlo codes suffer severe performance penalties when tallying a large
number of quantities. Care must be taken to ensure that a tally system scales
well with the total number of tally bins. In OpenMC, a mapping technique is used
that allows for a fast determination of what tally/bin combinations need to be
scored to a given particle's phase space coordinates. For each discrete filter
variable, a list is stored that contains the tally/bin combinations that could
be scored to for each value of the filter variable. If a particle is in cell
:math:`n`, the mapping would identify what tally/bin combinations specify cell
:math:`n` for the cell filter variable. In this manner, it is not necessary to
check the phase space variables against each tally. Note that this technique
only applies to discrete filter variables and cannot be applied to energy
bins. For energy filters, it is necessary to perform a binary search on the
specified energy grid.

-----------------------------------------
Volume-Integrated Flux and Reaction Rates
-----------------------------------------

One quantity we may wish to compute during the course of a Monte Carlo
simulation is the flux or a reaction rate integrated over a finite volume. The
volume may be a particular cell, a collection of cells, or the entire
geometry. There are various methods by which we can estimate reaction rates

Analog Estimator
----------------

The analog estimator is the simplest type of estimator for reaction rates. The
basic idea is that we simply count the number of actual reactions that take
place and use that as our estimate for the reaction rate. This can be written
mathematically as

.. math::
    :label: analog-estimator

    R_x = \frac{1}{W} \sum_{i \in A} w_i

where :math:`R_x` is the reaction rate for reaction :math:`x`, :math:`i` denotes
an index for each event, :math:`A` is the set of all events resulting in
reaction :math:`x`, and :math:`W` is the total starting weight of the particles,
and :math:`w_i` is the pre-collision weight of the particle as it enters event
:math:`i`. One should note that equation :eq:`analog-estimator` is
volume-integrated so if we want a volume-averaged quantity, we need to divided
by the volume of the region of integration.

Collision Estimator
-------------------

While the analog estimator is conceptually very simple and easy to implement, it
can suffer higher variance due to the fact low probability events will not occur
often enough to get good statistics if they are being tallied. Thus, it is
desirable to use a different estimator that allows us to score to the tally more
often. One such estimator is the collision estimator. Instead of tallying a
reaction only when it happens, the idea is to make a contribution to the tally
at every collision.

We can start by writing a formula for the collision estimate of the flux. Since
:math:`R = \Sigma_t \phi` where :math:`R` is the total reaction rate,
:math:`\Sigma_t` is the total macroscopic cross section, and :math:`\phi` is the
scalar flux, it stands to reason that we can estimate the flux by taking an
estimate of the total reaction rate and dividing it by the total macroscopic
cross section. This gives us the following formula:

.. math::
    :label: collision-estimator-flux

    \phi = \frac{1}{W} \sum_{i \in C} \frac{w_i}{\Sigma_t (E_i)}

where :math:`W` is again the total starting weight of the particles, :math:`C`
is the set of all events resulting in a collision with a nucleus, and
:math:`\Sigma_t (E)` is the total macroscopic cross section of the target
material at the incoming energy of the particle :math:`E_i`.

If we multiply both sides of equation :eq:`collision-estimator-flux` by the
macroscopic cross section for some reaction :math:`x`, then we get the collision
estimate for the reaction rate for that reaction:

.. math::
    :label: collision-estimator

    R_x = \frac{1}{W} \sum_{i \in C} \frac{w_i \Sigma_x (E_i)}{\Sigma_t (E_i)}

where :math:`\Sigma_x (E_i)` is the macroscopic cross section for reaction
:math:`x` at the incoming energy of the particle :math:`E_i`. In comparison to
equation :eq:`analog-estimator`, we see that the collision estimate will result
in a tally with a larger number of events that score to it with smaller
contributions (since we have multiplied it by :math:`\Sigma_x / \Sigma_t`).

Track-length Estimator
----------------------

One other method we can use to increase the number of events that scores to
tallies is to use an estimator the scores contributions to a tally at every
track for the particle rather than every collision. This is known as a
track-length estimator, sometimes also called a path-length estimator. We first
start with an expression for the volume integrated flux, which can be written as

.. math::
    :label: flux-integrated

    V \phi = \int d\mathbf{r} \int dE \int d\mathbf{\Omega} \int dt \,
    \psi(\mathbf{r}, \mathbf{\hat{\Omega}}, E, t).

where :math:`V` is the volume, :math:`\psi` is the angular flux,
:math:`\mathbf{r}` is the position of the particle, :math:`\mathbf{\hat{\Omega}}`
is the direction of the particle, :math:`E` is the energy of the particle, and
:math:`t` is the time. By noting that :math:`\psi(\mathbf{r},
\mathbf{\hat{\Omega}}, E, t) = v n(\mathbf{r}, \mathbf{\hat{\Omega}}, E, t)`
where :math:`n` is the angular neutron density, we can rewrite equation
:eq:`flux-integrated` as

.. math::
    :label: flux-integrated-2

    V \phi = \int d\mathbf{r} \int dE \int dt v \int d\mathbf{\Omega} \, n(\mathbf{r},
    \mathbf{\hat{\Omega}}, E, t))

Using the relations :math:`N(\mathbf{r}, E, t) = \int d\mathbf{\Omega}
n(\mathbf{r}, \mathbf{\hat{\Omega}}, E, t)` and :math:`d\ell = v \, dt` where
:math:`d\ell` is the differential unit of track length, we then obtain

.. math::
    :label: track-length-integral

    V \phi = \int d\mathbf{r} \int dE \int d\ell N(\mathbf{r}, E, t)

Equation :eq:`track-length-integral` indicates that we can use the length of a
particle's trajectory as an estimate for the flux, i.e. the track-length
estimator of the flux would be

.. math::
    :label: track-length-flux

    \phi = \frac{1}{W} \sum_{i \in T} w_i \ell_i

where :math:`T` is the set of all the particle's trajectories within the desired
volume and :math:`\ell_i` is the length of the :math:`i`-th trajectory. In the
same vein as equation :eq:`collision-estimator`, the track-length estimate of a
reaction rate is found by multiplying equation :eq:`track-length-flux` by a
macroscopic reaction cross section:

.. math::
    :label: track-length-estimator

    R_x = \frac{1}{W} \sum_{i \in T} w_i \ell_i \Sigma_x (E_i)

One important fact to take into consideration is that the use of a track-length
estimator precludes us from using any filter that requires knowledge of the
particle's state following a collision because by definition, it will not have
had a collision at every event. Thus, for tallies with outgoing-energy filters
(which require the post-collision energy) or for tallies of scattering moments
(which require the scattering cosine), we must use an analog estimator.

---------------
Surface Current
---------------

.. _MC21: http://www.osti.gov/bridge/servlets/purl/903083-HT5p1o/903083.pdf
