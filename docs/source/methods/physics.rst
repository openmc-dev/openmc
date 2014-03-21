.. _methods_physics:

=======
Physics
=======

-----------------------------------
Sampling Distance to Next Collision
-----------------------------------

As a particle travels through a homogeneous material, the probability
distribution function for the distance to its next collision :math:`\ell` is

.. math::
    :label: distance-pdf

    p(\ell) d\ell = \Sigma_t e^{-\Sigma_t \ell} d\ell

where :math:`\Sigma_t` is the total macroscopic cross section of the
material. Equation :eq:`distance-pdf` tells us that the further the distance is
to the next collision, the less likely the particle will travel that
distance. In order to sample the probability distribution function, we first
need to convert it to a cumulative distribution function

.. math::
    :label: distance-cdf

    \int_0^{\ell} d\ell' p(\ell') = \int_0^{\ell} d\ell' \Sigma_t e^{-\Sigma_t
    \ell'} = 1 - e^{-\Sigma_t \ell}.

By setting the cumulative distribution function equal to :math:`\xi`, a random
number on the unit interval, and solving for the distance :math:`\ell`, we
obtain a formula for sampling the distance to next collision:

.. math::
    :label: sample-distance-1

    \ell = -\frac{\ln (1 - \xi)}{\Sigma_t}.

Since :math:`\xi` is uniformly distributed on :math:`[0,1)`, this implies that
:math:`1 - \xi` is also uniformly distributed on :math:`[0,1)` as well. Thus,
the formula usually used to calculate the distance to next collision is

.. math::
    :label: sample-distance-2

    \ell = -\frac{\ln \xi}{\Sigma_t}

----------------------------------------------------
:math:`(n,\gamma)` and Other Disappearance Reactions
----------------------------------------------------

All absorption reactions other than fission do not produce any secondary
neutrons. As a result, these are the easiest type of reactions to handle. When a
collision occurs, the first step is to sample a nuclide within a material. Once
the nuclide has been sampled, then a specific reaction for that nuclide is
sampled. Since the total absorption cross section is pre-calculated at the
beginning of a simulation, the first step in sampling a reaction is to determine
whether a "disappearance" reaction occurs where no secondary neutrons are
produced. This is done by sampling a random number :math:`\xi` on the interval
:math:`[0,1)` and checking whether

.. math::
    :label: disappearance

    \xi \sigma_t (E) < \sigma_a (E) - \sigma_f (E)

where :math:`\sigma_t` is the total cross section, :math:`\sigma_a` is the
absorption cross section (this includes fission), and :math:`\sigma_f` is the
total fission cross section. If this condition is met, then the neutron is
killed and we proceed to simulate the next neutron from the source bank.

No secondary particles from disappearance reactions such as photons or
alpha-particles are produced or tracked. To truly capture the affects of gamma
heating in a problem, it would be necessary to explicitly track photons
originating from :math:`(n,\gamma)` and other reactions.

------------------
Elastic Scattering
------------------

Elastic scattering refers to the process by which a neutron scatters off a
nucleus and does not leave it in an excited. It is referred to as "elastic"
because in the center-of-mass system, the neutron does not actually lose
energy. However, in lab coordinates, the neutron does indeed lose
energy. Elastic scattering can be treated exactly in a Monte Carlo code thanks
to its simplicity.

Let us discuss how OpenMC handles two-body elastic scattering kinematics. The
first step is to determine whether the target nucleus has any associated
motion. Above a certain energy threshold (400 kT by default), all scattering is
assumed to take place with the target at rest. Below this threshold though, we
must account for the thermal motion of the target nucleus. Methods to sample the
velocity of the target nucleus are described later in section
:ref:`freegas`. For the time being, let us assume that we have sampled the
target velocity :math:`\mathbf{v}_t`. The velocity of the center-of-mass system
is calculated as

.. math::
    :label: velocity-com

    \mathbf{v}_{cm} = \frac{\mathbf{v}_n + A \mathbf{v}_t}{A + 1}

where :math:`\mathbf{v}_n` is the velocity of the neutron and :math:`A` is the
atomic mass of the target nucleus measured in neutron masses (commonly referred
to as the *atomic weight ratio*). With the velocity of the center-of-mass
calculated, we can then determine the neutron's velocity in the center-of-mass
system:

.. math::
    :label: velocity-neutron-com

    \mathbf{V}_n = \mathbf{v}_n - \mathbf{v}_{cm}

where we have used uppercase :math:`\mathbf{V}` to denote the center-of-mass
system. The direction of the neutron in the center-of-mass system is

.. math::
    :label: angle-neutron-com

    \mathbf{\Omega}_n = \frac{\mathbf{V}_n}{|| \mathbf{V}_n ||}.

At low energies, elastic scattering will be isotropic in the center-of-mass
system, but for higher energies, there may be p-wave and higher order scattering
that leads to anisotropic scattering. Thus, in general, we need to sample a
cosine of the scattering angle which we will refer to as :math:`\mu`. For
elastic scattering, the secondary angle distribution is always given in the
center-of-mass system and is sampled according to the procedure outlined in
:ref:`sample-angle`. After the cosine of the angle of scattering has been
sampled, we need to determine the neutron's new direction
:math:`\mathbf{\Omega}'_n` in the center-of-mass system. This is done with the
procedure in :ref:`transform-coordinates`. The new direction is multiplied by
the speed of the neutron in the center-of-mass system to obtain the new velocity
vector in the center-of-mass:

.. math::
    :label: velocity-neutron-com-2

    \mathbf{V}'_n = || \mathbf{V}_n || \mathbf{\Omega}'_n.

Finally, we transform the velocity in the center-of-mass system back to lab
coordinates:

.. math::
    :label: velocity-neutron-lab

    \mathbf{v}'_n = \mathbf{V}'_n + \mathbf{v}_{cm}

In OpenMC, the angle and energy of the neutron are stored rather than the
velocity vector itself, so the post-collision angle and energy can be inferred
from the post-collision velocity of the neutron in the lab system.

For tallies that require the scattering cosine, it is important to store the
scattering cosine in the lab system. If we know the scattering cosine in the
center-of-mass, the scattering cosine in the lab system can be calculated as

.. math::
    :label: cosine-lab

    \mu_{lab} = \frac{1 + A\mu}{\sqrt{A^2 + 2A\mu + 1}}.

However, equation :eq:`cosine-lab` is only valid if the target was at rest. When
the target nucleus does have thermal motion, the cosine of the scattering angle
can be determined by simply taking the dot product of the neutron's initial and
final direction in the lab system.

.. _inelastic-scatter:

--------------------
Inelastic Scattering
--------------------

The major algorithms for inelastic scattering were described in previous
sections. First, a scattering cosine is sampled using the algorithms in
:ref:`sample-angle`. Then an outgoing energy is sampled using the algorithms in
:ref:`sample-energy`. If the outgoing energy and scattering cosine were given in
the center-of-mass system, they are transformed to laboratory coordinates using
the algorithm described in :ref:`transform-coordinates`. Finally, the direction
of the particle is changed also using the procedure in
:ref:`transform-coordinates`.

Although inelastic scattering leaves the target nucleus in an excited state, no
secondary photons from nuclear de-excitation are tracked in OpenMC.

------------------------
:math:`(n,xn)` Reactions
------------------------

These types of reactions are just treated as inelastic scattering and as such
are subject to the same procedure as described in
:ref:`inelastic-scatter`. Rather than tracking multiple secondary neutrons, the
weight of the outgoing neutron is multiplied by the number of secondary
neutrons, e.g. for :math:`(n,2n)`, only one outgoing neutron is tracked but its
weight is doubled.

.. _fission:

-------
Fission
-------

While fission is normally considered an absorption reaction, as far as it
concerns a Monte Carlo simulation it actually bears more similarities to
inelastic scattering since fission results in secondary neutrons in the exit
channel. Other absorption reactions like :math:`(n,\gamma)` or
:math:`(n,\alpha)`, on the contrary, produce no neutrons. There are a few other
idiosyncrasies in treating fission. In an eigenvalue calculation, secondary
neutrons from fission are only "banked" for use in the next generation rather
than being tracked as secondary neutrons from elastic and inelastic scattering
would be. On top of this, fission is sometimes broken into first-chance fission,
second-chance fission, etc. An ACE table either lists the partial fission
reactions with secondary energy distributions for each one, or a total fission
reaction with a single secondary energy distribution.

When a fission reaction is sampled in OpenMC (either total fission or, if data
exists, first- or second-chance fission), the following algorithm is used to
create and store fission sites for the following generation. First, the average
number of prompt and delayed neutrons must be determined to decide whether the
secondary neutrons will be prompt or delayed. This is important because delayed
neutrons have a markedly different spectrum from prompt neutrons, one that has a
lower average energy of emission. The total number of neutrons emitted
:math:`\nu_t` is given as a function of incident energy in the ACE format. Two
representations exist for :math:`\nu_t`. The first is a polynomial of order
:math:`N` with coefficients :math:`c_0,c_1,\dots,c_N`. If :math:`\nu_t` has this
format, we can evaluate it at incoming energy :math:`E` by using the equation

.. math::
    :label: nu-polynomial

    \nu_t (E) = \sum_{i = 0}^N c_i E^i.

The other representation is just a tabulated function with a specified
interpolation law. The number of prompt neutrons released per fission event
:math:`\nu_p` is also given as a function of incident energy and can be
specified in a polynomial or tabular format. The number of delayed neutrons
released per fission event :math:`\nu_d` can only be specified in a tabular
format. In practice, we only need to determine :math:`nu_t` and
:math:`nu_d`. Once these have been determined, we can calculated the delayed
neutron fraction

.. math::
    :label: beta

    \beta = \frac{\nu_d}{\nu_t}.

We then need to determine how many total neutrons should be emitted from
fission. If no survival biasing is being used, then the number of neutrons
emitted is

.. math::
    :label: fission-neutrons

    \nu = \frac{w \nu_t}{k_{eff}}

where :math:`w` is the statistical weight and :math:`k_{eff}` is the effective
multiplication factor from the previous generation. The number of neutrons
produced is biased in this manner so that the expected number of fission
neutrons produced is the number of source particles that we started with in the
generation. Since :math:`\nu` is not an integer, we use the following procedure
to obtain an integral number of fission neutrons to produce. If :math:`\xi >
\nu - \lfloor \nu \rfloor`, then we produce :math:`\lfloor \nu \rfloor`
neutrons. Otherwise, we produce :math:`\lfloor \nu \rfloor + 1` neutrons. Then,
for each fission site produced, we sample the outgoing angle and energy
according to the algorithms given in :ref:`sample-angle` and
:ref:`sample-energy` respectively. If the neutron is to be born delayed, then
there is an extra step of sampling a delayed neutron precursor group since they
each have an associated secondary energy distribution.

The sampled outgoing angle and energy of fission neutrons along with the
position of the collision site are stored in an array called the fission
bank. In a subsequent generation, these fission bank sites are used as starting
source sites.

-----------------------------------------
Secondary Angles and Energy Distributions
-----------------------------------------

For any reactions with secondary neutrons, it is necessary to sample secondary
angle and energy distributions. This includes elastic and inelastic scattering,
fission, and :math:`(n,xn)` reactions. In some cases, the angle and energy
distributions may be specified separately, and in other cases, they may be
specified as a correlated angle-energy distribution. In the following sections,
we will outline the methods used to sample secondary distributions as well as
how they are used to modify the state of a particle.

.. _sample-angle:

Sampling Secondary Angle Distributions
--------------------------------------

For elastic scattering, it is only necessary to specific a secondary angle
distribution since the outgoing energy can be determined analytically. Other
reactions may also have separate secondary angle and secondary energy
distributions that are uncorrelated. In these cases, the secondary angle
distribution is represented as either

- An Isotropic angular distribution,
- An equiprobable distribution with 32 bins, or
- A tabular distribution.

Isotropic Angular Distribution
++++++++++++++++++++++++++++++

In the first case, no data needs to be stored on the ACE table, and the cosine
of the scattering angle is simply calculated as

.. math::
    :label: isotropic-angle

    \mu = 2\xi - 1

where :math:`\mu` is the cosine of the scattering angle and :math:`\xi` is a
random number sampled uniformly on :math:`[0,1)`.

Equiprobable Angle Bin Distribution
+++++++++++++++++++++++++++++++++++

For a 32 equiprobable bin distribution, we select a random number :math:`\xi` to
sample a cosine bin :math:`i` such that

.. math::
    :label: equiprobable-bin

    i = 1 + \lfloor 32\xi \rfloor.

The same random number can then also be used to interpolate between neighboring
:math:`\mu` values to get the final scattering cosine:

.. math::
    :label: equiprobable-cosine

    \mu = \mu_i + (32\xi - i) (\mu_{i+1} - \mu_i)

where :math:`\mu_i` is the :math:`i`-th scattering cosine.

.. _angle-tabular:

Tabular Angular Distribution
++++++++++++++++++++++++++++

As the `MCNP Manual`_ points out, using an equiprobable bin distribution works
well for high-probability regions of the scattering cosine probability, but for
low-probability regions it is not very accurate. Thus, a more accurate method is
to represent the scattering cosine with a tabular distribution. In this case, we
have a table of cosines and their corresponding values for a probability
distribution function and cumulative distribution function. For each incoming
neutron energy :math:`E_i`, let us call :math:`p_{i,j}` the j-th value in the
probability distribution function and :math:`c_{i,j}` the j-th value in the
cumulative distribution function. We first find the interpolation factor on the
incoming energy grid:

.. math::
    :label: interpolation-factor

    f = \frac{E - E_i}{E_{i+1} - E_i}

where :math:`E` is the incoming energy of the particle. Then, statistical
interpolation is performed to choose between using the cosines and distribution
functions corresponding to energy :math:`E_i` and :math:`E_{i+1}`. Let
:math:`\ell` be the chosen table where :math:`\ell = i` if :math:`\xi_1 > f` and
:math:`\ell = i + 1` otherwise, where :math:`\xi_1` is a random number. Another
random number :math:`\xi_2` is used to sample a scattering cosine bin :math:`j`
using the cumulative distribution function:

.. math::
    :label: sample-cdf

    c_{\ell,j} < \xi_2 < c_{\ell,j+1}

The final scattering cosine will depend on whether histogram or linear-linear
interpolation is used. In general, we can write the cumulative distribution
function as

.. math::
    :label: cdf

    c(\mu) = \int_{-1}^\mu p(\mu') d\mu'

where :math:`c(\mu)` is the cumulative distribution function and :math:`p(\mu)`
is the probability distribution function. Since we know that
:math:`c(\mu_{\ell,j}) = c_{\ell,j}`, this implies that for :math:`\mu >
\mu_{\ell,j}`,

.. math::
    :label: cdf-2

    c(\mu) = c_{\ell,j} + \int_{\mu_{\ell,j}}^{\mu} p(\mu') d\mu'

For histogram interpolation, we have that :math:`p(\mu') = p_{\ell,j}` for
:math:`\mu_{\ell,j} \le \mu' < \mu_{\ell,j+1}`. Thus, after integrating
:eq:`cdf-2` we have that

.. math::
    :label: cumulative-dist-histogram

    c(\mu) = c_{\ell,j} + (\mu - \mu_{\ell,j}) p_{\ell,j} = \xi_2

Solving for the scattering cosine, we obtain the final form for histogram
interpolation:

.. math::
    :label: cosine-histogram

    \mu = \mu_{\ell,j} + \frac{\xi_2 - c_{\ell,j}}{p_{\ell,j}}.

For linear-linear interpolation, we represent the function :math:`p(\mu')` as a
first-order polynomial in :math:`\mu'`. If we interpolate between successive
values on the probability distribution function, we know that

.. math::
    :label: pdf-interpolation

    p(\mu') - p_{\ell,j} = \frac{p_{\ell,j+1} - p_{\ell,j}}{\mu_{\ell,j+1} -
    \mu_{\ell,j}} (\mu' - \mu_{\ell,j})

Solving for :math:`p(\mu')` in equation :eq:`pdf-interpolation` and inserting it
into equation :eq:`cdf-2`, we obtain

.. math::
    :label: cdf-linlin

    c(\mu) = c_{\ell,j} + \int_{\mu_{\ell,j}}^{\mu} \left [ \frac{p_{\ell,j+1} -
    p_{\ell,j}}{\mu_{\ell,j+1} - \mu_{\ell,j}} (\mu' - \mu_{\ell,j}) +
    p_{\ell,j} \right ] d\mu'.

Let us now make a change of variables using

.. math::
    :label: introduce-eta

    \eta = \frac{p_{\ell,j+1} - p_{\ell,j}}{\mu_{\ell,j+1} - \mu_{\ell,j}}
    (\mu' - \mu_{\ell,j}) + p_{\ell,j}.

Equation :eq:`cdf-linlin` then becomes

.. math::
    :label: cdf-linlin-eta

    c(\mu) = c_{\ell,j} + \frac{1}{m} \int_{p_{\ell,j}}^{m(\mu - \mu_{\ell,j}) +
    p_{\ell,j}} \eta \, d\eta

where we have used

.. math::
    :label: slope

    m = \frac{p_{\ell,j+1} - p_{\ell,j}}{\mu_{\ell,j+1} - \mu_{\ell,j}}.

Integrating equation :eq:`cdf-linlin-eta`, we have

.. math::
    :label: cdf-linlin-integrated

    c(\mu) = c_{\ell,j} + \frac{1}{2m} \left ( \left [ m (\mu - \mu_{\ell,j} ) +
    p_{\ell,j} \right ]^2 - p_{\ell,j}^2 \right ) = \xi_2

Solving for :math:`\mu`, we have the final form for the scattering cosine using
linear-linear interpolation:

.. math::
    :label: cosine-linlin

    \mu = \mu_{\ell,j} + \frac{1}{m} \left ( \sqrt{p_{\ell,j}^2 + 2 m (\xi_2 -
    c_{\ell,j} )} - p_{\ell,j} \right )

.. _sample-energy:

Sampling Secondary Energy and Correlated Angle/Energy Distributions
-------------------------------------------------------------------

For a reaction with secondary neutrons, it is necessary to determine the
outgoing energy of the neutrons. For any reaction other than elastic scattering,
the outgoing energy must be determined based on tabulated or parameterized
data. The `ENDF-6 Format`_ specifies a variety of ways that the secondary energy
distribution can be represented. ENDF File 5 contains uncorrelated energy
distribution where ENDF File 6 contains correlated energy-angle
distributions. The ACE format specifies its own representations based loosely on
the formats given in ENDF-6. In this section, we will describe how the outgoing
energy of secondary particles is determined based on each ACE law.

One of the subtleties in the ACE format is the fact that a single reaction can
have multiple secondary energy distributions. This is mainly useful for
reactions with multiple neutrons in the exit channel such as :math:`(n,2n)` or
:math:`(n,3n)`. In these types of reactions, each neutron is emitted
corresponding to a different excitation level of the compound nucleus, and thus
in general the neutrons will originate from different energy distributions. If
multiple energy distributions are present, they are assigned probabilities that
can then be used to randomly select one.

Once a secondary energy distribution has been sampled, the procedure for
determining the outgoing energy will depend on which ACE law has been specified
for the data.

.. _ace-law-1:

ACE Law 1 - Tabular Equiprobable Energy Bins
++++++++++++++++++++++++++++++++++++++++++++

In the tabular equiprobable bin representation, an array of equiprobable
outgoing energy bins is given for a number of incident energies. While the
representation itself is simple, the complexity lies in how one interpolates
between incident as well as outgoing energies on such a table. If one performs
simple interpolation between tables for neighboring incident energies, it is
possible that the resulting energies would violate laws governing the
kinematics, i.e. the outgoing energy may be outside the range of available
energy in the reaction.

To avoid this situation, the accepted practice is to use a process known as
scaled interpolation [Doyas]_. First, we find the tabulated incident energies
which bound the actual incoming energy of the particle, i.e. find :math:`i` such
that :math:`E_i < E < E_{i+1}` and calculate the interpolation factor :math:`f`
via :eq:`interpolation-factor`. Then, we interpolate between the minimum and
maximum energies of the outgoing energy distributions corresponding to
:math:`E_i` and :math:`E_{i+1}`:

.. math::
    :label: ace-law-1-minmax

    E_{min} = E_{i,1} + f ( E_{i+1,1} - E_i ) \\
    E_{max} = E_{i,M} + f ( E_{i+1,M} - E_M )

where :math:`E_{min}` and :math:`E_{max}` are the minimum and maximum outgoing
energies of a scaled distribution, :math:`E_{i,j}` is the j-th outgoing energy
corresponding to the incoming energy :math:`E_i`, and :math:`M` is the number of
outgoing energy bins. Next, statistical interpolation is performed to choose
between using the outgoing energy distributions corresponding to energy
:math:`E_i` and :math:`E_{i+1}`. Let :math:`\ell` be the chosen table where
:math:`\ell = i` if :math:`\xi_1 > f` and :math:`\ell = i + 1` otherwise, and
:math:`\xi_1` is a random number. Now, we randomly sample an equiprobable
outgoing energy bin :math:`j` and interpolate between successive values on the
outgoing energy distribution:

.. math::
    :label: ace-law-1-intermediate

    \hat{E} = E_{\ell,j} + \xi_2 (E_{\ell,j+1} - E_{\ell,j})

where :math:`\xi_2` is a random number sampled uniformly on :math:`[0,1)`. Since
this outgoing energy may violate reaction kinematics, we then scale it to the
minimum and maximum energies we calculated earlier to get the final outgoing
energy:

.. math::
    :label: ace-law-1-energy

    E' = E_{min} + \frac{\hat{E} - E_{\ell,1}}{E_{\ell,M} - E_{\ell,1}}
    (E_{max} - E_{min})

ACE Law 3 - Inelastic Level Scattering
++++++++++++++++++++++++++++++++++++++

It can be shown (see Foderaro_) that in inelastic level scattering, the outgoing
energy of the neutron :math:`E'` can be related to the Q-value of the reaction
and the incoming energy:

.. math::
    :label: level-scattering

    E' = \left ( \frac{A}{A+1} \right )^2 \left ( E - \frac{A + 1}{A} Q \right )

where :math:`A` is the mass of the target nucleus measured in neutron masses.

.. _ace-law-4:

ACE Law 4 - Continuous Tabular Distribution
+++++++++++++++++++++++++++++++++++++++++++

This representation is very similar to :ref:`ace-law-1` except that instead of
equiprobable outgoing energy bins, the outgoing energy distribution for each
incoming energy is represented with a probability distribution function. For
each incoming neutron energy :math:`E_i`, let us call :math:`p_{i,j}` the j-th
value in the probability distribution function, :math:`c_{i,j}` the j-th value
in the cumulative distribution function, and :math:`E_{i,j}` the j-th outgoing
energy.

We proceed first as we did for ACE Law 1, determining the bounding energies of
the particle's incoming energy such that :math:`E_i < E < E_{i+1}` and
calculating an interpolation factor :math:`f` with equation
:eq:`interpolation-factor`. Next, statistical interpolation is performed to
choose between using the outgoing energy distributions corresponding to energy
:math:`E_i` and :math:`E_{i+1}`. Let :math:`\ell` be the chosen table where
:math:`\ell = i` if :math:`\xi_1 > f` and :math:`\ell = i + 1` otherwise, and
:math:`\xi_1` is a random number. Then, we sample an outgoing energy bin
:math:`j` using the cumulative distribution function:

.. math::
    :label: ace-law-4-sample-cdf

    c_{\ell,j} < \xi_2 < c_{\ell,j+1}

where :math:`\xi_2` is a random number sampled uniformly on :math:`[0,1)`. At
this point, we need to interpolate between the successive values on the outgoing
energy distribution using either histogram or linear-linear interpolation. The
formulas for these can be derived along the same lines as those found in
:ref:`angle-tabular`. For histogram interpolation, the interpolated outgoing
energy on the :math:`\ell`-th distribution is

.. math::
    :label: energy-histogram

    \hat{E} = E_{\ell,j} + \frac{\xi_2 - c_{\ell,j}}{p_{\ell,j}}.

If linear-linear interpolation is to be used, the outgoing energy on the
:math:`\ell`-th distribution is

.. math::
    :label: energy-linlin

    \hat{E} = E_{\ell,j} + \frac{E_{\ell,j+1} - E_{\ell,j}}{p_{\ell,j+1} -
    p_{\ell,j}} \left ( \sqrt{p_{\ell,j}^2 + 2 \frac{p_{\ell,j+1} -
    p_{\ell,j}}{E_{\ell,j+1} - E_{\ell,j}} ( \xi_2 - c_{\ell,j} )} - p_{\ell,j}
    \right ).

Since this outgoing energy may violate reaction kinematics, we then scale it to
minimum and maximum energies interpolated between the neighboring outgoing
energy distributions to get the final outgoing energy:

.. math::
    :label: ace-law-4-energy

    E' = E_{min} + \frac{\hat{E} - E_{\ell,1}}{E_{\ell,M} - E_{\ell,1}}
    (E_{max} - E_{min})

where :math:`E_{min}` and :math:`E_{max}` are defined the same as in equation
:eq:`ace-law-1-minmax`.

.. _maxwell:

ACE Law 7 - Maxwell Fission Spectrum
++++++++++++++++++++++++++++++++++++

One representation of the secondary energies for neutrons from fission is the
so-called Maxwell spectrum. A probability distribution for the Maxwell spectrum
can be written in the form

.. math::
    :label: maxwell-spectrum

    p(E') dE' = c E'^{1/2} e^{-E'/T(E)} dE'

where :math:`E` is the incoming energy of the neutron and :math:`T` is the
so-called nuclear temperature, which is a function of the incoming energy of the
neutron. The ACE format contains a list of nuclear temperatures versus incoming
energies. The nuclear temperature is interpolated between neighboring incoming
energies using a specified interpolation law. Once the temperature :math:`T` is
determined, we then calculate a candidate outgoing energy based on rule C64 in
the `Monte Carlo Sampler`_:

.. math::
    :label: maxwell-E-candidate

    E' = -T \left [ \log (\xi_1) + \log (\xi_2) \cos^2 \left ( \frac{\pi
    \xi_3}{2} \right ) \right ]

where :math:`\xi_1, \xi_2, \xi_3` are random numbers sampled on the unit
interval. The outgoing energy is only accepted if

.. math::
    :label: maxwell-restriction

    0 \le E' \le E - U

where :math:`U` is called the restriction energy and is specified on the ACE
table. If the outgoing energy is rejected, it is resampled using equation
:eq:`maxwell-E-candidate`.

ACE Law 9 - Evaporation Spectrum
++++++++++++++++++++++++++++++++

Evaporation spectra are primarily used in compound nucleus processes where a
secondary particle can "evaporate" from the compound nucleus if it has
sufficient energy. The probability distribution for an evaporation spectrum can
be written in the form

.. math::
    :label: evaporation-spectrum

    p(E') dE' = c E' e^{-E'/T(E)} dE'

where :math:`E` is the incoming energy of the neutron and :math:`T` is the
nuclear temperature, which is a function of the incoming energy of the
neutron. The ACE format contains a list of nuclear temperatures versus incoming
energies. The nuclear temperature is interpolated between neighboring incoming
energies using a specified interpolation law. Once the temperature :math:`T` is
determined, we then calculate a candidate outgoing energy based on rule C45 in
the `Monte Carlo Sampler`_:

.. math::
    :label: evaporation-E

    E' = -T \log (\xi_1 \xi_2)

where :math:`\xi_1, \xi_2` are random numbers sampled on the unit
interval. The outgoing energy is only accepted according to a specified
restriction energy as in equation :eq:`maxwell-restriction`.

ACE Law 11 - Energy-Dependent Watt Spectrum
+++++++++++++++++++++++++++++++++++++++++++

The probability distribution for a Watt fission spectrum can be written in the
form

.. math::
    :label: watt-spectrum

    p(E') dE' = c e^{-E'/a(E)} \sinh \sqrt{b(E) \, E'} dE'

where :math:`a` and :math:`b` are parameters for the distribution and are given
as tabulated functions of the incoming energy of the neutron. These two
parameters are interpolated on the incoming energy grid using a specified
interpolation law. Once the parameters have been determined, we sample a
Maxwellian spectrum with nuclear temperature :math:`a` using the algorithm
described in :ref:`maxwell` to get an energy :math:`W`. Then, the outgoing
energy is calculated as

.. math::
    :label: watt-E

    E' = W + \frac{a^2 b}{4} + (2\xi - 1) \sqrt{a^2 b W}

where :math:`\xi` is a random number sampled on the interval :math:`[0,1)`. The
outgoing energy is only accepted according to a specified restriction energy
:math:`U` as defined in equation :eq:`maxwell-restriction`.

This algorithm can be found in Forrest Brown's lectures_ on Monte Carlo methods
and is an unpublished sampling scheme based on the original Watt spectrum
derivation [Watt]_.

ACE Law 44 - Kalbach-Mann Correlated Scattering
+++++++++++++++++++++++++++++++++++++++++++++++

This law is very similar to ACE Law 4 except now the outgoing angle of the
neutron is correlated to the outgoing energy and is not sampled from a separate
distribution. For each incident neutron energy :math:`E_i` tabulated, there is
an array of precompound factors :math:`R_{i,j}` and angular distribution slopes
:math:`A_{i,j}` corresponding to each outgoing energy bin :math:`j` in addition
to the outgoing energies and distribution functions as in ACE Law 4.

The calculation of the outgoing energy of the neutron proceeds exactly the same
as in the algorithm described in :ref:`ace-law-4`. In that algorithm, we found
an interpolation factor :math:`f`, statistically sampled an incoming energy bin
:math:`\ell`, and sampled an outgoing energy bin :math:`j` based on the
tabulated cumulative distribution function. Once the outgoing energy has been
determined with equation :eq:`ace-law-4-energy`, we then need to calculate the
outgoing angle based on the tabulated Kalbach-Mann parameters. These parameters
themselves are subject to either histogram or linear-linear interpolation on the
outgoing energy grid. For histogram interpolation, the parameters are

.. math::
    :label: KM-parameters-histogram

    R = R_{\ell,j} \\
    A = A_{\ell,j}.

If linear-linear interpolation is specified, the parameters are

.. math::
    :label: KM-parameters-linlin

    R = R_{\ell,j} + \frac{\hat{E} - E_{\ell,j}}{E_{\ell,j+1} - E_{\ell,j}} (
    R_{\ell,j+1} - R_{\ell,j} ) \\
    A = A_{\ell,j} + \frac{\hat{E} - E_{\ell,j}}{E_{\ell,j+1} - E_{\ell,j}} (
    A_{\ell,j+1} - A_{\ell,j} )

where :math:`\hat{E}` is defined in equation :eq:`energy-linlin`. With the
parameters determined, the probability distribution function for the cosine of
the scattering angle is

.. math::
    :label: KM-pdf-angle

    p(\mu) d\mu = \frac{A}{2 \sinh (A)} \left [ \cosh (A\mu) + R \sinh (A\mu)
    \right ] d\mu.

The rules for sampling this probability distribution function can be derived
based on rules C39 and C40 in the `Monte Carlo Sampler`_. First, we sample two
random numbers :math:`\xi_3, \xi_4` on the unit interval. If :math:`\xi_3 > R`
then the outgoing angle is

.. math::
    :label: KM-angle-1

    \mu = \frac{1}{A} \ln \left ( T + \sqrt{T^2 + 1} \right )

where :math:`T = (2 \xi_4 - 1) \sinh (A)`. If :math:`\xi_3 \le R`, then the
outgoing angle is

.. math::
    :label: KM-angle-2

    \mu = \frac{1}{A} \ln \left ( \xi_4 e^A + (1 - \xi_4) e^{-A} \right ).

.. _ace-law-61:

ACE Law 61 - Correlated Energy and Angle Distribution
+++++++++++++++++++++++++++++++++++++++++++++++++++++

This law is very similar to ACE Law 44 in the sense that the outgoing angle of
the neutron is correlated to the outgoing energy and is not sampled from a
separate distribution. In this case though, rather than being determined from an
analytical distribution function, the cosine of the scattering angle is
determined from a tabulated distribution. For each incident energy :math:`i` and
outgoing energy :math:`j`, there is a tabulated angular distribution.

The calculation of the outgoing energy of the neutron proceeds exactly the same
as in the algorithm described in :ref:`ace-law-4`. In that algorithm, we found
an interpolation factor :math:`f`, statistically sampled an incoming energy bin
:math:`\ell`, and sampled an outgoing energy bin :math:`j` based on the
tabulated cumulative distribution function. Once the outgoing energy has been
determined with equation :eq:`ace-law-4-energy`, we then need to decide which
angular distribution to use. If histogram interpolation was used on the outgoing
energy bins, then we use the angular distribution corresponding to incoming
energy bin :math:`\ell` and outgoing energy bin :math:`j`. If linear-linear
interpolation was used on the outgoing energy bins, then we use the whichever
angular distribution was closer to the sampled value of the cumulative
distribution function for the outgoing energy. The actual algorithm used to
sample the chosen tabular angular distribution has been previously described in
:ref:`angle-tabular`.

ACE Law 66 - N-Body Phase Space Distribution
++++++++++++++++++++++++++++++++++++++++++++

Reactions in which there are more than two products of similar masses are
sometimes best treated by using what's known as an N-body phase
distribution. This distribution has the following probability density function
for outgoing energy of the :math:`i`-th particle in the center-of-mass system:

.. math::
    :label: n-body-pdf

    p_i(E') dE' = C_n \sqrt{E'} (E_i^{max} - E')^{(3n/2) - 4} dE'

where :math:`n` is the number of outgoing particles, :math:`C_n` is a
normalization constant, :math:`E_i^{max}` is the maximum center-of-mass energy
for particle :math:`i`, and :math:`E'` is the outgoing energy. The algorithm for
sampling the outgoing energy is based on algorithms R28, C45, and C64 in the
`Monte Carlo Sampler`_. First we calculate the maximum energy in the
center-of-mass using the following equation:

.. math::
    :label: n-body-emax

    E_i^{max} = \frac{A_p - 1}{A_p} \left ( \frac{A}{A+1} E + Q \right )

where :math:`A_p` is the total mass of the outgoing particles in neutron masses,
:math:`A` is the mass of the original target nucleus in neutron masses, and
:math:`Q` is the Q-value of the reaction. Next we sample a value :math:`x` from
a Maxwell distribution with a nuclear temperature of one using the algorithm
outlined in :ref:`maxwell`. We then need to determine a value :math:`y` that
will depend on how many outgoing particles there are. For :math:`n = 3`, we
simply sample another Maxwell distribution with unity nuclear temperature. For
:math:`n = 4`, we use the equation

.. math::
    :label: n-body-y4

    y = -\ln ( \xi_1 \xi_2 \xi_3 )

where :math:`\xi_i` are random numbers sampled on the interval
:math:`[0,1)`. For :math:`n = 5`, we use the equation

.. math::
    :label: n-body-y5

    y = -\ln ( \xi_1 \xi_2 \xi_3 \xi_4 ) - \ln ( \xi_5 ) \cos^2 \left (
    \frac{\pi}{2} \xi_6 \right )

After :math:`x` and :math:`y` have been determined, the outgoing energy is then
calculated as

.. math::
    :label: n-body-energy

    E' = \frac{x}{x + y} E_i^{max}

There are two important notes to make regarding the N-body phase space
distribution. First, the documentation (and code) for MCNP5-1.60 has a mistake
in the algorithm for :math:`n = 4`. That being said, there are no existing
nuclear data evaluations which use an N-body phase space distribution with
:math:`n = 4`, so the error would not affect any calculations. In the
ENDF/B-VII.0 nuclear data evaluation, only one reaction uses an N-body phase
space distribution at all, the :math:`(n,2n)` reaction with H-2.

.. _transform-coordinates:

-------------------------------------
Transforming a Particle's Coordinates
-------------------------------------

Once the cosine of the scattering angle :math:`\mu` has been sampled either from
a angle distribution or a correlated angle-energy distribution, we are still
left with the task of transforming the particle's coordinates. If the outgoing
energy and scattering cosine were given in the center-of-mass system, then we
first need to transform these into the laboratory system. The relationship
between the outgoing energy in center-of-mass and laboratory is

.. math::
    :label: energy-com-to-lab

    E' = E'_{cm} + \frac{E + 2\mu_{cm} (A + 1) \sqrt{EE'_{cm}}}{(A+1)^2}.

where :math:`E'_{cm}` is the outgoing energy in the center-of-mass system,
:math:`\mu_{cm}` is the scattering cosine in the center-of-mass system,
:math:`E'` is the outgoing energy in the laboratory system, and :math:`E` is the
incident neutron energy. The relationship between the scattering cosine in
center-of-mass and laboratory is

.. math::
    :label: angle-com-to-lab

    \mu = \mu_{cm} \sqrt{\frac{E'_{cm}}{E'}} + \frac{1}{A + 1}
    \sqrt{\frac{E}{E'}}

where :math:`\mu` is the scattering cosine in the laboratory system. The
scattering cosine still only tells us the cosine of the angle between the
original direction of the particle and the new direction of the particle. If we
express the pre-collision direction of the particle as :math:`\mathbf{\Omega} =
(u,v,w)` and the post-collision direction of the particle as
:math:`\mathbf{\Omega}' = (u',v',w')`, it is possible to relate the pre- and
post-collision components. We first need to uniformly sample an azimuthal angle
:math:`\phi` in :math:`[0, 2\pi)`. After the azimuthal angle has been sampled,
the post-collision direction is calculated as

.. math::
    :label: post-collision-angle

    u' = \mu u + \frac{\sqrt{1 - \mu^2} ( uw \cos\phi - v \sin\phi )}{\sqrt{1 -
    w^2}} \\

    v' = \mu v + \frac{\sqrt{1 - \mu^2} ( vw \cos\phi + u \sin\phi )}{\sqrt{1 -
    w^2}} \\

    w' = \mu w - \sqrt{1 - \mu^2} \sqrt{1 - w^2} \cos\phi.

.. _freegas:

------------------------------------------
Effect of Thermal Motion on Cross Sections
------------------------------------------

When a neutron scatters off of a nucleus, it may often be assumed that the
target nucleus is at rest. However, the target nucleus will have motion
associated with its thermal vibration, even at absolute zero (This is due to the
zero-point energy arising from quantum mechanical considerations). Thus, the
velocity of the neutron relative to the target nucleus is in general not the
same as the velocity of the neutron entering the collision.

The effect of the thermal motion on the interaction probability can be written
as

.. math::
    :label: doppler-broaden

    v_n \bar{\sigma} (v_n, T) = \int d\mathbf{v}_T v_r \sigma(v_r)
    M (\mathbf{v}_T)

where :math:`v_n` is the magnitude of the velocity of the neutron,
:math:`\bar{\sigma}` is an effective cross section, :math:`T` is the temperature
of the target material, :math:`\mathbf{v}_T` is the velocity of the target
nucleus, :math:`v_r = || \mathbf{v}_n - \mathbf{v}_T ||` is the magnitude of the
relative velocity, :math:`\sigma` is the cross section at 0 K, and :math:`M
(\mathbf{v}_T)` is the probability distribution for the target nucleus velocity
at temperature :math:`T` (a Maxwellian). In a Monte Carlo code, one must account
for the effect of the thermal motion on both the integrated cross section as
well as secondary angle and energy distributions. For integrated cross sections,
it is possible to calculate thermally-averaged cross sections by applying a
kernel Doppler broadening algorithm to data at 0 K (or some temperature lower
than the desired temperature). The most ubiquitous algorithm for this purpose is
the [SIGMA1]_ method developed by Red Cullen and subsequently refined by
others. This method is used in the NJOY_ and PREPRO_ data processing codes.

The effect of thermal motion on secondary angle and energy distributions can be
accounted for on-the-fly in a Monte Carlo simulation. We must first qualify
where it is actually used however. All threshold reactions are treated as being
independent of temperature, and therefore they are not Doppler broadened in NJOY
and no special procedure is used to adjust the secondary angle and energy
distributions. The only non-threshold reactions with secondary neutrons are
elastic scattering and fission. For fission, it is assumed that the neutrons are
emitted isotropically (this is not strictly true, but is nevertheless a good
approximation). This leaves only elastic scattering that needs a special thermal
treatment for secondary distributions.

Fortunately, it is possible to directly sample the velocity of the target
nuclide and then use it directly in the kinematic calculations. However, this
calculation is a bit more nuanced than it might seem at first glance. One might
be tempted to simply sample a Maxwellian distribution for the velocity of the
target nuclide.  Careful inspection of equation :eq:`doppler-broaden` however
tells us that target velocities that produce relative velocities which
correspond to high cross sections will have a greater contribution to the
effective reaction rate. This is most important when the velocity of the
incoming neutron is close to a resonance. For example, if the neutron's velocity
corresponds to a trough in a resonance elastic scattering cross section, a very
small target velocity can cause the relative velocity to correspond to the peak
of the resonance, thus making a disproportionate contribution to the reaction
rate. The conclusion is that if we are to sample a target velocity in the Monte
Carlo code, it must be done in such a way that preserves the thermally-averaged
reaction rate as per equation :eq:`doppler-broaden`.

The method by which most Monte Carlo codes sample the target velocity for use in
elastic scattering kinematics is outlined in detail by [Gelbard]_. The
derivation here largely follows that of Gelbard. Let us first write the reaction
rate as a function of the velocity of the target nucleus:

.. math::
    :label: reaction-rate

    R(\mathbf{v}_T) = || \mathbf{v}_n - \mathbf{v}_T || \sigma ( ||
    \mathbf{v}_n - \mathbf{v}_T || ) M ( \mathbf{v}_T )

where :math:`R` is the reaction rate. Note that this is just the right-hand side
of equation :eq:`doppler-broaden`. Based on the discussion above, we want to
construct a probability distribution function for sampling the target velocity
to preserve the reaction rate -- this is different from the overall probability
distribution function for the target velocity, :math:`M ( \mathbf{v}_T )`. This
probability distribution function can be found by integrating equation
:eq:`reaction-rate` to obtain a normalization factor:

.. math::
    :label: target-pdf-1

    p( \mathbf{v}_T ) d\mathbf{v}_T = \frac{R(\mathbf{v}_T) d\mathbf{v}_T}{\int
    d\mathbf{v}_T \, R(\mathbf{v}_T)}

Let us call the normalization factor in the denominator of equation
:eq:`target-pdf-1` :math:`C`.

It is normally assumed that :math:`\sigma (v_r)` is constant over the range of
relative velocities of interest. This is a good assumption for almost all cases
since the elastic scattering cross section varies slowly with velocity for light
nuclei, and for heavy nuclei where large variations can occur due to resonance
scattering, the moderating effect is rather small. Nonetheless, this assumption
may cause incorrect answers in systems with low-lying resonances that can cause
a significant amount of up-scatter that would be ignored by this assumption
(e.g. U-238 in commercial light-water reactors). Nevertheless, with this
assumption, we write :math:`\sigma (v_r) = \sigma_s` which simplifies
:eq:`target-pdf-1` to

.. math::
    :label: target-pdf-2

    p( \mathbf{v}_T ) d\mathbf{v}_T = \frac{\sigma_s}{C} || \mathbf{v}_n -
    \mathbf{v}_T || M ( \mathbf{v}_T ) d\mathbf{v}_T

The Maxwellian distribution in velocity is

.. math::
    :label: maxwellian-velocity

    M (\mathbf{v}_T) = \left ( \frac{m}{2\pi kT} \right )^{3/2} \exp \left (
    \frac{-m || \mathbf{v}_T^2 ||}{2kT} \right )

where :math:`m` is the mass of the target nucleus and :math:`k` is Boltzmann's
constant. Notice here that the term in the exponential is dependent only on the
speed of the target, not on the actual direction. Thus, we can change the
Maxwellian into a distribution for speed rather than velocity. The differential
element of velocity is

.. math::
    :label: differential-velocity

    d\mathbf{v}_T = v_T^2 dv_T d\mu d\phi

Let us define the Maxwellian distribution in speed as

.. math::
    :label: maxwellian-speed

    M (v_T) dv_T = \int_{-1}^1 d\mu \int_{0}^{2\pi} d\phi \, dv_T \, v_T^2
    M(\mathbf{v}_T) = \sqrt{ \frac{2}{\pi} \left ( \frac{m}{kT} \right )^3}
    v_T^2 \exp \left ( \frac{-m v_T}{2kT} \right ) dv_T.

To simplify things a bit, we'll define a parameter

.. math::
    :label: maxwellian-beta

    \beta = \sqrt{\frac{m}{2kT}}.

Substituting equation :eq:`maxwellian-beta` into equation
:eq:`maxwellian-speed`, we obtain

.. math::
    :label: maxwellian-speed2

    M (v_T) dv_T = \frac{4}{\sqrt{\pi}} \beta^3 v_T^2 \exp \left ( -\beta^2
    v_T^2 \right ) dv_T.

Now, changing variables in equation :eq:`target-pdf-2` by using the result from
equation :eq:`maxwellian-speed`, our new probability distribution function is

.. math::
    :label: target-pdf-3

    p( v_T, \mu ) dv_T d\mu = \frac{4\sigma_s}{\sqrt{\pi}C'} || \mathbf{v}_n -
    \mathbf{v}_T || \beta^3 v_T^2 \exp \left ( -\beta^2 v_T^2 \right ) dv_T d\mu

Again, the Maxwellian distribution for the speed of the target nucleus has no
dependence on the angle between the neutron and target velocity vectors. Thus,
only the term :math:`|| \mathbf{v}_n - \mathbf{v}_T ||` imposes any constraint
on the allowed angle. Our last task is to take that term and write it in terms
of magnitudes of the velocity vectors and the angle rather than the vectors
themselves. We can establish this relation based on the law of cosines which
tells us that

.. math::
    :label: lawcosine

    2 v_n v_T \mu = v_n^2 + v_T^2 - v_r^2.

Thus, we can infer that

.. math::
    :label: change-terms

    || \mathbf{v}_n - \mathbf{v}_T || = || \mathbf{v}_r || = v_r = \sqrt{v_n^2 +
       v_T^2 - 2v_n v_T \mu}.

Inserting equation :eq:`change-terms` into :eq:`target-pdf-3`, we obtain

.. math::
    :label: target-pdf-4

    p( v_T, \mu ) dv_T d\mu = \frac{4\sigma_s}{\sqrt{\pi}C'} \sqrt{v_n^2 +
       v_T^2 - 2v_n v_T \mu} \beta^3 v_T^2 \exp \left ( -\beta^2 v_T^2 \right )
       dv_T d\mu

This expression is still quite formidable and does not lend itself to any
natural sampling scheme. We can divide this probability distribution into two
parts as such:

.. math::
    :label: divide-pdf

    p(v_T, \mu) &= f_1(v_T, \mu) f_2(v_T) \\

    f_1(v_T, \mu) &= \frac{4\sigma_s}{\sqrt{\pi} C'} \frac{ \sqrt{v_n^2 +
       v_T^2 - 2v_n v_T \mu}}{v_n + v_T} \\

    f_2(v_T) &= (v_n + v_T) \beta^3 v_T^2 \exp \left ( -\beta^2 v_T^2 \right ).

In general, any probability distribution function of the form :math:`p(x) =
f_1(x) f_2(x)` with :math:`f_1(x)` bounded can be sampled by sampling
:math:`x'` from the distribution

.. math::
    :label: freegas-f2

    q(x) dx = \frac{f_2(x) dx}{\int f_2(x) dx}

and accepting it with probability

.. math::
    :label: freegas-accept

    p_{accept} = \frac{f_1(x')}{\max f_1(x)}

The reason for dividing and multiplying the terms by :math:`v_n + v_T` is to
ensure that the first term is bounded. In general, :math:`|| \mathbf{v}_n -
\mathbf{v}_T ||` can take on arbitrarily large values, but if we divide it by
its maximum value :math:`v_n + v_T`, then it ensures that the function will be
bounded. We now must come up with a sampling scheme for equation
:eq:`freegas-f2`. To determine :math:`q(v_T)`, we need to integrate :math:`f_2`
in equation :eq:`divide-pdf`. Doing so we find that

.. math::
    :label: integrate-f2

    \int_0^{\infty} dv_T (v_n + v_T) \beta^3 v_T^2 \exp \left ( -\beta^2 v_T^2
    \right ) = \frac{1}{4\beta} \left ( \sqrt{\pi} \beta v_n + 2 \right ).

Thus, we need to sample the probability distribution function

.. math::
    :label: freegas-f2-2

    q(v_T) dv_T = \left ( \frac{4\beta^2 v_n v_T^2}{\sqrt{\pi} \beta v_n + 2} +
    \frac{4\beta^4 v_T^3}{\sqrt{\pi} \beta v_n + 2} \right ) exp \left (
    -\beta^2 v_T^2 \right ).

Now, let us do a change of variables with the following definitions

.. math::
    :label: beta-to-x

    x = \beta v_T \\
    y = \beta v_n.

Substituting equation :eq:`beta-to-x` into equation :eq:`freegas-f2-2` along
with :math:`dx = \beta dv_T` and doing some crafty rearranging of terms yields

.. math::
    :label: freegas-f2-3

    q(x) dx = \left [ \left ( \frac{\sqrt{\pi} y}{\sqrt{\pi} y + 2} \right )
    \frac{4}{\sqrt{\pi}} x^2 e^{-x^2} + \left ( \frac{2}{\sqrt{\pi} y + 2}
    \right ) 2x^3 e^{-x^2} \right ] dx.

It's important to make note of the following two facts. First, the terms outside
the parentheses are properly normalized probability distribution functions that
can be sampled directly. Secondly, the terms inside the parentheses are always
less than unity. Thus, the sampling scheme for :math:`q(x)` is as follows. We
sample a random number :math:`\xi_1` on the interval :math:`[0,1)` and if

.. math::
    :label: freegas-alpha

    \xi_1 < \frac{2}{\sqrt{\pi} y + 2}

then we sample the probability distribution :math:`2x^3 e^{-x^2}` for :math:`x`
using rule C49 in the `Monte Carlo Sampler`_ which we can then use to determine
the speed of the target nucleus :math:`v_T` from equation
:eq:`beta-to-x`. Otherwise, we sample the probability distribution
:math:`\frac{4}{\sqrt{\pi}} x^2 e^{-x^2}` for :math:`x` using rule C61 in the
`Monte Carlo Sampler`_.

With a target speed sampled, we must then decide whether to accept it based on
the probability in equation :eq:`freegas-accept`. The cosine can be sampled
isotropically as :math:`\mu = 2\xi_2 - 1` where :math:`\xi_2` is a random number
on the unit interval. Since the maximum value of :math:`f_1(v_T, \mu)` is
:math:`4\sigma_s / \sqrt{\pi} C'`, we then sample another random number
:math:`\xi_3` and accept the sampled target speed and cosine if

.. math::
    :label: freegas-accept-2

    \xi_3 < \frac{\sqrt{v_n^2 + v_T^2 - 2 v_n v_T \mu}}{v_n + v_T}.

If is not accepted, then we repeat the process and resample a target speed and
cosine until a combination is found that satisfies equation
:eq:`freegas-accept-2`.

.. _sab_tables:

------------
|sab| Tables
------------

For neutrons with thermal energies, generally less than 4 eV, the kinematics of
scattering can be affected by chemical binding and crystalline effects of the
target molecule. If these effects are not accounted for in a simulation, the
reported results may be highly inaccurate. There is no general analytic
treatment for the scattering kinematics at low energies, and thus when nuclear
data is processed for use in a Monte Carlo code, special tables are created that
give cross sections and secondary angle/energy distributions for thermal
scattering that account for thermal binding effects. These tables are mainly
used for moderating materials such as light or heavy water, graphite, hydrogen
in ZrH, beryllium, etc.

The theory behind |sab| is rooted in quantum mechanics and is quite
complex. Those interested in first principles derivations for formulae relating
to |sab| tables should be referred to the excellent books by [Williams]_ and
[Squires]_. For our purposes here, we will focus only on the use of already
processed data as it appears in the ACE format.

Each |sab| table can contain the following:

- Thermal inelastic scattering cross section;
- Thermal elastic scattering cross section;
- Correlated energy-angle distributions for thermal inelastic and elastic
  scattering.

Note that when we refer to "inelastic" and "elastic" scattering now, we are
actually using these terms with respect to the *scattering system*. Thermal
inelastic scattering means that the scattering system is left in an excited
state; no particular nucleus is left in an excited state as would be the case
for inelastic level scattering. In a crystalline material, the excitation of the
scattering could correspond to the production of phonons. In a molecule, it
could correspond to the excitation of rotational or vibrational modes.

Both thermal elastic and thermal inelastic scattering are generally divided into
incoherent and coherent parts. Coherent elastic scattering refers to scattering
in crystalline solids like graphite or beryllium. These cross sections are
characterized by the presence of *Bragg edges* that relate to the crystal
structure of the scattering material. Incoherent elastic scattering refers to
scattering in hydrogenous solids such as polyethylene. As it occurs in ACE data,
thermal inelastic scattering includes both coherent and incoherent effects and
is dominant for most other materials including hydrogen in water.

Calculating Integrated Cross Sections
-------------------------------------

The first aspect of using |sab| tables is calculating cross sections to replace
the data that would normally appear on the incident neutron data, which do not
account for thermal binding effects. For incoherent elastic and inelastic
scattering, the cross sections are stored as linearly interpolable functions on
a specified energy grid. For coherent elastic data, the cross section can be
expressed as

.. math::
    :label: coherent-elastic-xs

    \sigma(E) = \frac{\sigma_c}{E} \sum_{E_i < E} f_i e^{-4WE_i}

where :math:`\sigma_c` is the effective bound coherent scattering cross section,
:math:`W` is the effective Debye-Waller coefficient, :math:`E_i` are the
energies of the Bragg edges, and :math:`f_i` are related to crystallographic
structure factors. Since the functional form of the cross section is just 1/E
and the proportionality constant changes only at Bragg edges, the
proportionality constants are stored and then the cross section can be
calculated analytically based on equation :eq:`coherent-elastic-xs`.

Outgoing Angle for Coherent Elastic Scattering
----------------------------------------------

Another aspect of using |sab| tables is determining the outgoing energy and
angle of the neutron after scattering. For incoherent and coherent elastic
scattering, the energy of the neutron does not actually change, but the angle
does change. For coherent elastic scattering, the angle will depend on which
Bragg edge scattered the neutron. The probability that edge :math:`i` will
scatter then neutron is given by

.. math::
    :label: coherent-elastic-probability

    \frac{f_i e^{-4WE_i}}{\sum_j f_j e^{-4WE_j}}.

After a Bragg edge has been sampled, the cosine of the angle of scattering is
given analytically by

.. math::
    :label: coherent-elastic-angle

    \mu = 1 - \frac{E_i}{E}

where :math:`E_i` is the energy of the Bragg edge that scattered the neutron.

Outgoing Angle for Incoherent Elastic Scattering
------------------------------------------------

For incoherent elastic scattering, the probability distribution for the cosine
of the angle of scattering is represent as a series of equally-likely discrete
cosines :math:`\mu_{i,j}` for each incoming energy :math:`E_i` on the thermal
elastic energy grid. First the outgoing angle bin :math:`j` is sampled. Then, if
the incoming energy of the neutron satisfies :math:`E_i < E < E_{i+1}` the final
cosine is

.. math::
    :label: incoherent-elastic-angle

    \mu = \mu_{i,j} + f (\mu_{i+1,j} - \mu_{i,j})

where the interpolation factor is defined as

.. math::
    :label: sab-interpolation-factor

    f = \frac{E - E_i}{E_{i+1} - E_i}.

Outgoing Energy and Angle for Inelastic Scattering
--------------------------------------------------

Each |sab| table provides a correlated angle-energy secondary distribution for
neutron thermal inelastic scattering.  There are three representations used
in the ACE thermal scattering data: equiprobable discrete outgoing
energies, non-uniform yet still discrete outgoing energies, and continuous
outgoing energies with corresponding probability and cumulative distribution
functions provided in tabular format.  These three representations all
represent the angular distribution in a common format, using a series of
discrete equiprobable outgoing cosines.

Equi-Probable Outgoing Energies
+++++++++++++++++++++++++++++++

If the thermal data was processed with :math:`iwt = 1` in NJOY, then the
outgoing energy spectra is represented in the ACE data as a set of discrete and
equiprobable outgoing energies.  The procedure to determine the outgoing energy
and angle is as such. First, the interpolation factor is determined from
equation :eq:`sab-interpolation-factor`.  Then, an outgoing energy bin is
sampled from a uniform distribution and then interpolated between values
corresponding to neighboring incoming energies:

.. math::
    :label: inelastic-energy

    E = E_{i,j} + f (E_{i+1,j} - E_{i,j})

where :math:`E_{i,j}` is the j-th outgoing energy corresponding to the i-th
incoming energy. For each combination of incoming and outgoing energies, there
is a series equiprobable outgoing cosines. An outgoing cosine bin is sampled
uniformly and then the final cosine is interpolated on the incoming energy grid:

.. math::
    :label: inelastic-angle

    \mu = \mu_{i,j,k} + f (\mu_{i+1,j,k} - \mu_{i,j,k})

where :math:`\mu_{i,j,k}` is the k-th outgoing cosine corresponding to the j-th
outgoing energy and the i-th incoming energy.

Skewed Equi-Probable Outgoing Energies
++++++++++++++++++++++++++++++++++++++

If the thermal data was processed with :math:`iwt=0` in NJOY, then the
outgoing energy spectra is represented in the ACE data according to the
following: the first and last outgoing energies have a relative probability of
1, the second and second-to-last energies have a relative probability of 4, and
all other energies have a relative probability of 10.  The procedure to
determine the outgoing energy and angle is similar to the method discussed
above, except that the sampled probability distribution is now skewed
accordingly.

Continuous Outgoing Energies
++++++++++++++++++++++++++++

If the thermal data was processed with :math:`iwt=2` in NJOY, then the
outgoing energy spectra is represented by a continuous outgoing energy spectra
in tabular form with linear-linear interpolation.  The sampling of the outgoing
energy portion of this format is very similar to :ref:`ACE Law 61<ace-law-61>`,
but the sampling of the correlated angle is performed as it was in the other
two representations discussed in this sub-section.  In the Law 61 algorithm,
we found an interpolation factor :math:`f`, statistically sampled an incoming
energy bin :math:`\ell`, and sampled an outgoing energy bin :math:`j` based on
the tabulated cumulative distribution function. Once the outgoing energy has
been determined with equation :eq:`ace-law-4-energy`, we then need to decide
which angular distribution data to use.  Like the linear-linear interpolation
case in Law 61, the angular distribution closest to the sampled value of the
cumulative distribution function for the outgoing energy is utilized.  The
actual algorithm utilized to sample the outgoing angle is shown in equation
:eq:`inelastic-angle`.

.. _probability_tables:

----------------------------------------------
Unresolved Resonance Region Probability Tables
----------------------------------------------

In the unresolved resonance energy range, resonances may be so closely spaced
that it is not possible for experimental measurements to resolve all
resonances. To properly account for self-shielding in this energy range, OpenMC
uses the probability table method [Levitt]_. For most thermal reactors, the use
of probability tables will not significantly affect problem results. However,
for some fast reactors and other problems with an appreciable flux spectrum in
the unresolved resonance range, not using probability tables may lead to
incorrect results.

Probability tables in the ACE format are generated from the UNRESR module in
NJOY following the method of Levitt. A similar method employed for the RACER and
MC21_ Monte Carlo codes is described in a paper by `Sutton and Brown`_. For the
discussion here, we will focus only on use of the probability table table as it
appears in the ACE format.

Each probability table for a nuclide contains the following information at a
number of incoming energies within the unresolved resonance range:

- Cumulative probabilities for cross section bands;
- Total cross section (or factor) in each band;
- Elastic scattering cross section (or factor) in each band;
- Fission cross section (or factor) in each band;
- :math:`(n,\gamma)` cross section (or factor) in each band; and
- Neutron heating number (or factor) in each band.

It should be noted that unresolved resonance probability tables affect only
integrated cross sections and no extra data need be given for secondary
angle/energy distributions. Secondary distributions for elastic and inelastic
scattering would be specified whether or not probability tables were present.

The procedure for determining cross sections in the unresolved range using
probability tables is as follows. First, the bounding incoming energies are
determined, i.e. find :math:`i` such that :math:`E_i < E < E_{i+1}`. We then
sample a cross section band :math:`j` using the cumulative probabilities for
table :math:`i`. This allows us to then calculate the elastic, fission, and
capture cross sections from the probability tables interpolating between
neighboring incoming energies. If interpolation is specified, then
the cross sections are calculated as

.. math::
    :label: ptables-linlin

    \sigma = \sigma_{i,j} + f (\sigma_{i+1,j} - \sigma{i,j})

where :math:`\sigma_{i,j}` is the j-th band cross section corresponding to the
i-th incoming neutron energy and :math:`f` is the interpolation factor defined
in the same manner as :eq:`sab-interpolation-factor`. If logarithmic
interpolation is specified, the cross sections are calculated as

.. math::
    :label: ptables-loglog

    \sigma = \exp \left ( \log \sigma_{i,j} + f \log
    \frac{\sigma_{i+1,j}}{\sigma_{i,j}} \right )

where the interpolation factor is now defined as

.. math::
    :label: log-interpolation-factor

    f = \frac{\log \frac{E}{E_i}}{\log \frac{E_{i+1}}{E_i}}.

A flag is also present in the probability table that specifies whether an
inelastic cross section should be calculated. If so, this is done from a normal
reaction cross section (either MT=51 or a special MT). Finally, if the
cross sections defined are above are specified to be factors and not true
cross sections, they are multiplied by the underlying smooth cross section in
the unresolved range to get the actual cross sections. Lastly, the total cross
section is calculated as the sum of the elastic, fission, capture, and inelastic
cross sections.

-----------------------------
Variance Reduction Techniques
-----------------------------

Survival Biasing
----------------

In problems with highly absorbing materials, a large fraction of neutrons may be
killed through absorption reactions, thus leading to tallies with very few
scoring events. To remedy this situation, an algorithm known as *survival
biasing* or *implicit absorption* (or sometimes *implicit capture*, even though
this is a misnomer) is commonly used.

In survival biasing, absorption reactions are prohibited from occurring and
instead, at every collision, the weight of neutron is reduced by probability of
absorption occurring, i.e.

.. math::
    :label: survival-biasing-weight

    w' = w \left ( 1 - \frac{\sigma_a (E)}{\sigma_t (E)} \right )

where :math:`w'` is the weight of the neutron after adjustment and :math:`w` is
the weight of the neutron before adjustment. A few other things need to be
handled differently if survival biasing is turned on. Although fission reactions
never actually occur with survival biasing, we still need to create fission
sites to serve as source sites for the next generation in the method of
successive generations. The algorithm for sampling fission sites is the same as
that described in :ref:`fission`. The only difference is in equation
:eq:`fission-neutrons`. We now need to produce

.. math::
    :label: fission-neutrons-survival

    \nu = \frac{w}{k} \frac{\nu_t \sigma_f(E)}{\sigma_t (E)}

fission sites, where :math:`w` is the weight of the neutron before being
adjusted. One should note this is just the expected number of neutrons produced
*per collision* rather than the expected number of neutrons produced given that
fission has already occurred.

Additionally, since survival biasing can reduce the weight of the neutron to
very low values, it is always used in conjunction with a weight cutoff and
Russian rouletting. Two user adjustable parameters :math:`w_c` and :math:`w_s`
are given which are the weight below which neutrons should undergo Russian
roulette and the weight should they survive Russian roulette. The algorithm for
Russian rouletting is as follows. After a collision if :math:`w < w_c`, then the
neutron is killed with probability :math:`1 - w/w_s`. If it survives, the weight
is set equal to :math:`w_s`. One can confirm that the average weight following
Russian roulette is simply :math:`w`, so the game can be considered "fair". By
default, the cutoff weight in OpenMC is :math:`w_c = 0.25` and the survival
weight is :math:`w_s = 1.0`. These parameters vary from one Monte Carlo code to
another.

----------
References
----------

.. [Doyas] Richard J. Doyas and Sterrett T. Perkins, "Interpolation of Tabular
   Secondary Neutron and Photon Energy Distributions," *Nucl. Sci. Eng.*,
   **50**, 390-392 (1972).

.. [Gelbard] Ely M. Gelbard, "Epithermal Scattering in VIM," FRA-TM-123, Argonne
   National Laboratory (1979).

.. [Levitt] Leo B. Levitt, "The Probability Table Method for Treating Unresolved
   Neutron Resonances in Monte Carlo Calculations," *Nucl. Sci. Eng.*, **49**,
   pp. 450-457 (1972).

.. [SIGMA1] Dermett E. Cullen and Charles R. Weisbin, "Exact Doppler Broadening
   of Tabulated Cross Sections," *Nucl. Sci. Eng.*, **60**, pp. 199-229 (1976).

.. [Squires] G. L. Squires, *Introduction to the Theory of Thermal Neutron
   Scattering*, Cambridge University Press (1978).

.. [Watt] B. E. Watt, "Energy Spectrum of Neutrons from Thermal Fission of
   U235," *Phys. Rev.*, **87** (6), 1037-1041 (1952).

.. [Williams] M. M. R. Williams, *The Slowing Down and Thermalization of
   Neutrons*, North-Holland Publishing Co., Amsterdam (1966). **Note:** This
   book can be obtained for free from the OECD_.

.. |sab| replace:: S(:math:`\alpha,\beta,T`)

.. _Foderaro: http://hdl.handle.net/1721.1/1716

.. _OECD: http://www.oecd-nea.org/dbprog/MMRW-BOOKS.html

.. _NJOY: http://t2.lanl.gov/codes.shtml

.. _PREPRO: http://www-nds.iaea.org/ndspub/endf/prepro/

.. _ENDF-6 Format: http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

.. _Monte Carlo Sampler: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721_3rdmcsampler.pdf

.. _MC21: http://www.osti.gov/bridge/servlets/purl/903083-HT5p1o/903083.pdf

.. _Sutton and Brown: http://www.osti.gov/bridge/product.biblio.jsp?osti_id=307911

.. _lectures: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-05-4983.pdf

.. _MCNP Manual: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/MCNP5_Manual_Volume_I_LA-UR-03-1987.pdf
