.. _methods_physics:

=======
Physics
=======

-----------------------------------------
Secondary Angles and Energy Distributions
-----------------------------------------

For any reactions with secondary neutrons, it is necessary to sample secondary
angle and energy distributions. This includes elastic and inelastic scattering,
fission, and (n,xn) reactions. In some cases, the distributions may be specified
separately, and in other cases, they may be specified as a correlated
angle-energy distribution. In this section, we will outline the methods used to
sample secondary distributions as well as how they are used to modify the state
of a particle.

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

where :math:`\xi` is a random number sampled uniformly on :math:`[0,1)`.

Equiprobable Angle Bin Distribution
+++++++++++++++++++++++++++++++++++

For a 32 equiprobable bin distribution, the procedure to determine the
scattering cosine is as follows. First, we select a random number :math:`\xi` to
sample a cosine bin :math:`i` such that

.. math::
    :label: equiprobable-bin

    i = 1 + \lfloor 32\xi \rfloor

The same random number can then also be used to interpolate between neighboring
:math:`\mu` values to get the final scattering cosine:

.. math::
    :label: equiprobable-cosine

    \mu = \mu_i + (32\xi - i) (\mu_{i+1} - \mu_i)

.. _angle-tabular:

Tabular Angular Distribution
++++++++++++++++++++++++++++

As the MCNP Manual points out, using an equiprobable bin distribution works well
for high-probability regions of the scattering cosine probability, but for
low-probability regions it is not very accurate. Thus, a more typical treatment
is to represent the scattering cosine with a tabular distribution. In this case,
we have a table of cosines and their corresponding values for a probability
distributions function and cumulative distribution function. For each incoming
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
:math:`\ell` be the chosen table where :math:`\ell = i` if :math:`\xi > f` and
:math:`\ell = i + 1` otherwise where :math:`\xi` is a random number. A different
random number is used to sample a scattering cosine bin :math:`j` using the
cumulative distribution function:

.. math::
    :label: sample-cdf

    c_{\ell,j} < \xi < c_{\ell,j+1}

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

For histogram interpolation, we have that :math:`p(\mu') = p_{\ell,j}`. Thus,
after integration we have that

.. math::
    :label: cumulative-dist-histogram

    c(\mu) = c_{\ell,j} + (\mu - \mu_{\ell,j}) p_{\ell,j} = \xi

Solving for the scattering cosine, we obtain the final form for histogram
interpolation:

.. math::
    :label: cosine-histogram

    \mu = \mu_{\ell,j} + \frac{\xi - c_{\ell,j}}{p_{\ell,j}}

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
    p_{\ell,j} \right ] d\mu'

Let us now make a change of variables using

.. math::
    :label: introduce-eta

    \eta = \frac{p_{\ell,j+1} - p_{\ell,j}}{\mu_{\ell,j+1} - \mu_{\ell,j}}
    (\mu' - \mu_{\ell,j})

Equation :eq:`cdf-linlin` then becomes

.. math::
    :label: cdf-linlin-eta

    c(\mu) = c_{\ell,j} + \frac{1}{m} \int_{p_{\ell,j}}^{m(\mu - \mu_{\ell,j}) +
    p_{\ell,j}} \eta \, d\eta

where we have used

.. math::
    :label: slope

    m = \frac{p_{\ell,j+1} - p_{\ell,j}}{\mu_{\ell,j+1} - \mu_{\ell,j}}

Integrating equation :eq:`cdf-linlin-eta`, we have

.. math::
    :label: cdf-linlin-integrated

    c(\mu) = c_{\ell,j} + \frac{1}{2m} \left ( \left [ m (\mu - \mu_{\ell,j} ) +
    p_{\ell,j} \right ]^2 - p_{\ell,j}^2 \right ) = \xi

Solving for :math:`\mu`, we have the final form for the scattering cosine using
linear-linear interpolation:

.. math::
    :label: cosine-linlin

    \mu = \mu_{\ell,j} + \frac{1}{m} \left ( \sqrt{p_{\ell,j}^2 + 2 m (\xi -
    c_{\ell,j} )} - p_{\ell,j} \right )

.. _sample-energy:

Sampling Secondary Energy and Correlated Angle/Energy Distributions
-------------------------------------------------------------------

For a reaction with secondary neutrons, it is necessary to determine the
outgoing energy of the neutrons. For anything other than elastic scattering, the
outgoing energy must be determined based on tabulated or parameterized data. The
`ENDF-6 Format`_ specifies a variety of ways that the secondary energy
distribution can be represented. ENDF File 5 contains uncorrelated energy
distribution where ENDF File 6 contains correlated energy-angle
distributions. The ACE format specifies its own representations based loosely on
the formats given in ENDF-6. In this section, we will describe how the outgoing
energy of secondary particles is determined based on each ACE law.

One of the subtleties in the ACE format is the fact that a single reaction can
have multiple secondary energy distributions. This is mainly useful for
reactions with multiple neutrons in the exit channel such as (n,2n) or
(n,3n). In these types of reactions, each neutron is emitted corresponding to a
different excitation level of the compound nucleus, and thus in general the
neutrons will originate from different energy distributions. The first step in
sampling a secondary energy is to sample between multiple energy distributions
if more than one is present.

Once a secondary energy distribution has been sampled, the procedure for
determining the outgoing energy will depend on which ACE law has been specified
for the data.

.. _ace-law-1:

ACE Law 1 - Tabular Equiprobable Energy Bins
++++++++++++++++++++++++++++++++++++++++++++

In the tabular equiprobable bin representation, an array of equiprobably
outgoing energy bins is given for a number of incident energies. While the
representation itself is simple, the complexity lies in how one interpolates
between incident as well as outgoing energies on such a table. If one does
simple interpolation between tables for neighboring incident energies, it is
possible for the resulting energies to violate laws governing the kinematics,
i.e. the outgoing energy may be outside the range of available energy in the
reaction.

To avoid this situation, the accepted practice is to use a process known as
scaled interpolation [Doyas]_. First, we find the tabulated incident energies
which bound the actual incoming energy of the particle, i.e. find :math:`i` such
that :math:`E_i < E < E_{i+1}` and calculate the interpolation factor :math:`f`
via :eq:`interpolation-factor`. Then, we intepolate between the minimum and
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
:math:`\ell = i` if :math:`\xi_1 > f` and :math:`\ell = i + 1` otherwise where
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

It can be shown [Foderaro]_ that in inelastic level scattering, the outgoing
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

Weproceed first as we did for ACE Law 1, determining the bounding energies of
the particle's incoming energy such that :math:`E_i < E < E_{i+1}` and
calculating an interpolationg factor :math:`f` with equation
:eq:`interpolation-factor`. Next, statistical interpolation is performed to
choose between using the outgoing energy distributions corresponding to energy
:math:`E_i` and :math:`E_{i+1}`. Let :math:`\ell` be the chosen table where
:math:`\ell = i` if :math:`\xi_1 > f` and :math:`\ell = i + 1` otherwise where
:math:`\xi_1` is a random number. Then, we sample an outgoing energy bin
:math:`j` using the cumulative distribution function:

.. math::
    :label: ace-law-4-sample-cdf

    c_{\ell,j} < \xi_2 < c_{\ell,j+1}

where :math:`\xi_2` is a random number sampled uniformly on :math:`[0,1)`. At
this point, we need to inteporlate between the successive values on the outgoing
energy distribution using either histogram or linear-linear interpolation. The
formulas for these can be derived along the same lines as those found in
:ref:`angle-tabular`. For histogram interpolation, the interpolated outgoing
energy on the :math:`\ell`-th distribution is

.. math::
    :label: energy-histogram

    \hat{E} = E_{\ell,j} + \frac{\xi_2 - c_{\ell,j}}{p_{\ell,j}}

If linear-linear interpolation is to be used, the outgoing energy on the
:math:`\ell`-th distribution is

.. math::
    :label: energy-linlin

    \hat{E} = E_{\ell,j} + \frac{E_{\ell,j+1} - E_{\ell,j}}{p_{\ell,j+1} -
    p_{\ell,j}} \left ( \sqrt{p_{\ell,j}^2 + 2 \frac{p_{\ell,j+1} -
    p_{\ell,j}}{E_{\ell,j+1} - E_{\ell,j}} ( \xi_2 - c_{\ell,j} )} - p_{\ell,j}
    \right )

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
as tabulated functions of the incoming energy of the neutron in the ACE
format. These two parameters are interpolated on the incoming energy grid using
a specified interpolation law. Once the parameters have been determined, we
sample a Maxwellian spectrum with nuclear temperature :math:`a` using the
algorithm described in :ref:`maxwell` to get an energy :math:`W`. Then, the
outgoing energy is calculated as

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
an array of precompoung factors :math:`R_{i,j}` and angular distribution slopes
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
    A = A_{\ell,j}

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
    \right ] d\mu

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

    \mu = \frac{1}{A} \ln \left ( \xi_4 e^A + (1 - \xi_4) e^{-A} \right )


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
distribution. First, the documentation (and code) for MCNP5 has a mistake in the
algorithm for :math:`n = 4`. That being said, there are no existing nuclear data
evaluations which use an N-body phase space distribution with :math:`n = 4`, so
the error would not affect any calculations. In the ENDF/B-VII.0 nuclear data
evaluation, only one reaction uses an N-body phase space distribution at all,
the (n,2n) reaction with H-2.

.. _rotate-angle:

Transforming a Particle's Coordinates
-------------------------------------

Once the cosine of the scattering angle :math:`\mu` has been sampled either from
a angle distribution or a correlated angle-energy distribution, we are still
left with the task of transforming the particle's coordinates. The scattering
cosine that we sampled only tells us the cosine of the angle between the
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

    w' = \mu w - \sqrt{1 - \mu^2} \sqrt{1 - w^2} \cos\phi

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
target velocity :math:`v_t`. The velocity of the center-of-mass system is
calculated as

.. math::
    :label: velocity-com

    \mathbf{v}_{cm} = \frac{\mathbf{v}_n + A \mathbf{v}_t}{A + 1}

where :math:`\mathbf{v}_n` is the velocity of the neutron and :math:`A` is the
atomic mass of the target nucleus measured in neutron masses (commonly referred
to as the atomic weight ratio). With the velocity of the center-of-mass
calculated, we can then determine the neutron's velocity in the center-of-mass
system:

.. math::
    :label: velocity-neutron-com

    \mathbf{V}_n = \mathbf{v}_n - \mathbf{v}_{cm}

where we have used uppercase :math:`\mathbf{V}` to denote the center-of-mass
system. The direction of the neutron in the center-of-mass system is

.. math::
    :label: angle-neutron-com

    \mathbf{\Omega}_n = \frac{\mathbf{V}_n}{|| \mathbf{V}_n ||}

At low energies, elastic scattering will be isotropic in the center-of-mass
system, but for higher energies, there may be p-wave and higher order scattering
that leads to anisotropic scattering. Thus, in general, we need to sample a
cosine of the scattering angle which we will refer to as :math:`\mu`. For
elastic scattering, the secondary angle distribution is always given in the
center-of-mass system and is sampled according to the procedure outlined in
:ref:`sample-angle`. After the cosine of the angle of scattering has been
sampled, we need to determine the neutron's new direction
:math:`\mathbf{\Omega}'_n` in the center-of-mass system. This is done with the
procedure in :ref:`rotate-angle`. The new direction is multiplied by the speed
of the neutron in the center-of-mass system to obtain the new velocity vector in
the center-of-mass:

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

For tally purposes, it is also important to keep track of the scattering cosine
in the lab system. If we know the scattering cosine in the center-of-mass, the
scattering cosine in the lab system can be calculated as

.. math::
    :label: cosine-lab

    \mu_{lab} = \frac{1 + A\mu}{\sqrt{A^2 + 2A\mu + 1}}.

However, this formula is only valid if the target was at rest. When the target
nucleus does have thermal motion, the cosine of the scattering angle can be
determined by simply taking the dot product of the neutron's initial and final
direction in the lab system.

.. _freegas:

------------------------------------------
Effect of Thermal Motion on Cross-Sections
------------------------------------------

When a neutron scatters off of a nucleus, many times it is assumed that the
target nucleus is at rest. However, if the material is at a temperature greater
than 0 K, it will have motion associated with the thermal vibration. Thus, the
velocity of the neutron relative to the target nucleus is in general not the
same as the velocity of the neutron entering the collision.

The affect of the thermal motion on the interaction probability can be written
as

.. math::
    :label: freegas1

    v_n \sigma (v_n, T) = \int_0^\infty d\mathbf{v}_T \sigma(v_r, 0)
    \mathbf{v}_r p(\mathbf{v}_T)
    
where :math:`v_n` is the magnitude of the velocity of the neutron,
:math:`\mathbf{v}_T` is the velocity of the target nucleus, :math:`\mathbf{v}_r`
is the relative velocity, and :math:`T` is the temperature of the target
material. In a Monte Carlo code, one must account for the effect of the thermal
motion on both the integrated cross-section as well as secondary angle and
energy distributions. For integrated cross-sections, it is possible to calculate
thermally-averaged cross-sections by applying a kernel Doppler broadening
algorithm to data at 0 K (or some temperature lower than the desired
temperature). The most ubiquitous algorithm for this purpose is the [SIGMA1]_
method developed by Red Cullen and subsequently refined by others. This method
is used in the NJOY_ and PREPRO_ data processing codes.

The effect of thermal motion on secondary angle and energy distributions can be
accounted for on-the-fly in a Monte Carlo simulation. We must first qualify
where it is actually used however. All threshold reactions are treated as being
independent of temperature, and therefore they are not Doppler broadened in NJOY
and no special procedure is used to adjust the secondary angle and energy
distributions. The only non-threshold reactions with secondary neutrons are
elastic scattering and fission. For fission, it is assumed that neutrons are
emitted isotropically (this is not strictly true, but is nevertheless a good
approximation). This leaves only elastic scattering that needs a special thermal
treatment for secondary distributions.

Fortunately, it is possible to directly sample the velocity of the target
nuclide and then use it directly in the kinematic calculations. However, this
calculation is a bit more nuanced than it might seem at first glance. One might
be tempted to simply sample a Maxwellian distribution for the velocity of the
target nuclide.  Careful inspection of equation :eq:`freegas1` however tells us
that target velocities that produce relative velocities which correspond to high
cross sections will have a greater contribution to the effective reaction
rate. This is most important when the velocity of the incoming neutron is close
to a resonance. For example, if the neutron's velocity corresponds to a trough
in a resonance elastic scattering cross-section, a very small target velocity
can cause the relative velocity to correspond to the peak of the resonance, thus
making a disproportionate contribution to the reaction rate. The conclusion is
that if we are to sample a target velocity in the Monte Carlo code, it must be
done in such a way that preserves the thermally-averaged reaction rate as per
equation :eq:`freegas`.

The method by which most Monte Carlo codes sample the target velocity for use in
elastic scattering kinematics is outlined in detail by [Gelbard]_. The
derivation here largely follows that of Gelbard. The first assumption we can
make is that the velocity distribution for the thermal motion is isotropic, i.e.

.. math::
    :label: freegas2

    p(\mathbf{v}_T) d\mathbf{v}_T = \frac{1}{4\pi} p(v_T) dv_T d\mu d\phi

With this assumption, we can now rewrite equation :eq:`freegas1` as

.. math::
    :label: freegas3

    v_n \sigma (v_n, T) = \frac{1}{2} \int_{-1}^1 d\mu \int\limits_{v_r > 0}
    dv_T v_r \sigma (v_r, 0) p(v_T)

after integrating over :math:`d\phi`. To change the outer variable of
integration from :math:`\mu` to :math:`v_r`, we can establish a relation between
these variables based on the law of cosines.

.. math::
    :label: lawcosine

    2 v_n v_T \mu = v_n^2 + v_T^2 - v_r^2

The probability distribution for the magnitude of the velocity of the target
nucleus and the angle between the neutron and target velocity is

.. math::
    :label: freegas4

    P(v_T, \mu) = \frac{\sigma (v_r, 0) v_r P(v_T)}{2 \sigma (v_n, T) v_n}

It is normally assumed that :math:`\sigma (v_r, 0)` is constant over the range
of relative velocities of interest. This is a good assumption for almost all
cases since the elastic scattering cross section varies slowly with velocity for
light nuclei, and for heavy nuclei where large variations can occur due to
resonance scattering, the moderating effect is rather small. Nonetheless, this
assumption can cause incorrect answers in systems with U-238 where the low-lying
resonances can cause a significant amount of up-scatter that would be ignored by
this assumption.

With this (sometimes incorrect) assumption, we see that the probability
distribution is proportional to

.. math::
    :label: freegas5

    P(v_T, \mu) \propto v_r P(v_T) = | v_n - v_T | P(v_T)

We can divide this probability distribution into two parts as such:

.. math::
    :label: freegas6

    P(v_T, \mu) &= f_1(v_T, \mu) f_2(v_T) \\
    f_1(v_T, \mu) &= \frac{| v_n - v_T |}{C (v_n + v_T)} \\
    f_2(v_T) &= (v_n + v_T) P(v_T)

where :math:`C = \int dv_T \sigma v_r P(v_T)`. In general, any probability
distribution function of the form :math:`p(x) = f_1(x) f_2(x)` with
:math:`f_1(x)` bounded can be sampled by sampling :math:`x_s` from the
distribution

.. math::
    :label: freegas7

    \frac{f_2(x)}{\int f_2(x) dx}

and accepting it with probability

.. math::
    :label: freegas8

    \frac{f_1(x_s)}{\max f_1(x)}

It is normally assumed that the velocity distribution of the target nucleus
assumes a Maxwellian distribution in velocity.

------------
|sab| Tables
------------

For neutrons with thermal energies, generally less than 4 eV, the kinematics of
scattering can be affected by chemical binding and crystalline effects of the
target molecule. If these effects are not accounted for in a simulation, the
reported results may be highly inaccurate. There is no general analytic
treatment for the scattering kinematics at low energies, and thus when nuclear
data is processed for use in a Monte Carlo code, special tables are created that
give altered cross-sections and secondary angle/energy distributions for thermal
scattering. These tables are mainly used for moderating materials such as light
or heavy water, graphite, hydrogen in ZrH, beryllium, etc.

The theory behind |sab| is rooted in quantum mechanics and is quite
complex. Those interested in first principles derivations for formulae relating
to |sab| tables should be referred to the excellent books by [Williams]_ and
[Squires]_. For our purposes here, we will focus only on the use of already
processed data as it appears in the ACE format.

Each |sab| table can contain the following:

- Thermal inelastic scattering cross section
- Thermal elastic scattering cross section
- Correlated energy-angle distributions for thermal inelastic and elastic
  scattering

Note that when we refer to "inelastic" and "elastic" scattering now, we are
actually using these terms with respect to the *scattering system*. Thermal
inelastic scattering means that the scattering system is left in an excited
state, not any particular nucleus as is the case in inelastic level
scattering. In a crystalline material, the excitation could be the production of
phonons. In a molecule, it could be the excitation of rotational or vibrational
modes.

Both thermal elastic and thermal inelastic scattering are generally divided into
incoherent and coherent parts. Coherent elastic scattering refers to scattering
in crystalline solids like graphite or beryllium. These cross-sections are
characterized by the presence of "Bragg edges" that relate to the crystal
structure of the scattering material. Incoherent elastic scattering refers to
scattering in hydrogenous solids such as polyethylene. As it occurs in ACE data,
thermal inelastic scattering includes both coherent and incoherent effects and
is dominant for most other materials including hydrogen in water.

Calculating Integrated Cross Sections
-------------------------------------

The first aspect of using |sab| tables is calculating cross-sections to replace
the data that would normally appear on the incident neutron data, which do not
account for thermal binding effects. For incoherent elastic and inelastic
scattering, the cross-sections are stored as linearly interpolable functions on
a specified energy grid. For coherent elastic data, the cross section can be
expressed as

.. math::
    :label: coherent-elastic-xs

    \sigma(E) = \frac{\sigma_c}{E} \sum_{E_i < E} f_i e^{-4WE_i}.

where :math:`\sigma_c` is the effective bound coherent scattering cross section,
:math:`W` is the effective Debye-Waller coefficient, :math:`E_i` are the
energies of the Bragg edges, and :math:`f_i` are related to crystallographic
structure factors. Since the functional form of the cross-section is just 1/E
and the proportionality constant changes only at Bragg edges, the
proportionality constants are stored and then the cross-section can be
calculated analytically based on equation :eq:`coherent-elastic-xs`.

Outgoing Angle for Coherent Elastic Scattering
----------------------------------------------

The other aspect of using |sab| tables is determining the outgoing energy and
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

On each |sab| table, there is a correlated angle-energy secondary distribution
for neutron thermal inelastic scattering. While the documentation for the ACE
format implies that there are a series of equiprobable outgoing energies, the
outgoing energies may have non-uniform probability distribution. In particular,
if the thermal data were processed with :math:`iwt = 0` in NJOY, then the first
and last outgoing energies have a relative probability of 1, the second and
second to last energies have a relative probability of 4, and all other energies
have a relative probability of 10. The procedure to determine the outgoing
energy and angle is as such. First, the interpolation factor is determined from
equation :eq:`sab-interpolation-factor`. Then, an outgoing energy bin is sampled
either from a uniform distribution or from a skewed distribution as
discussed. The outgoing energy is then interpolated between values corresponding
to neighboring incoming energies:

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

- Cumulative probabilities for cross section bands
- Total cross section (or factor) in each band
- Elastic scattering cross section (or factor) in each band
- Fission cross section (or factor) in each band
- :math:`(n,\gamma)` cross section (or factor) in each band
- Neutron heating number (or factor) in each band

It should be noted that unresolved resonance probability tables affect only
integrated cross sections and no extra data need be given for secondary
angle/energy distributions. Secondary distributions for elastic and inelastic
scattering would be specified whether or not probability tables were present.

The procedure for determining cross sections in the unresolved range using
probability tables is as follows. First, the bounding incoming energies are
determined, i.e. find :math:`i` such that :math:`E_i < E < E_{i+1}`. We then
sample a cross section band :math:`j` using the cumulative probabilities for
table :math:`i`. This allows us to then calculate the elastic, fission, and
capture cross-sections from the probability tables interpolating between
neighboring incoming energies. If interpolation is specified, then
the cross sections are calculated as

.. math::
    :label: ptables-linlin

    \sigma = \sigma_{i,j} + f (\sigma_{i+1,j} - \sigma{i,j})

where :math:`f` is the interpolation factor defined in the same manner as
:eq:`sab-interpolation-factor`. If logarithmic interpolation is specified, the
cross sections are calculated as

.. math::
    :label: ptables-loglog

    \sigma = \exp \left ( \log \sigma_{i,j} + f \log
    \frac{\sigma_{i+1,j}}{\sigma_{i,j}} \right )

where the interpolation factor is now defined as

.. math::
    :label: log-interpolation-factor

    f = \frac{\log \frac{E}{E_i}}{\log \frac{E_{i+1}}{E_i}}

A flag is also present in the probability table that specifies whether an
inelastic cross section should be calculated. If so, this is done from a normal
reaction cross section (either MT=51 or a special MT). Finally, if the
cross-sections defined are above are specified to be factors and not true
cross-sections, they are multiplied by the underlying smooth cross section in
the unresolved range to get the actual cross sections. Lastly, the total cross
section is calculated as the sum of the elastic, fission, capture, and inelastic
cross sections.

----------
References
----------

.. [Doyas] Richard J. Doyas and Sterrett T. Perkins, "Interpolation of Tabular
   Secondary Neutron and Photon Energy Distributions," *Nucl. Sci. Eng.*,
   **50**, 390-392 (1972).

.. [Foderaro] Anthony Foderaro, *The Elements of Neutron Interaction Theory*,
   MIT Press, Cambridge, Massachusetts (1971). **Note:** Students, faculty, and
   staff at MIT can obtian a PDF copy of this book for free from the `MIT
   Press`_.

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

.. |sab| replace:: S(:math:`\alpha,\beta`)

.. _OECD: http://www.oecd-nea.org/dbprog/MMRW-BOOKS.html

.. _NJOY: http://t2.lanl.gov/codes.shtml

.. _PREPRO: http://www-nds.iaea.org/ndspub/endf/prepro/

.. _ENDF-6 Format: http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

.. _Monte Carlo Sampler: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721_3rdmcsampler.pdf

.. _MC21: http://www.osti.gov/bridge/servlets/purl/903083-HT5p1o/903083.pdf

.. _Sutton and Brown: http://www.osti.gov/bridge/product.biblio.jsp?osti_id=307911

.. _MIT Press: http://hdl.handle.net/1721.1/1716

.. _lectures: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-05-4983.pdf
