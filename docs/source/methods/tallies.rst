.. _methods_tallies:

=======
Tallies
=======

Note that the methods discussed in this section are written specifically for
continuous-energy mode but equivalent apply to the multi-group mode if the
particle's energy is replaced with the particle's group

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
from scattering, higher scattering moments, :math:`(n,xn)` reaction rates,
absorption reaction rate, fission reaction rate, neutron production rate from
fission, and surface currents. The following variables can be used as filters:
universe, material, cell, birth cell, surface, mesh, pre-collision energy,
post-collision energy, polar angle, azimuthal angle, and the cosine of the
change-in-angle due to a scattering event.

With filters for pre- and post-collision energy and scoring functions for
scattering and fission production, it is possible to use OpenMC to generate
cross sections with user-defined group structures. These multigroup cross
sections can subsequently be used in deterministic solvers such as coarse mesh
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
only applies to discrete filter variables and cannot be applied to energy,
angle, or change-in-angle bins. For these filters, it is necessary to perform
a binary search on the specified energy grid.

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
by the volume of the region of integration. If survival biasing is employed, the
analog estimator cannot be used for any reactions with zero neutrons in the exit
channel.

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
    \psi(\mathbf{r}, \mathbf{\hat{\Omega}}, E, t)

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
    \mathbf{\hat{\Omega}}, E, t)).

Using the relations :math:`N(\mathbf{r}, E, t) = \int d\mathbf{\Omega}
n(\mathbf{r}, \mathbf{\hat{\Omega}}, E, t)` and :math:`d\ell = v \, dt` where
:math:`d\ell` is the differential unit of track length, we then obtain

.. math::
    :label: track-length-integral

    V \phi = \int d\mathbf{r} \int dE \int d\ell N(\mathbf{r}, E, t).

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

    R_x = \frac{1}{W} \sum_{i \in T} w_i \ell_i \Sigma_x (E_i).

One important fact to take into consideration is that the use of a track-length
estimator precludes us from using any filter that requires knowledge of the
particle's state following a collision because by definition, it will not have
had a collision at every event. Thus, for tallies with outgoing-energy filters
(which require the post-collision energy), scattering change-in-angle filters,
or for tallies of scattering moments (which require the scattering cosine of
the change-in-angle), we must use an analog estimator.

.. TODO: Add description of surface current tallies

----------
Statistics
----------

As was discussed briefly in :ref:`methods_introduction`, any given result from a
Monte Carlo calculation, colloquially known as a "tally", represents an estimate
of the mean of some `random variable`_ of interest. This random variable
typically corresponds to some physical quantity like a reaction rate, a net
current across some surface, or the neutron flux in a region. Given that all
tallies are produced by a `stochastic process`_, there is an associated
uncertainty with each value reported. It is important to understand how the
uncertainty is calculated and what it tells us about our results. To that end,
we will introduce a number of theorems and results from statistics that should
shed some light on the interpretation of uncertainties.

Law of Large Numbers
--------------------

The `law of large numbers`_ is an important statistical result that tells us
that the average value of the result a large number of repeated experiments
should be close to the `expected value`_. Let :math:`X_1, X_2, \dots, X_n` be an
infinite sequence of `independent, identically-distributed random variables`_
with expected values :math:`E(X_1) = E(X_2) = \mu`. One form of the law of large
numbers states that the sample mean :math:`\bar{X_n} = \frac{X_1 + \dots +
X_n}{n}` `converges in probability`_ to the true mean, i.e. for all
:math:`\epsilon > 0`

.. math::

    \lim\limits_{n\rightarrow\infty} P \left ( \left | \bar{X}_n - \mu \right |
    \ge \epsilon \right ) = 0.

.. _central-limit-theorem:

Central Limit Theorem
---------------------

The `central limit theorem`_ (CLT) is perhaps the most well-known and ubiquitous
statistical theorem that has far-reaching implications across many
disciplines. The CLT is similar to the law of large numbers in that it tells us
the limiting behavior of the sample mean. Whereas the law of large numbers tells
us only that the value of the sample mean will converge to the expected value of
the distribution, the CLT says that the distribution of the sample mean will
converge to a `normal distribution`_. As we defined before, let :math:`X_1, X_2,
\dots, X_n` be an infinite sequence of independent, identically-distributed
random variables with expected values :math:`E(X_i) = \mu` and variances
:math:`\text{Var} (X_i) = \sigma^2 < \infty`. Note that we don't require that
these random variables take on any particular distribution -- they can be
normal, log-normal, Weibull, etc. The central limit theorem states that as
:math:`n \rightarrow \infty`, the random variable :math:`\sqrt{n} (\bar{X}_n -
\mu)` `converges in distribution`_ to the standard normal distribution:

.. math::
    :label: central-limit-theorem

    \sqrt{n} \left ( \frac{1}{n} \sum_{i=1}^n X_i - \mu \right ) \xrightarrow{d}
    \mathcal{N} (0, \sigma^2)

Estimating Statistics of a Random Variable
------------------------------------------

Mean
++++

Given independent samples drawn from a random variable, the sample mean is
simply an estimate of the average value of the random variable. In a Monte Carlo
simulation, the random variable represents physical quantities that we want
tallied. If :math:`X` is the random variable with :math:`N` observations
:math:`x_1, x_2, \dots, x_N`, then an unbiased estimator for the population mean
is the sample mean, defined as

.. math::
    :label: sample-mean

    \bar{x} = \frac{1}{N} \sum_{i=1}^N x_i.

Variance
++++++++

The variance of a population indicates how spread out different members of the
population are. For a Monte Carlo simulation, the variance of a tally is a
measure of how precisely we know the tally value, with a lower variance
indicating a higher precision. There are a few different estimators for the
population variance. One of these is the second central moment of the
distribution also known as the biased sample variance:

.. math::
    :label: biased-variance

    s_N^2 = \frac{1}{N} \sum_{i=1}^N \left ( x_i - \bar{x} \right )^2 = \left (
    \frac{1}{N} \sum_{i=1}^N x_i^2 \right ) - \bar{x}^2.

This estimator is biased because its expected value is actually not equal to the
population variance:

.. math::
    :label: biased-variance-expectation

    E[s_N^2] = \frac{N - 1}{N} \sigma^2

where :math:`\sigma^2` is the actual population variance. As a result, this
estimator should not be used in practice. Instead, one can use `Bessel's
correction`_ to come up with an unbiased sample variance estimator:

.. math::
    :label: unbiased-variance

    s^2 = \frac{1}{N - 1} \sum_{i=1}^N \left ( x_i - \bar{x} \right )^2 =
    \frac{1}{N - 1} \left ( \sum_{i=1}^N x_i^2 - N\bar{x}^2 \right ).

This is the estimator normally used to calculate sample variance. The final form
in equation :eq:`unbiased-variance` is especially suitable for computation since
we do not need to store the values at every realization of the random variable
as the simulation proceeds. Instead, we can simply keep a running sum and sum of
squares of the values at each realization of the random variable and use that to
calculate the variance.

Variance of the Mean
++++++++++++++++++++

The previous sections discussed how to estimate the mean and variance of a
random variable using statistics on a finite sample. However, we are generally
not interested in the *variance of the random variable* itself; we are more
interested in the *variance of the estimated mean*. The sample mean is the
result of our simulation, and the variance of the sample mean will tell us how
confident we should be in our answers.

Fortunately, it is quite easy to estimate the variance of the mean if we are
able to estimate the variance of the random variable. We start with the
observation that if we have a series of uncorrelated random variables, we can
write the variance of their sum as the sum of their variances:

.. math::
    :label: bienayme-formula

    \text{Var} \left ( \sum_{i=1}^N X_i \right ) = \sum_{i=1}^N \text{Var} \left
    ( X_i \right )

This result is known as the Bienaymé formula. We can use this result to
determine a formula for the variance of the sample mean. Assuming that the
realizations of our random variable are again identical,
independently-distributed samples, then we have that

.. math::
    :label: sample-variance-mean

    \text{Var} \left ( \bar{X} \right ) = \text{Var} \left ( \frac{1}{N}
    \sum_{i=1}^N X_i \right ) = \frac{1}{N^2} \sum_{i=1}^N \text{Var} \left (
    X_i \right ) = \frac{1}{N^2} \left ( N\sigma^2 \right ) =
    \frac{\sigma^2}{N}.

We can combine this result with equation :eq:`unbiased-variance` to come up with
an unbiased estimator for the variance of the sample mean:

.. math::
    :label: sample-variance-mean-formula

    s_{\bar{X}}^2 = \frac{1}{N - 1} \left ( \frac{1}{N} \sum_{i=1}^N x_i^2 -
    \bar{x}^2 \right ).

At this point, an important distinction should be made between the estimator for
the variance of the population and the estimator for the variance of the
mean. As the number of realizations increases, the estimated variance of the
population based on equation :eq:`unbiased-variance` will tend to the true
population variance. On the other hand, the estimated variance of the mean will
tend to zero as the number of realizations increases. A practical interpretation
of this is that the longer you run a simulation, the better you know your
results. Therefore, by running a simulation long enough, it is possible to
reduce the stochastic uncertainty to arbitrarily low levels.

Confidence Intervals
++++++++++++++++++++

While the sample variance and standard deviation gives us some idea about the
variability of the estimate of the mean of whatever quantities we've tallied, it
does not help us interpret how confidence we should be in the results. To
quantity the reliability of our estimates, we can use `confidence intervals`_
based on the calculated sample variance.

A :math:`1-\alpha` confidence interval for a population parameter is defined as
such: if we repeat the same experiment many times and calculate the confidence
interval for each experiment, then :math:`1 - \alpha` percent of the calculated
intervals would encompass the true population parameter. Let :math:`x_1, x_2,
\dots, x_N` be samples from a set of independent, identically-distributed random
variables each with population mean :math:`\mu` and variance
:math:`\sigma^2`. The t-statistic is defined as

.. math::
    :label: t-statistic

    t = \frac{\bar{x} - \mu}{s/\sqrt{N}}

where :math:`\bar{x}` is the sample mean from equation :eq:`sample-mean` and
:math:`s` is the standard deviation based on equation
:eq:`unbiased-variance`. If the random variables :math:`X_i` are
normally-distributed, then the t-statistic has a `Student's t-distribution`_
with :math:`N-1` degrees of freedom. This implies that

.. math::
    :label: t-probability

    Pr \left ( -t_{1 - \alpha/2, N - 1} \le \frac{\bar{x} - \mu}{s/\sqrt{N}} \le
    t_{1 - \alpha/2, N - 1} \right ) = 1 - \alpha

where :math:`t_{1-\alpha/2, N-1}` is the :math:`1 - \alpha/2` percentile of a
t-distribution with :math:`N-1` degrees of freedom. Thus, the :math:`1 - \alpha`
two sided confidence interval for the sample mean is

.. math::
    :label: two-sided-ci

    \bar{x} \pm t_{1 - \alpha/2, N-1} \frac{s}{\sqrt{N}}.

One should be cautioned that equation :eq:`two-sided-ci` only applies if the
*underlying random variables* are normally-distributed. In general, this may not
be true for a tally random variable --- the central limit theorem guarantees
only that the sample mean is normally distributed, not the underlying random
variable. If batching is used, then the underlying random variable, which would
then be the averages from each batch, will be normally distributed as long as
the conditions of the central limit theorem are met.

Let us now outline the method used to calculate the percentile of the Student's
t-distribution. For one or two degrees of freedom, the percentile can be written
analytically. For one degree of freedom, the t-distribution becomes a standard
`Cauchy distribution`_ whose cumulative distribution function is

.. math::
    :label: cauchy-cdf

    c(x) = \frac{1}{\pi} \arctan x + \frac{1}{2}.

Thus, inverting the cumulative distribution function, we find the :math:`x`
percentile of the standard Cauchy distribution to be

.. math::
    :label: percentile-1

    t_{x,1} = \tan \left ( \pi \left ( x - \frac{1}{2} \right ) \right ).

For two degrees of freedom, the cumulative distribution function is the
second-degree polynomial

.. math::
    :label: t-2-polynomial

    c(x) = \frac{1}{2} + \frac{x}{2\sqrt{x^2 + 2}}

Solving for :math:`x`, we find the :math:`x` percentile to be

.. math::
    :label: percentile-2

    t_{x,2} = \frac{2\sqrt{2} (x - 1/2)}{\sqrt{1 - 4 (x - 1/2)^2}}

For degrees of freedom greater than two, it is not possible to obtain an
analytical formula for the inverse of the cumulative distribution function. We
must resort to either numerically solving for the inverse or to an
approximation. Approximations for percentiles of the t-distribution have been
found with high levels of accuracy. OpenMC uses the `following approximation`_:

.. math::
    :label: percentile-n

    t_{x,n} = \sqrt{\frac{n}{n-2}} \left ( z_x + \frac{1}{4} \frac{z_x^3 -
    3z_x}{n-2} + \frac{1}{96} \frac{5z_x^5 - 56z_x^3 + 75z_x}{(n-2)^2} +
    \frac{1}{384} \frac{3z_x^7 - 81z_x^5 + 417z_x^3 - 315z_x}{(n-2)^3} \right )

where :math:`z_x` is the :math:`x` percentile of the standard normal
distribution. In order to determine an arbitrary percentile of the standard
normal distribution, we use an `unpublished rational approximation`_. After
using the rational approximation, one iteration of Newton's method is applied to
improve the estimate of the percentile.

.. only:: html

   .. rubric:: References

.. _following approximation: https://doi.org/10.1080/03610918708812641

.. _Bessel's correction: https://en.wikipedia.org/wiki/Bessel's_correction

.. _random variable: https://en.wikipedia.org/wiki/Random_variable

.. _stochastic process: https://en.wikipedia.org/wiki/Stochastic_process

.. _independent, identically-distributed random variables: https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables

.. _law of large numbers: https://en.wikipedia.org/wiki/Law_of_large_numbers

.. _expected value: https://en.wikipedia.org/wiki/Expected_value

.. _converges in probability: https://en.wikipedia.org/wiki/Convergence_of_random_variables#Convergence_in_probability

.. _normal distribution: https://en.wikipedia.org/wiki/Normal_distribution

.. _converges in distribution: https://en.wikipedia.org/wiki/Convergence_of_random_variables#Convergence_in_distribution

.. _confidence intervals: https://en.wikipedia.org/wiki/Confidence_interval

.. _Student's t-distribution: https://en.wikipedia.org/wiki/Student%27s_t-distribution

.. _Cauchy distribution: https://en.wikipedia.org/wiki/Cauchy_distribution

.. _unpublished rational approximation: https://web.archive.org/web/20150926021742/http://home.online.no/~pjacklam/notes/invnorm/

.. _MC21: http://www.osti.gov/bridge/servlets/purl/903083-HT5p1o/903083.pdf
