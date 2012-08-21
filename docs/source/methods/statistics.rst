.. _methods_statistics:

==========
Statistics
==========

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

--------------------
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

---------------------
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

------------------------------------------
Estimating Statistics of a Random Variable
------------------------------------------

Mean
----

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
--------

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
--------------------

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
--------------------

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

One should be cautioned that equation :eq:`two-sided-ci` **only applies if the
underlying random variables are normally-distributed!** In general, this may not
be true for a tally random variable -- the central limit theorem guarantees only
that the sample mean is normally distributed, not the underlying random
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
found with high levels of accuracy. OpenMC uses the approximation from
[George]_:

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

------------------------
Random Number Generation
------------------------

In order to sample probability distributions, one must be able to produce random
numbers. The standard technique to do this is to generate numbers on the
interval :math:`[0,1)` from a deterministic sequence that has a properties that
make it appear to be random, e.g. being uniformly distributed and not exhibiting
correlation between successive terms. Since the numbers are not truly "random"
in the strict sense, they are typically referred to as pseudo-random numbers,
and the techniques used to generate them are pseudo-random number generators
(PRNGs). Numbers sampled on the unit interval can then be used transformed for
the purpose of sampling other continuous or discrete probability distributions.

There are a great number of algorithms for generating random numbers. One of the
simplest and commonly used algorithms is called a `linear congruential
generator`_. We start with some random number seed :math:`\xi_0` and a sequence
of random numbers is generated using the following recurrence relation:

.. math::
    :label: lcg

    \xi_{i+1} = g \xi_i + c \mod M

where :math:`g`, :math:`c`, and :math:`M` are constants. The choice of these
constants will have a profound effect on the quality and performance of the
generator, so they should not be chosen arbitrarily. As Donald Knuth said in his
seminal work *The Art of Computer Programming*, "random numbers should not be
generated with a method chosen at random. Some theory should be used."
Typically, :math:`M` is chosen to be a power of two as this enables :math:`x
\mod M` to be performed using the binary AND operator with a bit mask. The
constants for the linear congruential generator used by default in OpenMC are
:math:`g = 2806196910506780709`, :math:`c = 1`, and :math:`M = 2^{63}`.

One of the important capabilities for a random number generator is to be able to
skip ahead in the sequence of random numbers. Without this capability, it would
be very difficult to maintain reproducibility in a parallel calculation. If we
want to skip ahead :math:`N` random numbers and :math:`N` is large, the cost of
just sampling :math:`N` random numbers to get to that position may be
prohibitively expensive. Fortunately, algorithms have been developed that allow
us to skip ahead in :math:`O(\log N)` operations instead of :math:`O(N)`. One
algorithm to do so is described in a paper by Brown_. This algorithm relies on
the following relationship:

.. math::
    :label: lcg-skipahead

    \xi_{i+k} = g^k \xi_i + c \frac{g^k - 1}{g - 1} \mod M

Note that equation :eq:`lcg-skipahead` has the same form as equation :eq:`lcg`
so the idea is to determine the new multiplicative and additive constants in
:math:`O(\log N)` operations.

----------
References
----------

.. [George] E. E. Olusegun George and Meenakshi Sivaram, "A modification of the
   Fisher-Cornish approximation for the student t percentiles," Communication
   in Statistics - Simulation and Computation, 16 (4), pp. 1123-1132 (1987).

.. _linear congruential generator: http://en.wikipedia.org/wiki/Linear_congruential_generator

.. _Brown: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/anl_rn_arb-strides_1994.pdf

.. _Bessel's correction: http://en.wikipedia.org/wiki/Bessel's_correction

.. _random variable: http://en.wikipedia.org/wiki/Random_variable

.. _stochastic process: http://en.wikipedia.org/wiki/Stochastic_process

.. _independent, identically-distributed random variables: http://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables

.. _law of large numbers: http://en.wikipedia.org/wiki/Law_of_large_numbers

.. _expected value: http://en.wikipedia.org/wiki/Expected_value

.. _converges in probability: http://en.wikipedia.org/wiki/Convergence_of_random_variables#Convergence_in_probability

.. _normal distribution: http://en.wikipedia.org/wiki/Normal_distribution

.. _converges in distribution: http://en.wikipedia.org/wiki/Convergence_of_random_variables#Convergence_in_distribution

.. _confidence intervals: http://en.wikipedia.org/wiki/Confidence_interval

.. _Student's t-distribution: http://en.wikipedia.org/wiki/Student%27s_t-distribution

.. _Cauchy distribution: http://en.wikipedia.org/wiki/Cauchy_distribution

.. _unpublished rational approximation: http://home.online.no/~pjacklam/notes/invnorm/
