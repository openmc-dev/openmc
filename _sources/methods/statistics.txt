.. _methods_statistics:

==========
Statistics
==========

As was discussed briefly in :ref:`methods_introduction`, any given result from a
Monte Carlo calculation, colloquially known as a "tally", represents an estimate
of the mean of some random variable of interest. This random variable typically
corresponds to some physical quantity like a reaction rate, a net current across
some surface, or the neutron flux in a region. Given that all tallies are
produced by a stochastic process, there is an associated uncertainty with each
value reported. It is important to understand how the uncertainty is calculated
and what it tells us about our results. To that end, we will introduce a number
of theorems and results from statistics that should shed some light on the
interpretation of uncertainties.

--------------------
Law of Large Numbers
--------------------

Let :math:`X_1, X_2, \dots, X_n` be an infinite sequence of independent,
identically-distributed random variables with expected values :math:`E(X_1) =
E(X_2) = \mu`. The sample mean :math:`\bar{X_n} = \frac{X_1 + \dots + X_n}{n}`
converges in probability to the true mean, i.e. for all :math:`\epsilon > 0`

.. math::

    \lim\limits_{n\rightarrow\infty} P \left ( \left | \bar{X}_n - \mu \right |
    \ge \epsilon \right ) = 0.

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
generated with a method chosen at random". Some theory should be used."
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

.. _linear congruential generator: http://en.wikipedia.org/wiki/Linear_congruential_generator

.. _Brown: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/anl_rn_arb-strides_1994.pdf
