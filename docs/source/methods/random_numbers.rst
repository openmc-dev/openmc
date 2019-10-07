.. _methods_random_numbers:

========================
Random Number Generation
========================

In order to sample probability distributions, one must be able to produce random
numbers. The standard technique to do this is to generate numbers on the
interval :math:`[0,1)` from a deterministic sequence that has properties that
make it appear to be random, e.g. being uniformly distributed and not exhibiting
correlation between successive terms. Since the numbers produced this way are
not truly "random" in a strict sense, they are typically referred to as
pseudorandom numbers, and the techniques used to generate them are pseudorandom
number generators (PRNGs). Numbers sampled on the unit interval can then be
transformed for the purpose of sampling other continuous or discrete probability
distributions.

------------------------------
Linear Congruential Generators
------------------------------

There are a great number of algorithms for generating random numbers. One of the
simplest and commonly used algorithms is called a `linear congruential
generator`_. We start with a random number *seed* :math:`\xi_0` and a sequence
of random numbers can then be generated using the following recurrence relation:

.. math::
    :label: lcg

    \xi_{i+1} = g \xi_i + c \mod M

where :math:`g`, :math:`c`, and :math:`M` are constants. The choice of these
constants will have a profound effect on the quality and performance of the
generator, so they should not be chosen arbitrarily. As Donald Knuth stated in
his seminal work *The Art of Computer Programming*, "random numbers should not
be generated with a method chosen at random. Some theory should be used."
Typically, :math:`M` is chosen to be a power of two as this enables :math:`x
\mod M` to be performed using the bitwise AND operator with a bit mask. The
constants for the linear congruential generator used by default in OpenMC are
:math:`g = 2806196910506780709`, :math:`c = 1`, and :math:`M = 2^{63}` (see
`L'Ecuyer`_).

Skip-ahead Capability
---------------------

One of the important capabilities for a random number generator is to be able to
skip ahead in the sequence of random numbers. Without this capability, it would
be very difficult to maintain reproducibility in a parallel calculation. If we
want to skip ahead :math:`N` random numbers and :math:`N` is large, the cost of
sampling :math:`N` random numbers to get to that position may be prohibitively
expensive. Fortunately, algorithms have been developed that allow us to skip
ahead in :math:`O(\log_2 N)` operations instead of :math:`O(N)`. One algorithm
to do so is described in a paper by Brown_. This algorithm relies on the following
relationship:

.. math::
    :label: lcg-skipahead

    \xi_{i+k} = g^k \xi_i + c \frac{g^k - 1}{g - 1} \mod M

Note that equation :eq:`lcg-skipahead` has the same general form as equation :eq:`lcg`, so
the idea is to determine the new multiplicative and additive constants in
:math:`O(\log_2 N)` operations.

.. only:: html

   .. rubric:: References


.. _L'Ecuyer: https://doi.org/10.1090/S0025-5718-99-00996-5
.. _Brown: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/anl-rn-arb-stride.pdf
.. _linear congruential generator: https://en.wikipedia.org/wiki/Linear_congruential_generator
