.. _methods_introduction:

============
Introduction
============

The physical process by which a population of particles evolves over time is
governed by a number of `probability distributions`_. For instance, given a
particle traveling through some material, there is a probability distribution
for the distance it will travel until its next collision (an exponential
distribution). Then, when it collides with a nucleus, there is associated
probability of undergoing each possible reaction with that nucleus. While the
behavior of any single particle is unpredictable, the average behavior of a
large population of particles originating from the same source is well defined.

If the probability distributions that govern the transport of a particle are
known, the process of single particles randomly streaming and colliding with
nuclei can be simulated directly with computers using a technique known as
`Monte Carlo`_ simulation. If enough particles are simulated this way, the
average behavior can be determined to within arbitrarily small statistical
error, a fact guaranteed by the `central limit theorem`_. To be more precise,
the central limit theorem tells us that the variance of the sample mean of some
physical parameter being estimated with Monte Carlo will be inversely
proportional to the number of realizations, i.e. the number of particles we
simulate:

.. math::

    \sigma^2 \propto \frac{1}{N}.

where :math:`\sigma^2` is the variance of the sample mean and :math:`N` is the
number of realizations.

---------------------
Criticality Algorithm
---------------------

.. _probability distributions: http://en.wikipedia.org/wiki/Probability_distribution
.. _Monte Carlo: http://en.wikipedia.org/wiki/Monte_Carlo_method
.. _central limit theorem: http://en.wikipedia.org/wiki/Central_limit_theorem
