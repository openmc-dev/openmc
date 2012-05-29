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
