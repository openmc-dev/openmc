.. _methods_subcritical-multiplication:

=======================================
Subcritical Multiplication Calculations
=======================================

Subcritical multiplication problems are fixed source simulations with a large
amount of multiplication that are still subcritical with respect to
:math:`k_{eff}`. These problems are common in the design of accelerator driven
systems (ADS) which use an accelerator source to drive a subcritical
multiplication reaction in, e.g., spent fuel to transmute nuclear waste  and
even generate power [Bowman]_. An ADS is inherently safe and allows for a much
more flexible fuel composition, hence their popularity for waste transmutation.

For ADS's, the production of fission neutrons is central as these produced
neutrons amplify the external source. For the case of a proton spallation ADS,
source neutrons are produced from spallation reactions in a heavy metal target
bombarded by high-energy protons. These source neutrons then induce fission in
a subcritical core, producing additional neutrons. 


.. _methods_subcritical-multiplication-factors:

----------------------------------
Subcritical Multiplication Factors
----------------------------------
In a fixed source simulation, the total neutron production (per source particle) is used to define

.. math::
    :label: integral_multiplicity

    \begin{align*}
        M - 1 &= \frac{1}{N_{source}}\int d\boldsymbol{r} \int d\boldsymbol{\Omega} \int dE \nu \Sigma_f(\boldsymbol{r}, E) \psi(\boldsymbol{r}, \boldsymbol{\Omega}, E)\\
        &\coloneqq k + k^2 + k^3 + ... \\
        &= \frac{k}{1-k}\\
        \implies M &= \frac{1}{1-k}
    \end{align*}
   
Where :math:`M` is the subcritical multiplicity, and :math:`k` the subcritical
multiplication factor. The identification on the second line comes from the
picture of a single source neutron inducing several generations of fission
neutrons, producing on average :math:`k` neutrons in the first generation, 
which in turn produce :math:`k^2` neutrons in the second generation, and so on.
However, the above picture cannot be taken literally, because the neutrons
born from the external source will have a different importance to neutron
production than will fission neutrons, and we have the following alternative
picture for :math:`M` [Kobayashi]_:

.. math::
    :label: subcritical_k_factors
    
    \begin{align*}
        M-1 &= k_q + k_q k_s + k_q k_s^2 + ... \\
        &= \frac{k_q}{1-k_s} \\
        \implies \frac{k}{1-k} &= \frac{k_q}{1-k_s}\\
        k &= \frac{k_q}{1 - k_s + k_q}
    \end{align*}
Where :math:`k_q` is the multiplication factor of source neutrons, and :math:`k_s`
is the multiplication factor of fission neutrons, which together define an overall
subcritical multiplication factor :math:`k`. From the above it is clear that 
:math:`k < 1 \iff k_s < 1`, and :math:`k_s <1` for :math:`k_{eff}<1`, so a
subcritical system will correctly have :math:`k < 1`. It is however not the case
that :math:`k_s = k_{eff}`, because the angular flux of fission neutrons is not
necessarily the fundamental mode calculated in eigenvalue mode, nor is it
true that :math:`k = k_{eff}`. In fact, for deeply subcritical systems, 
:math:`k_{eff}` generally underestimates :math:`k` [Forget]_. In addition, we may
have :math:`k_q>1` despite :math:`k_s<1`, and in the limit, we may have :math:`k`
arbitrarily close to 1: an arbitrarily multiplying system that is still
subcritical with respect to fission neutrons. In fact, this is a primary design
consideration of ADS's, where :math:`k_s` is fixed :math:`<1` to ensure
subcritciality, while :math:`k_q` is maximized to achieve optimal multiplication.
It is therefore necessary to perform fixed source simulations to accurately determine 
subcritical multiplication and the flux distribution in ADS's.

.. _methods_subcritical-multiplication-methods:
----------------------------------
Subcritical Multiplication Methods
----------------------------------
The most straightforward method for performing a subcritical multiplication calculation is to perform a fixed source
simulation. For near critical systems, however, excessive secondary particle production can result in memory issues if
not properly managed. In very near critical systems, performing an eigenvalue simulation may be sufficient,
as the flux distribution will be close to the fundamental mode, and :math:`k` will be close to :math:`k_{eff}`, just note that
tally results must be multiplied by :math:`1/(1-k)`. For more deeply subcritical systems, there are several tailor-made methods
for subcritical multiplication calculations that limit excessive secondary particle production while equivalent results to fixed source simulations.
The method currently implemented in OpenMC is that of [Forget]_, which uses a modified power iteration to converge a source
distribution taking into account the external source and all generations of fission neutrons, which is then sampled from in generations
to accumulate statistics about the system. This method resembles the standard eigenvalue simulation (see :ref:`methods_eigenvalue`) in every way except 
that tallies have a pre-defined normalization by :math:`1/(1-k)` (rather than normalized to an arbitrary power distribution), the :math:`k`-eigenvalue is the 
subcritical multiplication factor (rather than the fundamental mode eigenvalue :math:`k_{eff}`), and the results are determined by the prescribed 
source distribution.


.. [Bowman] Bowman, Charles D. "Accelerator-driven systems for nuclear waste 
   transmutation." *Annual Review of Nuclear and Particle Science* 48.1 (1998): 
   505-556.
.. [Kobayashi] Kobayashi, Keisuke, and Kenji Nishihara. "Definition of 
   subcriticality using the importance function for the production of fission 
   neutrons." *Nuclear science and engineering* 136.2 (2000): 272-281.
.. [Forget] Forget, Benoit. "An Efficient Subcritical Multiplication Mode for 
   Monte Carlo Solvers." *Nuclear Science and Engineering* (2025): 1-11.