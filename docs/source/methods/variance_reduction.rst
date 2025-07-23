.. _methods_variance_reduction:

==================
Variance Reduction
==================

.. _methods_variance_reduction_intro:

------------
Introduction
------------

Transport problems can sometimes involve a significant degree of attenuation
between the source and a detector (tally) region, which can result in a flux
differential of ten orders of magnitude (or more) throughout the simulation
domain. As Monte Carlo uncertainties tend to be inversely proportional to the
physical flux density, it can be extremely difficult to accurately resolve
tallies in locations that are optically far from the source. This issue is
particularly common in fixed source simulations, where some tally locations may
not experience a single scoring event, even after billions of analog histories.

Variance reduction techniques aim to either flatten the global uncertainty
distribution, such that all regions of phase space have a fairly similar
uncertainty, or to reduce the uncertainty in specific locations (such as a
detector). There are three strategies available in OpenMC for variance reduction:
weight-windowing via the Monte Carlo MAGIC method or the FW-CADIS method, and
source biasing. Both weight windowing strategies work by developing a mesh that 
can be utilized by subsequent Monte Carlo solves to split particles heading 
towards areas of lower flux densities while terminating particles in higher flux 
regions. In contrast, source biasing modifies source site sampling behavior to 
preferentially track particles more likely to reach phase space regions of 
interest.


------------
MAGIC Method
------------

The Method of Automatic Generation of Importances by Calculation, or `MAGIC
method <https://doi.org/10.1016/j.fusengdes.2011.01.059>`_, is an iterative
technique that uses spatial flux information :math:`\phi(r)` obtained from a
normal Monte Carlo solve to produce weight windows :math:`w(r)` that can be
utilized by a subsequent iteration of Monte Carlo. While the first generation of
weight windows produced may only help to reduce variance slightly, use of these
weights to generate another set of weight windows results in a progressively
improving iterative scheme.

Equation :eq:`magic` defines how the lower bound of weight windows
:math:`w_{\ell}(r)` are generated with MAGIC using forward flux information.
Here, we can see that the flux at location :math:`r` is normalized by the
maximum flux in any group at that location. We can also see that the weights are
divided by a factor of two, which accounts for the typical :math:`5\times`
factor separating the lower and upper weight window bounds in OpenMC.

.. math::
    :label: magic

    w_{\ell}(r) = \frac{\phi(r)}{2\,\text{max}(\phi(r))}

A major advantage of this technique is that it does not require any special
transport machinery; it simply uses multiple Monte Carlo simulations to
iteratively improve a set of weight windows (which are typically defined on a
mesh covering the simulation domain). The downside to this method is that as the
flux differential increases between areas near and far from the source, it
requires more outer Monte Carlo iterations, each of which can be expensive in
itself. Additionally, computation of weight windows based on regular (forward)
neutron flux tally information does not produce the most numerically effective
set of weight windows. Nonetheless, MAGIC remains a simple and effective
technique for generating weight windows.

--------
FW-CADIS
--------

As discussed in the previous section, computation of weight windows based on
regular (forward) neutron flux tally information does not produce the most
numerically efficient set of weight windows. It is highly preferable to generate
weight windows based on spatial adjoint flux :math:`\phi^{\dag}(r)`
information. The adjoint flux is essentially the "reverse" simulation problem,
where we sample a random point and assume this is where a particle was absorbed,
and then trace it backwards (upscattering in energy), until we sample the point
where it was born from.

The Forward-Weighted Consistent Adjoint Driven Importance Sampling method, or
`FW-CADIS method <https://doi.org/10.13182/NSE12-33>`_, produces weight windows
for global variance reduction given adjoint flux information throughout the
entire domain. The weight window lower bound is defined in Equation
:eq:`fw_cadis`, and also involves a normalization step not shown here.

.. math::
    :label: fw_cadis

    w_{\ell}(r) = \frac{1}{2\phi^{\dag}(r)}

While the algorithm itself is quite simple, it requires estimates of the global
adjoint flux distribution, which is difficult to generate directly with Monte
Carlo transport. Thus, FW-CADIS typically uses an alternative solver (often
deterministic) that can be more readily adapted for generating adjoint flux
information, and which is often much cheaper than Monte Carlo given that a rough
solution is often sufficient for weight window generation.

The FW-CADIS implementation in OpenMC utilizes its own internal random ray
multigroup transport solver to generate the adjoint source distribution. No
coupling to any external transport is solver is necessary. The random ray solver
operates on the same geometry as the Monte Carlo solver, so no redefinition of
the simulation geometry is required. More details on how the adjoint flux is
computed are given in the :ref:`adjoint methods section <adjoint>`.

More information on the workflow is available in the :ref:`user guide
<variance_reduction>`, but generally production of weight windows with FW-CADIS
involves several stages (some of which are highly automated). These tasks
include generation of approximate multigroup cross section data for use by the
random ray solver, running of the random ray solver in normal (forward flux)
mode to generate a source for the adjoint solver, running of the random ray
solver in adjoint mode to generate adjoint flux tallies, and finally the
production of weight windows via the FW-CADIS method. As is discussed in the
user guide, most of these steps are automated together, making the additional
burden on the user fairly small.

The major advantage of this technique is that it typically produces much more
numerically efficient weight windows as compared to those generated with MAGIC,
sometimes with an order-of-magnitude improvement in the figure of merit
(Equation :eq:`variance_fom`), which accounts for both the variance and the
execution time. Another major advantage is that the cost of the random ray
solver is typically negligible compared to the cost of the subsequent Monte
Carlo solve itself, making it a very cheap method to deploy. The downside to
this method is that it introduces a second transport method into the mix (random
ray), such that there are more free input parameters for the user to know about
and adjust, potentially making the method more complex to use. However, as many
of the parameters have natural choices, much of this parameterization can be
handled automatically behind the scenes without the need for the user to be
aware of this.

.. math::
    :label: variance_fom

    \text{FOM} = \frac{1}{\text{Time} \times \sigma^2}

.. _methods_source_biasing:

--------------
Source Biasing
--------------

In contrast to the previous two methods which introduce population controls 
during transport, source biasing modifies the sampling of fixed source site 
distributions. The basic premise of the technique is that for each spatial, 
angular (direction), energy, or time distribution of a fixed source, an 
additional distribution can be specified, provided that the two share a common 
support. Samples are then drawn from this "bias" distribution, which can be 
chosen to preferentially direct particles towards phase space regions of 
interest. In order to avoid biasing the tally results, however, a weight 
adjustment is applied to each sampled site as described below.

Assume that the unbiased probability density function of a random variable 
:math:`X:x \rightarrow \mathbb{R}` is given by :math:`p(x)`, but that using the 
biasing distribution :math:`b(x)` will result in a greater number of particle 
trajectories reaching some phase space region of interest. Then a sample 
:math:`x_0` may be drawn from :math:`b(x)` while maintaining a fair game, 
provided that its weight is adjusted according to Equation :eq:`source_bias`:: 


.. math::
    :label: source_bias

    w = w_0 \times \frac{p(x_0)}{b(x_0)} 

Here, :math:`w_0` is the weight of an unbiased sample from :math:`p(x)`, 
typically unity.

When an independent source is sampled in OpenMC, the particle's coordinate in 
each variable of phase space :math:`(\mathbf{r},\mathbf{\Omega},E,t)` is 
successively drawn from an independent probability distribution. Multiple 
variables can be biased, in which case the resultant weight :math:`w` applied to 
the particle is the product of the weights assigned from all sampled 
distributions: space, angle, energy, and time, as shown in Equation 
:eq:`tot_wgt`.

.. math::
    :label: tot_wgt

    w = w_r \times w_{\Omega} \times w_E \times w_t 

Returning now to Equation :eq:`source_bias`, the requirement for common support 
becomes evident. If :math:`\mathrm{supp} b(x)` fully contains but is not 
identical to :math:`\mathrm{supp} p(x)`, then some samples from :math:`b(x)` 
will correspond to points where :math:`p(x) = 0`. Thus these source sites would 
be assigned a starting weight of 0, meaning the particles would be killed 
immediately upon transport, yet would still count towards the total specified 
in :attr:`Settings.particles` (hence biasing tallies and increasing variance). 
Conversely, if :math:`\mathrm{supp} b(x)` is fully contained by but not 
identical to :math:`\mathrm{supp} p(x)`, the contributions of some regions 
outside :math:`\mathrm{supp} b(x)` will not be counted towards the integral, 
potentially biasing the tally. The weight assigned to such points 
:math:`\mathbf{x_i}` would be undefined since 
:math:`b(\mathbf{x_i}) = \mathbf{0}`.

Finally, in addition to offering global or local variance reduction 
capabilities, source biasing usually requires fewer additional lines to 
implement than FW-CADIS and MAGIC weight windowing in simple applications. In 
comparison to these techniques, it also permits continuous biasing in spatial, 
angle, energy, and time dimensions, and does not require additional transport 
calculations. However, as particle transport proceeds as usual after a biased 
source is sampled, particle attenuation in optically thick regions outside the 
source volume will not be affected by source biasing. Instead, transport 
biasing techniques such as weight windows will be more useful in such 
scenarios.
