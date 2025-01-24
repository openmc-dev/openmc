.. _methods_variance_reduction:

==================
Variance Reduction
==================

.. _methods_variance_reduction_intro:

------------
Introduction
------------

Transport problems can sometimes involve a significant degree of attenuation
between the neutron source and a detector (tally) region, which can result in a
flux differential of ten orders of magnitude (or more) throughout the simulation
domain. As Monte Carlo uncertainties tend to be inversely proportional to the
physical flux density, it can be extremely difficult to accurately resolve
tallies in locations that are optically far from the source. This issue is
particularly common in fixed source simulations, where some tally locations may
not experience a single scoring event, even after billions of analog histories.

Variance reduction techniques aim to either flatten the global uncertainty
distribution, such that all regions of phase space have a fairly similar
uncertainty, or to reduce the uncertainty in specific locations (such as a
detector). There are two strategies available in OpenMC for variance reduction:
the Monte Carlo MAGIC method, and the FW-CADIS method. Both strategies work by
developing a weight window mesh, which can be utilized by subsequent Monte Carlo
solves to split particles heading towards areas of lower flux densities while
terminating particles in higher flux regions -- all while maintaining a fair
game.

------------
MAGIC Method
------------

The MAGIC method is an iterative technique that uses spatial flux information
:math:`\phi(r)` obtained from a normal Monte Carlo solve to produce weight
windows :math:`w(r)` that can be utilized by a subsequent iteration of Monte
Carlo. While the first generation of weight windows produced may only help to
reduce variance slightly, use of these weights to generate another set of weight
windows results in a progressively improving iterative scheme. 

Equation :eq:`magic` defines how the lower bound of weight windows
:math:`w_{\ell}(r)` are generated with MAGIC using forward flux information.
Here, we can see that the flux at location :math:`r` is normalized by the
maximum flux in any group at that location. We can also see that the weights are
divided by a factor of two, which accounts for the typical :math:`5\times`
factor separating the lower and upper weight window bounds in OpenMC.

.. math::
    :label: magic

    w_{\ell}(r) = \frac{\phi(r)}{2\text{max}(\phi(r))}

A major advantage of this technique is that it does not require any special
transport machinery -- it simply uses multiple Monte Carlo simulations to
iteratively improve a set of weight windows (which are typically defined on a
mesh covering the simulation domain). The downside to this method is that as the
flux differential increases between areas near and far from the source, it
requires more outer Monte Carlo iterations, each of which can be 
expensive in itself. Additionally, computation of weight windows based on
regular (forward) neutron flux tally information does not produce the most
numerically effective set of weight windows. Nonetheless, MAGIC remains a simple
and effective technique for generating weight windows.

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

The FW-CADIS method produces weight windows for global variance reduction given
adjoint flux information throughout the entire domain. The weight window lower
bound is defined in Equation :eq:`fw_cadis`, and also involves a normalization
step not shown here.

.. math::
    :label: fw_cadis

    w_{\ell}(r) = \frac{1}{2\phi^{\dag}(r)}

While the algorithm itself is quite simple, it requires estimates of the global
adjoint flux distribution, which is difficult to generate directly with Monte
Carlo transport. Thus, FW-CADIS typically uses an alternative solver (often
deterministic) that can be more readily adapted for generating adjoint flux
information, and which is often much cheaper than Monte Carlo given that
maximal-fidelity is not needed for weight window generation.

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
sometimes with an improvement on the variance vs. runtime figure of merit
(Equation :eq:`variance_fom`) of an order of magnitude. Another major advantage
is that the cost of the random ray solver is typically negligible compared to
the cost of the subsequent Monte Carlo solve itself, making it a very cheap
method to deploy. The downside to this method is that it introduces a second
transport method into the mix (random ray), such that there are more free input
parameters for the user to know about and adjust, potentially making the method
more complex to use. However, as many of the parameters have natural choices,
much of this parameterization can be handled automatically behind the scenes
without the need for the user to be aware of this.

.. math::
    :label: variance_fom

    \text{FOM} = \frac{1}{\text{Time} \times \sigma^2}