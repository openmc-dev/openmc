.. _methods_introduction:

============
Introduction
============

The physical process by which a population of particles evolves over time is
governed by a number of `probability distributions`_. For instance, given a
particle traveling through some material, there is a probability distribution
for the distance it will travel until its next collision (an exponential
distribution). Then, when it collides with a nucleus, there is an associated
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

------------------------
Overview of Program Flow
------------------------

OpenMC performs a Monte Carlo simulation one particle at a time -- at no point
is more than one particle being tracked on a single program instance. Before any
particles are tracked, the problem must be initialized. This involves the
following steps:

  - Read input files and building data structures for the geometry, materials,
    tallies, and other associated variables.

  - Initialize the pseudorandom number generator.

  - Read the contiuous-energy or multi-group cross section data specified in
    the problem.

  - If using a special energy grid treatment such as a union energy grid or
    lethargy bins, that must be initialized as well in a continuous-energy
    problem.

  - In a multi-group problem, individual nuclide cross section information is
    combined to produce material-specific cross section data.

  - In a fixed source problem, source sites are sampled from the specified
    source. In an eigenvalue problem, source sites are sampled from some
    initial source distribution or from a source file. The source sites
    consist of coordinates, a direction, and an energy.

Once initialization is complete, the actual transport simulation can
proceed. The life of a single particle will proceed as follows:

  1. The particle's properties are initialized from a source site previously
     sampled.

  2. Based on the particle's coordinates, the current cell in which the particle
     resides is determined.

  3. The energy-dependent cross sections for the material that the particle is
     currently in are determined. Note that this includes the total
     cross section, which is not pre-calculated.

  4. The distance to the nearest boundary of the particle's cell is determined
     based on the bounding surfaces to the cell.

  5. The distance to the next collision is sampled. If the total material
     cross section is :math:`\Sigma_t`, this can be shown to be

     .. math::

         d = -\frac{\ln \xi}{\Sigma_t}

     where :math:`\xi` is a `pseudorandom number`_ sampled from a uniform
     distribution on :math:`[0,1)`.

  6. If the distance to the nearest boundary is less than the distance to the next
     collision, the particle is moved forward to this boundary. Then, the process
     is repeated from step 2. If the distance to collision is closer than the
     distance to the nearest boundary, then the particle will undergo a collision.

  7. The material at the collision site may consist of multiple nuclides. First,
     the nuclide with which the collision will happen is sampled based on the
     total cross sections. If the total cross section of material :math:`i` is
     :math:`\Sigma_{t,i}`, then the probability that any nuclide is sampled is

     .. math::

         P(i) = \frac{\Sigma_{t,i}}{\Sigma_t}.

     Note that the above selection of collided nuclide only applies to
     continuous-energy simulations as multi-group simulations use nuclide
     data which has already been combined in to material-specific data.

  8. Once the specific nuclide is sampled, the random samples a reaction for
     that nuclide based on the microscopic cross sections. If the microscopic
     cross section for some reaction :math:`x` is :math:`\sigma_x` and the total
     microscopic cross section for the nuclide is :math:`\sigma_t`, then the
     probability that reaction :math:`x` will occur is

     .. math::

         P(x) = \frac{\sigma_x}{\sigma_t}.

     Since multi-group simulations use material-specific data, the above is
     performed with those material multi-group cross sections (i.e.,
     macroscopic cross sections for the material) instead of microscopic
     cross sections for the nuclide).

  9. If the sampled reaction is elastic or inelastic scattering, the outgoing
     energy and angle is sampled from the appropriate distribution.  In
     continuous-energy simulation, reactions of type :math:`(n,xn)` are treated
     as scattering and any additional particles which may be created are added
     to a secondary particle bank to be tracked later. In a multi-group
     simulation, this secondary bank is not used but the particle weight is
     increased accordingly.  The original particle then continues from step 3.
     If the reaction is absorption or fission, the particle dies and if
     necessary, fission sites are created and stored in the fission bank.

After all particles have been simulated, there are a few final tasks that must
be performed before the run is finished. This include the following:

  - With the accumulated sum and sum of squares for each tally, the sample mean
    and its variance is calculated.

  - All tallies and other results are written to disk.

  - If requested, a source file is written to disk.

  - Dynamically-allocated memory should be freed.

.. _probability distributions: https://en.wikipedia.org/wiki/Probability_distribution
.. _Monte Carlo: https://en.wikipedia.org/wiki/Monte_Carlo_method
.. _central limit theorem: https://en.wikipedia.org/wiki/Central_limit_theorem
.. _pseudorandom number: https://en.wikipedia.org/wiki/Pseudorandom_number_generator
