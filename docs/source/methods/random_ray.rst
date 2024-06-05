.. _methods_random_ray:

==========
Random Ray
==========

.. _methods_random_ray_intro:

-------------------
What is Random Ray?
-------------------

`Random ray <Tramm-2017a>`_ is a stochastic transport method, closely related to
the deterministic Method of Characteristics (MOC) [Askew-1972]_. Rather than
each ray representing a single neutron as in Monte Carlo, it represents a
characteristic line through the simulation geometry upon which the transport
equation can be written as an ordinary differential equation that can be solved
analytically (although with discretization required in energy, making it a
multigroup method). The behavior of the governing transport equation can be
approximated by solving along many characteristic tracks (rays) through the
system. Unlike particles in Monte Carlo, rays in random ray or MOC are not
affected by the material characteristics of the simulated problem---rays are
selected so as to explore the full simulation problem with a statistically equal
distribution in space and angle.

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/pHQq3FE4PDo?si=kPm9ngMBr95wLRGC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>

The above animation is an example of the random ray integration process at work,
showing a series of random rays being sampled and transported through the
geometry. In the following sections, we will discuss how the random ray solver
works.

----------------------------------------------
Why is a Random Ray Solver Included in OpenMC?
----------------------------------------------

* One area that Monte Carlo struggles with is maintaining numerical efficiency
  in regions of low physical particle flux. Random ray, on the other hand, has
  approximately even variance throughout the entire global simulation domain,
  such that areas with low neutron flux are no less well known that areas of
  high neutron flux. Absent weight windows in MC, random ray can be several
  orders of magnitude faster than multigroup Monte Carlo in classes of problems
  where areas with low physical neutron flux need to be resolved. While MC
  uncertainty can be greatly improved with variance reduction techniques, they
  add some user complexity, and weight windows can often be expensive to
  generate via MC transport alone (e.g., via the `MAGIC method
  <https://doi.org/10.1016/j.fusengdes.2011.01.059>`_). The random ray solver
  may be used in future versions of OpenMC as a fast way to generate weight
  windows for subsequent usage by the MC solver in OpenMC.

* In practical implementation terms, random ray is mechanically very similar to
  how Monte Carlo works, in terms of the process of ray tracing on constructive
  solid geometry (CSG) and handling stochastic convergence, etc. In the original
  1972 paper by Askew that introduces MOC (which random ray is a variant of), he
  stated:

    .. epigraph::

        "One of the features of the method proposed [MoC] is that ... the
        tracking process needed to perform this operation is common to the
        proposed method ... and to Monte Carlo methods. Thus a single tracking
        routine capable of recognizing a geometric arrangement could be utilized
        to service all types of solution, choice being made depending which was
        more appropriate to the problem size and required accuracy."

        -- Askew [Askew-1972]_

    This prediction holds up---the additional requirements needed in OpenMC to
    handle random ray transport turned out to be fairly small.

* It amortizes the code complexity in OpenMC for representing multigroup cross
  sections. There is a significant amount of interface code, documentation, and
  complexity in allowing OpenMC to generate and use multigroup XS data in its
  MGMC mode. Random ray allows the same multigroup data to be used, making full
  reuse of these existing capabilities.

-------------------------------
Random Ray Numerical Derivation
-------------------------------

In this section, we will derive the numerical basis for the random ray solver
mode in OpenMC. The derivation of random ray is also discussed in several papers
(`1 <Tramm-2017a>`_, `2 <Tramm-2017b>`_, `3 <Tramm-2018>`_), and some of those
derivations are reproduced here verbatim. Several extensions are also made to
add clarity, particularly on the topic of OpenMC's treatment of cell volumes in
the random ray solver.

~~~~~~~~~~~~~~~~~~~~~~~~~
Method of Characteristics
~~~~~~~~~~~~~~~~~~~~~~~~~

The Boltzmann neutron transport equation is a partial differential equation
(PDE) that describes the angular flux within a system. It is a balance equation,
with the streaming and absorption terms typically appearing on the left hand
side, which are balanced by the scattering source and fission source terms on
the right hand side.

.. math::
    :label: transport

    \begin{align*}
    \mathbf{\Omega} \cdot \mathbf{\nabla} \psi(\mathbf{r},\mathbf{\Omega},E) & + \Sigma_t(\mathbf{r},E) \psi(\mathbf{r},\mathbf{\Omega},E) = \\
    & \int_0^\infty d E^\prime \int_{4\pi} d \Omega^{\prime} \Sigma_s(\mathbf{r},\mathbf{\Omega}^\prime \rightarrow \mathbf{\Omega}, E^\prime \rightarrow E) \psi(\mathbf{r},\mathbf{\Omega}^\prime, E^\prime) \\
    & + \frac{\chi(\mathbf{r}, E)}{4\pi k_{eff}} \int_0^\infty dE^\prime \nu \Sigma_f(\mathbf{r},E^\prime) \int_{4\pi}d \Omega^\prime \psi(\mathbf{r},\mathbf{\Omega}^\prime,E^\prime)
    \end{align*}

In Equation :eq:`transport`, :math:`\psi` is the angular neutron flux. This
parameter represents the total distance traveled by all neutrons in a particular
direction inside of a control volume per second, and is often given in units of
:math:`1/(\text{cm}^{2} \text{s})`. As OpenMC does not support time dependence
in the random ray solver mode, we consider the steady state equation, where the
units of flux become :math:`1/\text{cm}^{2}`. The angular direction unit vector,
:math:`\mathbf{\Omega}`, represents the direction of travel for the neutron. The
spatial position vector, :math:`\mathbf{r}`,  represents the location within the
simulation. The neutron energy, :math:`E`, or speed in continuous space, is
often given in units of electron volts. The total macroscopic neutron cross
section is :math:`\Sigma_t`. This value represents the total probability of
interaction between a neutron traveling at a certain speed (i.e., neutron energy
:math:`E`) and a target nucleus (i.e., the material through which the neutron is
traveling) per unit path length, typically given in units of
:math:`1/\text{cm}`. Macroscopic cross section data is a combination of
empirical data and quantum mechanical modeling employed in order to generate an
evaluation represented either in pointwise form or resonance parameters for each
target isotope of interest in a material, as well as the density of the
material, and is provided as input to a simulation. The scattering neutron cross
section, :math:`\Sigma_s`, is similar to the total cross section but only
measures scattering interactions between the neutron and the target nucleus, and
depends on the change in angle and energy the neutron experiences as a result of
the interaction. Several additional reactions like (n,2n) and (n,3n) are
included in the scattering transfer cross section. The fission neutron cross
section, :math:`\Sigma_f`, is also similar to the total cross section but only
measures the fission interaction between a neutron and a target nucleus. The
energy spectrum for neutrons born from fission, :math:`\chi`, represents a known
distribution of outgoing neutron energies based on the material that fissioned,
which is taken as input data to a computation. The average number of neutrons
born per fission is :math:`\nu`. The eigenvalue of the equation,
:math:`k_{eff}`, represents the effective neutron multiplication factor. If the
right hand side of Equation :eq:`transport` is condensed into a single term,
represented by the total neutron source term :math:`Q(\mathbf{r}, \mathbf{\Omega},E)`,
the form given in Equation :eq:`transport_simple` is reached.

.. math::
    :label: transport_simple

    \overbrace{\mathbf{\Omega} \cdot \mathbf{\nabla} \psi(\mathbf{r},\mathbf{\Omega},E)}^{\text{streaming term}} + \overbrace{\Sigma_t(\mathbf{r},E) \psi(\mathbf{r},\mathbf{\Omega},E)}^{\text{absorption term}} = \overbrace{Q(\mathbf{r}, \mathbf{\Omega},E)}^{\text{total neutron source term}}

Fundamentally, MOC works by solving Equation :eq:`transport_simple` along a
single characteristic line, thus altering the full spatial and angular scope of
the transport equation into something that holds true only for a particular
linear path (or track) through the reactor. These tracks are linear for neutral
particles that are not subject to field effects. With our transport equation in
hand, we will now derive the solution along a track. To accomplish this, we
parameterize :math:`\mathbf{r}` with respect to some reference location
:math:`\mathbf{r}_0` such that :math:`\mathbf{r} = \mathbf{r}_0 + s\mathbf{\Omega}`. In this
manner, Equation :eq:`transport_simple` can be rewritten for a specific segment
length :math:`s` at a specific angle :math:`\mathbf{\Omega}` through a constant
cross section region of the reactor geometry as in Equation :eq:`char_long`.

.. math::
    :label: char_long

    \mathbf{\Omega} \cdot \mathbf{\nabla} \psi(\mathbf{r}_0 + s\mathbf{\Omega},\mathbf{\Omega},E) + \Sigma_t(\mathbf{r}_0 + s\mathbf{\Omega},E) \psi(\mathbf{r}_0 + s\mathbf{\Omega},\mathbf{\Omega},E) = Q(\mathbf{r}_0 + s\mathbf{\Omega}, \mathbf{\Omega},E)

As this equation holds along a one dimensional path, we can assume the
dependence of :math:`s` on :math:`\mathbf{r}_0` and :math:`\mathbf{\Omega}` such that
:math:`\mathbf{r}_0 + s\mathbf{\Omega}` simplifies to :math:`s`. When the differential
operator is also applied to the angular flux :math:`\psi`, we arrive at the
characteristic form of the Boltzmann Neutron Transport Equation given in
Equation :eq:`char`.

.. math::
    :label: char

    \frac{d}{ds} \psi(s,\mathbf{\Omega},E) + \Sigma_t(s,E) \psi(s,\mathbf{\Omega},E) = Q(s, \mathbf{\Omega},E)

An analytical solution to this characteristic equation can be achieved with the
use of an integrating factor:

.. math::
    :label: int_factor

    e^{ \int_0^s ds' \Sigma_t (s', E)}

to arrive at the final form of the characteristic equation shown in Equation
:eq:`full_char`.

.. math::
    :label: full_char

    \psi(s,\mathbf{\Omega},E) = \psi(\mathbf{r}_0,\mathbf{\Omega},E) e^{-\int_0^s ds^\prime \Sigma_t(s^\prime,E)} + \int_0^s ds^{\prime\prime} Q(s^{\prime\prime},\mathbf{\Omega}, E) e^{-\int_{s^{\prime\prime}}^s ds^\prime \Sigma_t(s^\prime,E)}

With this characteristic form of the transport equation, we now have an
analytical solution along a linear path through any constant cross section
region of a system. While the solution only holds along a linear track, no
discretizations have yet been made.

Similar to many other solution approaches to the Boltzmann neutron transport
equation, the MOC approach also uses a "multigroup" approximation in order to
discretize the continuous energy spectrum of neutrons traveling through the
system into fixed set of energy groups :math:`G`, where each group :math:`g \in
G` has its own specific cross section parameters. This makes the difficult
non-linear continuous energy dependence much more manageable as group wise cross
section data can be precomputed and fed into a simulation as input data. The
computation of multigroup cross section data is not a trivial task and can
introduce errors in the simulation. However, this is an active field of research
common to all multigroup methods, and there are numerous generation methods
available that are capable of reducing the biases introduced by the multigroup
approximation. Commonly used methods include the subgroup self-shielding method
and use of fast (unconverged) Monte Carlo simulations to produce cross section
estimates. It is important to note that Monte Carlo methods are capable of
treating the energy variable of the neutron continuously, meaning that they do
not need to make this approximation and are therefore not subject to any
multigroup errors.

Following the multigroup discretization, another assumption made is that a large
and complex problem can be broken up into small constant cross section regions,
and that these regions have group dependent, flat, isotropic sources (fission
and scattering), :math:`Q_g`. Anisotropic as well as higher order sources are
also possible with MOC-based methods but are not used yet in OpenMC for
simplicity. With these key assumptions, the multigroup MOC form of the neutron
transport equation can be written as in Equation :eq:`moc_final`.

.. math::
    :label: moc_final

    \psi_g(s, \mathbf{\Omega}) = \psi_g(\mathbf{r_0}, \mathbf{\Omega}) e^{-\int_0^s ds^\prime \Sigma_{t_g}(s^\prime)} + \int_0^s ds^{\prime\prime} Q_g(s^{\prime\prime},\mathbf{\Omega}) e^{-\int_{s^{\prime\prime}}^s ds^\prime \Sigma_{t_g}(s^\prime)}

The CSG definition of the system is used to create spatially defined source
regions (each region being denoted as :math:`i`). These neutron source regions
are often approximated as being constant
(flat) in source intensity but can also be defined using a higher order source
(linear, quadratic, etc.) that allows for fewer source regions to be required to
achieve a specified solution fidelity. In OpenMC, the approximation of a
spatially constant isotropic fission and scattering source :math:`Q_{i,g}` in
cell :math:`i` leads
to simple exponential attenuation along an individual characteristic of length
:math:`s` given by Equation :eq:`fsr_attenuation`.

.. math::
    :label: fsr_attenuation

    \psi_g(s) = \psi_g(0) e^{-\Sigma_{t,i,g} s} + \frac{Q_{i,g}}{\Sigma_{t,i,g}} \left( 1 - e^{-\Sigma_{t,i,g} s} \right)

For convenience, we can also write this equation in terms of the incoming and
outgoing angular flux (:math:`\psi_g^{in}` and :math:`\psi_g^{out}`), and
consider a specific tracklength for a particular ray :math:`r` crossing cell
:math:`i` as :math:`\ell_r`, as in:

.. math::
    :label: fsr_attenuation_in_out

    \psi_g^{out} = \psi_g^{in} e^{-\Sigma_{t,i,g} \ell_r} + \frac{Q_{i,g}}{\Sigma_{t,i,g}} \left( 1 - e^{-\Sigma_{t,i,g} \ell_r} \right) .

We can then define the average angular flux of a single ray passing through the
cell as:

.. math::
    :label: average

    \overline{\psi}_{r,i,g} = \frac{1}{\ell_r} \int_0^{\ell_r} \psi_{g}(s)ds .

We can then substitute in Equation :eq:`fsr_attenuation` and solve, resulting
in:

.. math::
    :label: average_solved

    \overline{\psi}_{r,i,g} = \frac{Q_{i,g}}{\Sigma_{t,i,g}} - \frac{\psi_{r,g}^{out} - \psi_{r,g}^{in}}{\ell_r \Sigma_{t,i,g}} .

By rearranging Equation :eq:`fsr_attenuation_in_out`, we can then define
:math:`\Delta \psi_{r,g}` as the change in angular flux for ray :math:`r`
passing through region :math:`i` as:

.. math::
    :label: delta_psi

    \Delta \psi_{r,g} = \psi_{r,g}^{in} - \psi_{r,g}^{out} = \left(\psi_{r,g}^{in} - \frac{Q_{i,g}}{\Sigma_{t,i,g}} \right) \left( 1 - e^{-\Sigma_{t,i,g} \ell_r} \right) .

Equation :eq:`delta_psi` is a useful expression as it is easily computed with
the known inputs for a ray crossing through the region.

By substituting :eq:`delta_psi` into :eq:`average_solved`, we can arrive at a
final expression for the average angular flux for a ray crossing a region as:

.. math::
    :label: average_psi_final

    \overline{\psi}_{r,i,g} = \frac{Q_{i,g}}{\Sigma_{t,i,g}} + \frac{\Delta \psi_{r,g}}{\ell_r \Sigma_{t,i,g}}

~~~~~~~~~~~
Random Rays
~~~~~~~~~~~

In the previous subsection, the governing characteristic equation along a 1D
line through the system was written, such that an analytical solution for the
ODE can be computed. If enough characteristic tracks (ODEs) are solved, then the
behavior of the governing PDE can be numerically approximated. In traditional
deterministic MOC, the selection of tracks is chosen deterministically, where
azimuthal and polar quadratures are defined along with even track spacing in
three dimensions. This is the point at which random ray diverges from
deterministic MOC numerically. In the random ray method, rays are randomly
sampled from a uniform distribution in space and angle and tracked along a
predefined distance through the geometry before terminating. **Importantly,
different rays are sampled each power iteration, leading to a fully stochastic
convergence process.** This results in a need to utilize both inactive and
active batches as in the Monte Carlo method.

While Monte Carlo implicitly converges the scattering source fully within each
iteration, random ray (and MOC) solvers are not typically written to fully
converge the scattering source within a single iteration. Rather, both the
fission and scattering sources are updated each power iteration, thus requiring
enough outer iterations to reach a stationary distribution in both the fission
source and scattering source. So, even in a low dominance ratio problem like a
2D pincell, several hundred inactive batches may still be required with random
ray to allow the scattering source to fully develop, as neutrons undergoing
hundreds of scatters may constitute a non-trivial contribution to the fission
source. We note that use of a two-level second iteration scheme is sometimes
used by some MOC or random ray solvers so as to fully converge the scattering
source with many inner iterations before updating the fission source in the
outer iteration. It is typically more efficient to use the single level
iteration scheme, as there is little reason to spend so much work converging the
scattering source if the fission source is not yet converged.

Overall, the difference in how random ray and Monte Carlo converge the
scattering source means that in practice, random ray typically requires more
inactive iterations than are required in Monte Carlo. While a Monte Carlo
simulation may need 100 inactive iterations to reach a stationary source
distribution for many problems, a random ray solve will likely require 1,000
iterations or more. Source convergence metrics (e.g., Shannon entropy) are thus
recommended when performing random ray simulations to ascertain when the source
has fully developed.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Converting Angular Flux to Scalar Flux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thus far in our derivation, we have been able to write analytical equations that
solve for the change in angular flux of a ray crossing a flat source region
(Equation :eq:`delta_psi`) as well as the ray's average angular flux through
that region (Equation :eq:`average_psi_final`). To determine the source for the
next power iteration, we need to assemble our estimates of angular fluxes from
all the sampled rays into scalar fluxes within each FSR.

We can define the scalar flux in region :math:`i` as:

.. math::
    :label: integral

    \phi_i = \frac{\int_{V_i} \int_{4\pi} \psi(r, \Omega) d\Omega d\mathbf{r}}{\int_{V_i} d\mathbf{r}} .

The integral in the numerator:

.. math::
    :label: numerator

    \int_{V_i} \int_{4\pi} \psi(r, \Omega) d\Omega d\mathbf{r} .

is not known analytically, but with random ray, we are going the numerically
approximate it by discretizing over a finite number of tracks (with a finite
number of locations and angles) crossing the domain. We can then use the
characteristic method to determine the total angular flux along that line.

Conceptually, this can be thought of as taking a volume-weighted sum of angular
fluxes for all :math:`N_i` rays that happen to pass through cell :math:`i` that
iteration. When written in discretized form (with the discretization happening
in terms of individual ray segments :math:`r` that pass through region
:math:`i`), we arrive at:

.. math::
    :label: discretized

    \phi_{i,g} = \frac{\int_{V_i} \int_{4\pi} \psi(r, \Omega) d\Omega d\mathbf{r}}{\int_{V_i} d\mathbf{r}} = \overline{\overline{\psi}}_{i,g} \approx \frac{\sum\limits_{r=1}^{N_i} \ell_r w_r \overline{\psi}_{r,i,g}}{\sum\limits_{r=1}^{N_i} \ell_r w_r} .

Here we introduce the term :math:`w_r`, which represents the "weight" of the ray
(its 2D area), such that the volume that a ray is responsible for can be
determined by multiplying its length :math:`\ell` by its weight :math:`w`. As
the scalar flux vector is a shape function only, we are actually free to
multiply all ray weights :math:`w` by any constant such that the overall shape
is still maintained, even if the magnitude of the shape function changes. Thus,
we can simply set :math:`w_r` to be unity for all rays, such that:

.. math::
    :label: weights

    \text{Volume of cell } i = V_i \approx \sum\limits_{r=1}^{N_i} \ell_r w_r = \sum\limits_{r=1}^{N_i} \ell_r .

We can then rewrite our discretized equation as:

.. math::
    :label: discretized_2

    \phi_{i,g} \approx \frac{\sum\limits_{r=1}^{N_i} \ell_r w_r \overline{\psi}_{r,i,g}}{\sum\limits_{r=1}^{N_i} \ell_r w_r} = \frac{\sum\limits_{r=1}^{N_i} \ell_r \overline{\psi}_{r,i,g}}{\sum\limits_{r=1}^{N_i} \ell_r} .

Thus, the scalar flux can be inferred if we know the volume weighted sum of the
average angular fluxes that pass through the cell. Substituting
:eq:`average_psi_final` into :eq:`discretized_2`, we arrive at:

.. math::
    :label: scalar_full

    \phi_{i,g} = \frac{\int_{V_i} \int_{4\pi} \psi(r, \Omega) d\Omega d\mathbf{r}}{\int_{V_i} d\mathbf{r}} = \overline{\overline{\psi}}_{i,g} = \frac{\sum\limits_{r=1}^{N_i} \ell_r \overline{\psi}_{r,i,g}}{\sum\limits_{r=1}^{N_i} \ell_r} = \frac{\sum\limits_{r=1}^{N_i} \ell_r \frac{Q_{i,g}}{\Sigma_{t,i,g}} + \frac{\Delta \psi_{r,g}}{\ell_r \Sigma_{t,i,g}}}{\sum\limits_{r=1}^{N_i} \ell_r},

which when partially simplified becomes:

.. math::
    :label: scalar_four_vols

    \phi =  \frac{Q_{i,g} \sum\limits_{r=1}^{N_i} \ell_r}{\Sigma_{t,i,g} \sum\limits_{r=1}^{N_i} \ell_r} + \frac{\sum\limits_{r=1}^{N_i} \ell_r \frac{\Delta \psi_i}{\ell_r}}{\Sigma_{t,i,g} \sum\limits_{r=1}^{N_i} \ell_r} .

Note that there are now four (seemingly identical) volume terms in this equation.

~~~~~~~~~~~~~~
Volume Dilemma
~~~~~~~~~~~~~~

At first glance, Equation :eq:`scalar_four_vols` appears ripe for cancellation
of terms. Mathematically, such cancellation allows us to arrive at the following
"naive" estimator for the scalar flux:

.. math::
    :label: phi_naive

    \phi_{i,g}^{naive} = \frac{Q_{i,g} }{\Sigma_{t,i,g}} + \frac{\sum\limits_{r=1}^{N_i} \Delta \psi_{r,g}}{\Sigma_{t,i,g} \sum\limits_{r=1}^{N_i} \ell_r} .

This derivation appears mathematically sound at first glance but unfortunately
raises a serious issue as discussed in more depth by `Tramm et al.
<Tramm-2020>`_ and `Cosgrove and Tramm <Cosgrove-2023>`_. Namely, the second
term:

.. math::
    :label: ratio_estimator

     \frac{\sum\limits_{r=1}^{N_i} \Delta \psi_{r,g}}{\Sigma_{t,i,g} \sum\limits_{r=1}^{N_i} \ell_r}

features stochastic variables (the sums over random ray lengths and angular
fluxes) in both the numerator and denominator, making it a stochastic ratio
estimator, which is inherently biased. In practice, usage of the naive estimator
does result in a biased, but "consistent"  estimator (i.e., it is biased, but
the bias tends towards zero as the sample size increases). Experimentally, the
right answer can be obtained with this estimator, though a very fine ray density
is required to eliminate the bias.

How might we solve the biased ratio estimator problem? While there is no obvious
way to alter the numerator term (which arises from the characteristic
integration approach itself), there is potentially more flexibility in how we
treat the stochastic term in the denominator, :math:`\sum\limits_{r=1}^{N_i}
\ell_r` . From Equation :eq:`weights` we know that this term can be directly
inferred from the volume of the problem, which does not actually change between
iterations. Thus, an alternative treatment for this "volume" term in the
denominator is to replace the actual stochastically sampled total track length
with the expected value of the total track length. For instance, if the true
volume of the FSR is known (as is the total volume of the full simulation domain
and the total tracklength used for integration that iteration), then we know the
true expected value of the tracklength in that FSR. That is, if a FSR accounts
for 2% of the overall volume of a simulation domain, then we know that the
expected value of tracklength in that FSR will be 2% of the total tracklength
for all rays that iteration. This is a key insight, as it allows us to the
replace the actual tracklength that was accumulated inside that FSR each
iteration with the expected value.

If we know the analytical volumes, then those can be used to directly compute
the expected value of the tracklength in each cell. However, as the analytical
volumes are not typically known in OpenMC due to the usage of user-defined
constructive solid geometry, we need to source this quantity from elsewhere. An
obvious choice is to simply accumulate the total tracklength through each FSR
across all iterations (batches) and to use that sum to compute the expected
average length per iteration, as:

.. math::
    :label: sim_estimator

       \sum\limits^{}_{i} \ell_i \approx \frac{\sum\limits^{B}_{b}\sum\limits^{N_i}_{r} \ell_{b,r} }{B}

where :math:`b` is a single batch in :math:`B` total batches simulated so far.

In this manner, the expected value of the tracklength will become more refined
as iterations continue, until after many iterations the variance of the
denominator term becomes trivial compared to the numerator term, essentially
eliminating the presence of the stochastic ratio estimator. A "simulation
averaged" estimator is therefore:

.. math::
    :label: phi_sim

    \phi_{i,g}^{simulation} = \frac{Q_{i,g} }{\Sigma_{t,i,g}} + \frac{\sum\limits_{r=1}^{N_i} \Delta \psi_{r,g}}{\Sigma_{t,i,g} \frac{\sum\limits^{B}_{b}\sum\limits^{N_i}_{r} \ell_{b,r} }{B}}

In practical terms, the "simulation averaged" estimator is virtually
indistinguishable numerically from use of the true analytical volume to estimate
this term. Note also that the term "simulation averaged" refers only to the
volume/length treatment, the scalar flux estimate itself is computed fully again
each iteration.

There are some drawbacks to this method. Recall, this denominator volume term
originally stemmed from taking a volume weighted integral of the angular flux,
in which case the denominator served as a normalization term for the numerator
integral in Equation :eq:`integral`. Essentially, we have now used a different
term for the volume in the numerator as compared to the normalizing volume in
the denominator. The inevitable mismatch (due to noise) between these two
quantities results in a significant increase in variance. Notably, the same
problem occurs if using a tracklength estimate based on the analytical volume,
as again the numerator integral and the normalizing denominator integral no
longer match on a per-iteration basis.

In practice, the simulation averaged method does completely remove the bias,
though at the cost of a notable increase in variance. Empirical testing reveals
that on most problems, the simulation averaged estimator does win out overall in
numerical performance, as a much coarser quadrature can be used resulting in
faster runtimes overall. Thus, OpenMC uses the simulation averaged estimator in
its random ray mode.

~~~~~~~~~~~~~~~
Power Iteration
~~~~~~~~~~~~~~~

Given a starting source term, we now have a way of computing an estimate of the
scalar flux in each cell by way of transporting rays randomly through the
domain, recording the change in angular flux for the rays into each cell as they
make their traversals, and summing these contributions up as in Equation
:eq:`phi_sim`. How then do we turn this into an iterative process such that we
improve the estimate of the source and scalar flux over many iterations, given
that our initial starting source will just be a guess?

The source :math:`Q^{n}` for iteration :math:`n` can be inferred
from the scalar flux from the previous iteration :math:`n-1` as:

.. math::
    :label: source_update

    Q^{n}(i, g) = \frac{\chi}{k^{n-1}_{eff}} \nu \Sigma_f(i, g) \phi^{n-1}(g) + \sum\limits^{G}_{g'} \Sigma_{s}(i,g,g') \phi^{n-1}(g')

where :math:`Q^{n}(i, g)` is the total source (fission + scattering) in region
:math:`i` and energy group :math:`g`. Notably, the in-scattering source in group
:math:`g` must be computed by summing over the contributions from all groups
:math:`g' \in G`.

In a similar manner, the eigenvalue for iteration :math:`n` can be computed as:

.. math::
    :label: eigenvalue_update

    k^{n}_{eff} = k^{n-1}_{eff} \frac{F^n}{F^{n-1}},

where the total spatial- and energy-integrated fission rate :math:`F^n` in
iteration :math:`n` can be computed as:

.. math::
    :label: fission_source

    F^n = \sum\limits^{M}_{i} \left( V_i \sum\limits^{G}_{g} \nu \Sigma_f(i, g) \phi^{n}(g) \right)

where :math:`M` is the total number of FSRs in the simulation. Similarly, the
total spatial- and energy-integrated fission rate :math:`F^{n-1}` in iteration
:math:`n-1` can be computed as:

.. math::
    :label: fission_source_prev

    F^{n-1} = \sum\limits^{M}_{i} \left( V_i \sum\limits^{G}_{g} \nu \Sigma_f(i, g) \phi^{n-1}(g) \right)

Notably, the volume term :math:`V_i` appears in the eigenvalue update equation.
The same logic applies to the treatment of this term as was discussed earlier.
In OpenMC, we use the "simulation averaged" volume derived from summing over all
ray tracklength contributions to a FSR over all iterations and dividing by the
total integration tracklength to date. Thus, Equation :eq:`fission_source`
becomes:

.. math::
    :label: fission_source_volumed

    F^n = \sum\limits^{M}_{i} \left( \frac{\sum\limits^{B}_{b}\sum\limits^{N_i}_{r} \ell_{b,r} }{B} \sum\limits^{G}_{g} \nu \Sigma_f(i, g) \phi^{n}(g) \right)

and a similar substitution can be made to update Equation
:eq:`fission_source_prev` . In OpenMC, the most up-to-date version of the volume
estimate is used, such that the total fission source from the previous iteration
(:math:`n-1`) is also recomputed each iteration.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ray Starting Conditions and Inactive Length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another key area of divergence between deterministic MOC and random ray is the
starting conditions for rays. In deterministic MOC, the angular flux spectrum
for rays are stored at any reflective or periodic boundaries so as to provide a
starting condition for the next iteration. As there are many tracks, storage of
angular fluxes can become costly in terms of memory consumption unless there are
only vacuum boundaries present.

In random ray, as the starting locations of rays are sampled anew each
iteration, the initial angular flux spectrum for the ray is unknown. While a
guess can be made by taking the isotropic source from the FSR the ray was
sampled in, direct usage of this quantity would result in significant bias and
error being imparted on the simulation.

Thus, an `on-the-fly approximation method <Tramm-2017a>`_ was developed (known
as the "dead zone"), where the first several mean free paths of a ray are
considered to be "inactive" or "read only". In this sense, the angular flux is
solved for using the MOC equation, but the ray does not "tally" any scalar flux
back to the FSRs that it travels through. After several mean free paths have
been traversed, the ray's angular flux spectrum typically becomes dominated by
the accumulated source terms from the cells it has traveled through, while the
(incorrect) starting conditions have been attenuated away. In the animation in
the :ref:`introductory section on this page <methods_random_ray_intro>`, the
yellow portion of the ray lengths is the dead zone. As can be seen in this
animation, the tallied :math:`\sum\limits_{r=1}^{N_i} \Delta \psi_{r,g}` term
that is plotted is not affected by the ray when the ray is within its inactive
length. Only when the ray enters its active mode does the ray contribute to the
:math:`\sum\limits_{r=1}^{N_i} \Delta \psi_{r,g}` sum for the iteration.

~~~~~~~~~~~~~~~~~~~~~
Ray Ending Conditions
~~~~~~~~~~~~~~~~~~~~~

To ensure that a uniform density of rays is integrated in space and angle
throughout the simulation domain, after exiting the initial inactive "dead zone"
portion of the ray, the rays are run for a user-specified distance. Typically, a
choice of at least several times the length of the inactive "dead zone" is made
so as to amortize the cost of the dead zone. For example, if a dead zone of 30
cm is selected, then an active length of 300 cm might be selected so that the
cost of the dead zone is at most 10% of the overall runtime.

--------------------
Simplified Algorithm
--------------------

A simplified set of functions that execute a single random ray power iteration
are given below. Not all global variables are defined in this illustrative
example, but the high level components of the algorithm are shown. A number of
significant simplifications are made for clarity---for example, no inactive
"dead zone" length is shown, geometry operations are abstracted, no parallelism
(or thread safety) is expressed, a naive exponential treatment is used, and rays
are not halted at their exact termination distances, among other subtleties.
Nonetheless, the below algorithms may be useful for gaining intuition on the
basic components of the random ray process. Rather than expressing the algorithm
in abstract pseudocode, C++ is used to make the control flow easier to
understand.

The first block below shows the logic for a single power iteration (batch):

.. code-block:: C++

    double power_iteration(double k_eff) {

        // Update source term (scattering + fission)
        update_neutron_source(k_eff);

        // Reset scalar fluxes to zero
        fill<float>(global::scalar_flux_new, 0.0f);

        // Transport sweep over all random rays for the iteration
        for (int i = 0; i < nrays; i++) {
            RandomRay ray;
            initialize_ray(ray);
            transport_single_ray(ray);
        }

        // Normalize scalar flux and update volumes
        normalize_scalar_flux_and_volumes();

        // Add source to scalar flux, compute number of FSR hits
        add_source_to_scalar_flux();

        // Compute k-eff using updated scalar flux
        k_eff = compute_k_eff(k_eff);

        // Set phi_old = phi_new
        global::scalar_flux_old.swap(global::scalar_flux_new);

        return k_eff;
    }

The second function shows the logic for transporting a single ray within the
transport loop:

.. code-block:: C++

    void transport_single_ray(RandomRay& ray) {

        // Reset distance to zero
        double distance = 0.0;

        // Continue transport of ray until active length is reached
        while (distance < user_setting::active_length) {
            // Ray trace to find distance to next surface (i.e., segment length)
            double s = distance_to_nearest_boundary(ray);

            // Attenuate flux (and accumulate source/attenuate) on segment
            attenuate_flux(ray, s);

            // Advance particle to next surface
            ray.location = ray.location + s * ray.direction;

            // Move ray across the surface
            cross_surface(ray);

            // Add segment length "s" to total distance traveled
            distance += s;
        }
    }

The final function below shows the logic for solving for the characteristic MOC
equation (and accumulating the scalar flux contribution of the ray into the
scalar flux value for the FSR).

.. code-block:: C++

    void attenuate_flux(RandomRay& ray, double s) {

        // Determine which flat source region (FSR) the ray is currently in
        int fsr = get_fsr_id(ray.location);

        // Determine material type
        int material = get_material_type(fsr);

        // MOC incoming flux attenuation + source contribution/attenuation equation
         for (int e = 0; e < global::n_energy_groups; e++) {
            float sigma_t = global::macro_xs[material].total;
            float tau = sigma_t * s;
            float delta_psi = (ray.angular_flux[e] - global::source[fsr][e] / sigma_t) * (1 - exp(-tau));
            ray.angular_flux_[e] -= delta_psi;
            global::scalar_flux_new[fsr][e] += delta_psi;
        }

        // Record total tracklength in this FSR (to compute volume)
        global::volume[fsr] += s;
    }

------------------------
How are Tallies Handled?
------------------------

Most tallies, filters, and scores that you would expect to work with a
multigroup solver like random ray should work. For example, you can define 3D
mesh tallies with energy filters and flux, fission, and nu-fission scores, etc.
There are some restrictions though. For starters, it is assumed that all filter
mesh boundaries will conform to physical surface boundaries (or lattice
boundaries) in the simulation geometry. It is acceptable for multiple cells
(FSRs) to be contained within a filter mesh cell (e.g., pincell-level or
assembly-level tallies should work), but it is currently left as undefined
behavior if a single simulation cell is able to score to multiple filter mesh
cells. In the future, the capability to fully support mesh tallies may be added
to OpenMC, but for now this restriction needs to be respected.

---------------------------
Fundamental Sources of Bias
---------------------------

Compared to continuous energy Monte Carlo simulations, the known sources of bias
in random ray particle transport are:

    - **Multigroup Energy Discretization:** The multigroup treatment of flux and
      cross sections incurs a significant bias, as a reaction rate (:math:`R_g =
      V \phi_g \Sigma_g`) for an energy group :math:`g` can only be conserved
      for a given choice of multigroup cross section :math:`\Sigma_g` if the
      flux (:math:`\phi_g`) is known a priori. If the flux was already known,
      then there would be no point to the simulation, resulting in a fundamental
      need for approximating this quantity. There are numerous methods for
      generating relatively accurate multigroup cross section libraries that can
      each be applied to a narrow design area reliably, although there are
      always limitations and/or complexities that arise with a multigroup energy
      treatment. This is by far the most significant source of simulation bias
      between Monte Carlo and random ray for most problems. While the other
      areas typically have solutions that are highly effective at mitigating
      bias, error stemming from multigroup energy discretization is much harder
      to remedy.
    - **Flat Source Approximation:**. In OpenMC, a "flat" (0th order) source
      approximation is made, wherein the scattering and fission sources within a
      cell are assumed to be spatially uniform. As the source in reality is a
      continuous function, this leads to bias, although the bias can be reduced
      to acceptable levels if the flat source regions are sufficiently small.
      The bias can also be mitigated by assuming a higher-order source (e.g.,
      linear or quadratic), although OpenMC does not yet have this capability.
      In practical terms, this source of bias can become very large if cells are
      large (with dimensions beyond that of a typical particle mean free path),
      but the subdivision of cells can often reduce this bias to trivial levels.
    - **Anisotropic Source Approximation:** In OpenMC, the source is not only
      assumed to be flat but also isotropic, leading to bias. It is possible for
      MOC (and likely random ray) to treat anisotropy explicitly, but this is
      not currently supported in OpenMC. This source of bias is not significant
      for some problems, but becomes more problematic for others. Even in the
      absence of explicit treatment of anistropy, use of transport-corrected
      multigroup cross sections can often mitigate this bias, particularly for
      light water reactor simulation problems.
    - **Angular Flux Initial Conditions:** Each time a ray is sampled, its
      starting angular flux is unknown, so a guess must be made (typically the
      source term for the cell it starts in). Usage of an adequate inactive ray
      length (dead zone) mitigates this error. As the starting guess is
      attenuated at a rate of :math:`\exp(-\Sigma_t \ell)`, this bias can driven
      below machine precision in a low cost manner on many problems.

.. _Tramm-2017a: https://doi.org/10.1016/j.jcp.2017.04.038
.. _Tramm-2017b: https://doi.org/10.1016/j.anucene.2017.10.015
.. _Tramm-2018: https://dspace.mit.edu/handle/1721.1/119038
.. _Tramm-2020: https://doi.org/10.1051/EPJCONF/202124703021
.. _Cosgrove-2023: https://doi.org/10.1080/00295639.2023.2270618

.. only:: html

   .. rubric:: References

.. [Askew-1972] Askew, “A Characteristics Formulation of the Neutron Transport
    Equation in Complicated Geometries.” Technical Report AAEW-M 1108, UK Atomic
    Energy Establishment (1972).
