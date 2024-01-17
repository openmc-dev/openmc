.. _methods_random_ray:

===============
Random Ray
===============

-----------------------------------
What is Random Ray?
-----------------------------------

Random ray [`Tramm 2017a`_] is a stochastic transport method, closely related to the deterministic Method of Characteristics (MOC) [Askew]_. Rather than each ray representing a single neutron as in Monte Carlo, it represents a characteristic line through the reactor upon which the transport equation can be written as an ordinary differential equation that can be solved analytically (although with discretization required in energy space, making it a multigroup method). The behavior of the governing transport equation can be approximated by solving along many characteristic tracks (rays) through the reactor. Unlike particles in Monte Carlo, rays in random ray or MOC are not affected by the material characteristics of the simulated problem -- rays are selected so as to explore the full simulation problem with a statistically equal distribution in space and angle.

-----------------------------------------------
Why is a Random Ray Solver Included in OpenMC?
-----------------------------------------------

There are a few good reasons:

* The random ray solver complements the capabilities of MC nicely. One area that MC struggles with is maintaining accuracy in regions of low physical particle flux. Random ray, on the other hand, has approximately even variance throughout the entire global simulation domain, such that areas with low neutron flux are no less well known that areas of high neutron flux. Absent weight windows in MC, random ray can be several orders of magnitude faster than multigroup Monte Carlo in classes of problems where areas with low physical neutron flux need to be resolved. While MC uncertainty can be greatly improved with variance reduction techniques, they add some user complexity, and weight windows can often be expensive to generate via MC transport alone (e.g., MAGIC). While not yet implemented in this PR, I also plan to add capability so that the random ray solver can be used to generate fast and high quality weight windows via FW-CADIS that can then be used to accelerate convergence in MC. Early work in the field has shown significant speedup in weight window generation and weight window quality with random ray and FW-CADIS as compared to MAGIC.

* In practical implementation terms, random ray is mechanically very similar to how Monte Carlo works, in terms of the process of ray tracing on CSG geometry and handling stochastic convergence, etc. In the original 1972 paper by Askew that introduces MOC (which random ray is a variant of), he stated:

    .. epigraph:: 
    
        "One of the features of the method proposed [MoC] is that ... the tracking process needed to perform this operation is common to the proposed method ... and to Monte Carlo methods. Thus a single tracking routine capable of recognizing a geometric arrangement could be utilized to service all types of solution, choice being made depending which was more appropriate to the problem size and required accuracy."

        -- Askew

    This prediction holds up -- the additional requirements needed in OpenMC to handle random ray transport turned out to be fairly small, only in the range of 800 lines of code for the first implementation I did. Additionally, most of the solver was able to be added in a neatly "silo'd" fashion that doesn't complicate the MC portion of OpenMC. The current PR has increased in size to around 2400 lines of code, but much of that is due to the new example, regression test, and other boilerplate interface code etc. Thus, for relatively few additional lines of code, users will have capabilities to more efficiently handle certain types of simulation problems that are more computationally challenging to solve with MC, and will definitely have a faster transport method for generating weight windows.

* It amortizes the code complexity in OpenMC for representing multigroup cross sections. There is a significant amount of interface code, documentation, and complexity in allowing OpenMC to generate and use multigroup XS data in its MGMC mode. Random ray allows the same multigroup data to be used, making full reuse of these existing capabilities.

------------------------------------------
Random Ray Numerical Derivation
------------------------------------------

The derivation of Random Ray is discussed extensively in [`Tramm 2017a`_], but we will reproduce portions of this dissertation varbatim in this section for convenience.

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Method of Characteristics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Boltzmann neutron transport equation is a partial differential equation (PDE) that describes the angular flux within a system. It is a balance equation, with the streaming and absorption terms typically appearing on the left hand side, which are balanced by the scattering source and fission source terms on the right hand side. 

.. math::
    :label: transport

    \vec{\Omega} \cdot \vec{\nabla} \psi(\vec{r},\vec{\Omega},E) + \Sigma_t(\vec{r},E) \psi(\vec{r},\vec{\Omega},E) = \int_0^\infty d E^\prime \int_{4\pi} d \Omega^{\prime} \Sigma_s(\vec{r},\vec{\Omega}^\prime \rightarrow \vec{\Omega}, E^\prime \rightarrow E) \psi(\vec{r},\vec{\Omega}^\prime, E^\prime) + \frac{\chi(\vec{r}, E)}{4\pi k_{eff}} \int_0^\infty dE^\prime \nu \Sigma_f(\vec{r},E^\prime) \int_{4\pi}d \Omega^\prime \psi(\vec{r},\vec{\Omega}^\prime,E^\prime)

Equation :eq:`transport` defines the neutron transport equation where :math:`\psi`` is the angular neutron flux. This parameter represents the total distance travelled by all neutrons in a particular direction inside of a control volume per second, and is often given in units of :math:`1/(cm^{2} s)`. The angular direction unit vector, :math:`\vec{\Omega}`, represents the direction of travel for the neutron. The spatial position vector, :math:`\vec{r}`,  represents the location within the simulation. The neutron energy, :math:`E`, or speed in continuous space, is often given in units of electron volts. The total macroscopic neutron cross section is :math:`\Sigma_t`` . This value represents the total probability of interaction between a neutron travelling at a certain speed (i.e., neutron energy :math:`E`) and a target nucleus (i.e., the material through which the neutron is travelling) per unit path length, typically given in units of :math:`1/cm`. Macroscopic cross section data is a combination of empirical data and quantum mechanical modelling employed in order to generate an evaluation represented either in pointwise form or resonance parameters for each target isotope of interest in a material, as well as the density of the material, and is provided as input to a simulation. The scattering neutron cross section, :math:`\Sigma_s`, is similar to the total cross section but only measures scattering interactions between the neutron and the target nucleus, and includes extra dependencies on the change in angle and energy the neutron experiences as a result of the interaction. Several additional reactions like (n,2n) and (n,3n) are included in the scattering transfer cross section. The fission neutron cross section, :math:`\Sigma_f`, is also similar to the total cross section but only measures the fission interaction between a neutron and a target nucleus. The energy spectrum for neutrons born from fission, :math:`\chi`, represents a known distribution of outgoing neutron energies based on the material that fissioned, which is taken as input data to a computation. The average number of neutrons born per fission is :math:`\nu`. The eigenvalue of the equation, :math:`k_{eff}`, represents the effective neutron multiplication factor. If the right hand side of Equation :eq:`transport` is condensed into a single term, represented by the total neutron source term :math:`Q(\vec{r}, \vec{\Omega},E)`, the form given in Equation :eq:`transport_simple` is reached."

.. math::
    :label: transport_simple

    \overbrace{\vec{\Omega} \cdot \vec{\nabla} \psi(\vec{r},\vec{\Omega},E)}^{\text{streaming term}} + \overbrace{\Sigma_t(\vec{r},E) \psi(\vec{r},\vec{\Omega},E)}^{\text{absorption term}} = \overbrace{Q(\vec{r}, \vec{\Omega},E)}^{\text{total neutron source term}}

Fundamentally, MOC works by solving Equation :eq:`transport_simple` along a single characteristic line, thus altering the full spatial and angular scope of the transport equation into something that holds true only for a particular linear path (or track) through the reactor. These tracks are linear for neutrons as they are neutral particles and are therefore not subject to field effects. To accomplish this, we parameterize :math:`\vec{r}` with respect to some reference location :math:`\vec{r}_0` such that :math:`\vec{r} = \vec{r}_0 + s\vec{\Omega}`. In this manner, Equation :eq:`transport_simple` can be rewritten for a specific segment length :math:`s` at a specific angle :math:`\vec{\Omega}` through a constant cross section region of the reactor geometry as in Equation :eq:`char_long`.

.. math::
    :label: char_long

    \vec{\Omega} \cdot \vec{\nabla} \psi(\vec{r}_0 + s\vec{\Omega},\vec{\Omega},E) + \Sigma_t(\vec{r}_0 + s\vec{\Omega},E) \psi(\vec{r}_0 + s\vec{\Omega},\vec{\Omega},E) = Q(\vec{r}_0 + s\vec{\Omega}, \vec{\Omega},E)

As this equation holds along a 1 dimensional path, we can assume the dependence of :math:`s` on :math:`\vec{r}_0` and :math:`\vec{\Omega}`` such that :math:`\vec{r}_0 + s\vec{\Omega}` simplifies to :math:`s`. When the differential operator is also applied to the angular flux :math:`\psi`, we arrive at the characteristic form of the Boltzmann Neutron Transport Equation given in Equation :eq:`char`.

.. math::
    :label: char

    \frac{d}{ds} \psi(s,\vec{\Omega},E) + \Sigma_t(s,E) \psi(s,\vec{\Omega},E) = Q(s, \vec{\Omega},E)

An analytical solution to this characteristic equation can be achieved with the use of an integrating factor:

.. math::
    :label: int_factor

    e^{ \int_0^s ds' \Sigma^T (s', E)}

to arrive at the final form of the characteristic equation shown in Equation :eq:`full_char`.

.. math::
    :label: full_char

    \psi(s,\vec{\Omega},E) = \psi(\vec{r}_0,\vec{\Omega},E) e^{-\int_0^s ds^\prime \Sigma_t(s^\prime,E)} + \int_0^s ds^{\prime\prime} Q(s^{\prime\prime},\vec{\Omega}, E) e^{-\int_{s^{\prime\prime}}^s ds^\prime \Sigma_t(s^\prime,E)}

With this characteristic form of the transport equation, we now have an analytical solution along a linear path through any constant cross section region of a reactor, with only the continuous energy dependence remaining as an issue to be addressed.  Similar to many other solution approaches to the Boltzmann neutron transport equation, the MOC approach also uses a "multi-group" approximation in order to discretize the continuous energy spectrum of neutrons travelling through the reactor into fixed set of energy groups :math:`G`, where each group :math:`g \in G` has its own specific cross section parameters. This makes the difficult non-linear continuous energy dependence much more manageable as group wise cross section data can be precomputed and fed into a simulation as input data. The computation of multi-group cross section data is not a trivial task and can introduce errors in the simulation. However, this is an active field of research common to all multi-group methods, and there are numerous generation methods available that are capable of minimizing the biases introduced by the multi-group approximation. Commonly used methods include the subgroup self-shielding method and use of smaller Monte Carlo simulations to produce cross section data. It is important to note that Monte Carlo methods are capable of treating the energy variable of the neutron continuously, meaning that they do not need to make this approximation and are therefore not subject to any multi-group errors.

Following the multi-group assumption, another assumption made is that a large and complex problem can be broken up into small constant cross section regions, and that these regions have group dependent, flat, isotropic sources (fission + scattering), :math:`Q_g`. Anisotropic as well as higher order sources are also possible with MOC-based methods, but are not used yet in OpenMC for simplicity. With these key assumptions, the multi-group MOC form of the neutron transport equation can be written as in Equation :eq:`moc_final`.

.. math::
    :label: moc_final

    \psi_g(s, \vec{\Omega}) = \psi_g(\vec{r_0}, \vec{\Omega}) e^{-\int_0^s ds^\prime \Sigma_{t_g}(s^\prime)} + \int_0^s ds^{\prime\prime} Q_g(s^{\prime\prime},\vec{\Omega}) e^{-\int_{s^{\prime\prime}}^s ds^\prime \Sigma_{t_g}(s^\prime)}

The constructive solid geometry (CSG) definition of the reactor is used to create spatially defined source regions. These neutron source regions are often approximated as being constant (flat) in intensity of source, but can also be defined using a higher order source (linear, quadratic, etc.) that allows for fewer source regions to be required to achieve a specified solution fidelity. In OpenMC, the normal approximation of a spatially constant isotropic fission and scattering source :math:`q_0` leads to simple exponential attenuation along an individual characteristic of length :math:`s` given by Equation :eq:`fsr_attenuation`.

.. math::
    :label: fsr_attenuation

    \psi_g(s) = \psi_g(0) e^{-\Sigma_{t,g} s} + \frac{q_0}{\Sigma_{t,g}} \left( 1 - e^{-\Sigma_{t,g} s} \right)

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Random Rays
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the previous subsection, the govering characteristic equation along a 1D line through the reactor was written, such that an analytical solution for the ODE can be computed. If enough characteristic tracks (ODEs) are solved, then the behavior of the governing PDE can be numerically approximated. In traditional deterministic MOC, the selection of tracks has historically been a deterministic one, where azimuthal and polar quadratures are defined along with even track spacing in 3 dimensions. This is the point at which random ray diverges from deterministic MOC numerically. In Random Ray, rays are randomly sampled from a uniform distribution in space and angle and tracked along a set distance through the geometry before terminating. Importantly, different rays are sampled each power iteration, leading to a fully stochastic convergence process. I.e., inactive and active batches must be used, just as in Monte Carlo. While Monte Carlo implicitly converges the scattering source fully within each iteration, random ray (and MOC) solvers are not typically written to fully converge the scattering source within a single iteration. Rather, both the fission and scattering sources are updated each power iteration, thus requiring enough outer iterations so as to reach a stationary distribution in both the fission source and scattering source. I.e., even in a low dominance ration problem like a 2D pincell, several hundred inactive batches may still be required with random ray so as to allow the scattering source to fully develop, as neutrons undergoing hundreds of scatters may constitue a non-trivial contribution to the fission source.

Fundamentally, this distinction means that random ray typically requires more inactive iterations than are required in Monte Carlo, as the scattering source must also be developed. While a Monte Carlo simulation may only need 20-50 inactive iterations to reach a stationary source distribution for a full core light water reactor, a random ray solve will likely require 1,000 iterations or more. Source convergence metrics (e.g., Shannon Entropy) are thus highly useful tools when performing Random Ray simulations so as to help judge when the source has fully developed.

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ray Starting Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another key area of divergence between deterministic MOC and random ray is the starting conditions for rays. In deterministic MOC, the angular flux spectrum for rays are stored at any reflective or periodic boundaries so as to provide a starting condition for the next iteration. As there are many tracks, storage of angular fluxes can become costly in terms of memory consumption unless there are only vacuum boundaries present.

In random ray, as the starting locations of rays are sampled anew each iteration, the initial angular flux spectrum for the ray is unknown. While a guess can be made by taking the isotropic source from the FSR the ray was sampled in, direct usage of this quantity would result in significant bias and error being imparted on the simulation. f

Thus, an on-the-fly approximation method was developed (known as the "dead zone"), where the first several mean free paths of a ray are considered to be "inactive" or "read only". In this sense, the angular flux is solved for using the MOC equation, but the ray does not "tally" any scalar flux back to the FSRs that it travels through. After several mean free paths have been traversed, the ray's angular flux spectrum typically becomes dominated by the accumulated source terms from the cells it has travelled through, while the (incorrect) starting conditions have been attenuated away. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ray Ending Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To ensure that a uniform density of rays is integrated in space and angle throughout the simulation domain, after exiting the initial inactive "dead zone" portion of the ray, the rays are run for a user-specified distance. Typically, a choice of at least several times the length of the inactive "dead zone" is made so as to amortize the cost of the dead zone. E.g., if a dead zone of 30 cm is selected, then an active length of 300 cm might be selected so as to ensure the cost of the dead zone is below 10% of the overall runtime. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~
Power Iteration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simplified set of functions that execute a single random ray power iteration are given below. Not all global variables are defined in this illustrative example, but the high level components of the algorithm are shown. A number of significant simplifications are made for clarity -- for example, no inactive "dead zone" length is shown, geometry operations are abstracted, no parallelism (or thread safety) is expressed, among other subtleties.

The first block below shows the logic for the overall power iteration:

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

The second function shows the logic for transporting a single ray within the transport loop:

.. code-block:: C++

    double transport_single_ray(RandomRay& ray) {

        // Update source term (scattering + fission)
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

The final function below shows the logic for solving for the characteristic MOC equation (and accumulating the scalar flux contribution of the ray into the scalar flux value for the FSR).

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
            float delta_psi = (ray.angular_flux[e] - global::source[fsr][e]) * (1 - exp(-tau));
            ray.angular_flux_[e] -= delta_psi;
            global::scalar_flux_new[fsr][e] += delta_psi;
        }

        // Record total tracklength in this FSR (to compute volume)
        global::volume[fsr] += s;
    }

-----------------------------
How are Tallies Handled?
-----------------------------

Most tallies, filters, and scores that you would expect to work with a multigroup solver like random ray should work. E.g., you can define 3D mesh tallies with energy filters and flux, fission, and nu-fission scores, etc. There are some restrictions though. For starters, it is assumed that all filter mesh boundaries will conform to physical surface boundaries (or lattice boundaries) in the simulation geometry. It is acceptable for multiple cells (FSRs) to be contained within a filter mesh cell (e.g., pincell-level or assembly-level tallies should work), but it is currently left as undefined behavior if a single simulation cell is able to score to multiple filter mesh cells. In the future, we plan to add the capability to fully support mesh tallies, but for now this restriction needs to be respected.

-----------------------------
How is Plotting Handled?
-----------------------------

.. only:: html

   .. rubric:: References

.. [Askew] Askew, “A Characteristics Formulation of the Neutron Transport Equation in Complicated Geometries.” Technical Report AAEW-M 1108, UK Atomic Energy Establishment (1972).   

.. _Tramm 2017a: https://doi.org/10.1016/j.jcp.2017.04.038

.. _Tramm 2017b: https://doi.org/10.1016/j.anucene.2017.10.015

.. _Cosgrove 2023: https://doi.org/10.1080/00295639.2023.2270618

.. _Tramm 2018: http://hdl.handle.net/1721.1/119038