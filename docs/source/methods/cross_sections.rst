.. _methods_cross_sections:

=============================
Cross Section Representations
=============================

----------------------
Continuous-Energy Data
----------------------

The data governing the interaction of neutrons with
various nuclei for continous-energy problems are represented using the ACE
format which is used by MCNP_ and Serpent_. ACE-format data can be generated
with the NJOY_ nuclear data processing system which converts raw
`ENDF/B data`_ into linearly-interpolable data as required by most Monte Carlo
codes. The use of a standard cross section format allows for a direct comparison
of OpenMC with other codes since the same cross section libraries can be used.

The ACE format contains continuous-energy cross sections for the following types
of reactions: elastic scattering, fission (or first-chance fission,
second-chance fission, etc.), inelastic scattering, :math:`(n,xn)`,
:math:`(n,\gamma)`, and various other absorption reactions. For those reactions
with one or more neutrons in the exit channel, secondary angle and energy
distributions may be provided. In addition, fissionable nuclides have total,
prompt, and/or delayed :math:`\nu` as a function of energy and neutron precursor
distributions. Many nuclides also have probability tables to be used for
accurate treatment of self-shielding in the unresolved resonance range. For
bound scatterers, separate tables with :math:`S(\alpha,\beta,T)` scattering law
data can be used.

Energy Grid Methods
-------------------

The method by which continuous energy cross sections for each nuclide in a
problem are stored as a function of energy can have a substantial effect on the
performance of a Monte Carlo simulation. Since the ACE format is based on
linearly-interpolable cross sections, each nuclide has cross sections tabulated
over a wide range of energies. Some nuclides may only have a few points
tabulated (e.g. H-1) whereas other nuclides may have hundreds or thousands of
points tabulated (e.g. U-238).

At each collision, it is necessary to sample the probability of having a
particular type of interaction whether it be elastic scattering, :math:`(n,2n)`,
level inelastic scattering, etc. This requires looking up the microscopic cross
sections for these reactions for each nuclide within the target material. Since
each nuclide has a unique energy grid, it would be necessary to search for the
appropriate index for each nuclide at every collision. This can become a very
time-consuming process, especially if there are many nuclides in a problem as
there would be for burnup calculations. Thus, there is a strong motive to
implement a method of reducing the number of energy grid searches in order to
speed up the calculation.

Logarithmic Mapping
+++++++++++++++++++

To speed up energy grid searches, OpenMC uses logarithmic mapping technique
[Brown]_ to limit the range of energies that must be searched for each
nuclide. The entire energy range is divided up into equal-lethargy segments, and
the bounding energies of each segment are mapped to bounding indices on each of
the nuclide energy grids. By default, OpenMC uses 8000 equal-lethargy segments
as recommended by Brown.

Other Methods
+++++++++++++

A good survey of other energy grid techniques, including unionized energy grids,
can be found in a paper by Leppanen_.

Windowed Multipole Representation
---------------------------------

In addition to the usual pointwise representation of cross sections, OpenMC
offers support for an experimental data format called windowed multipole (WMP).
This data format requires less memory than pointwise cross sections, and it
allows on-the-fly Doppler broadening to arbitrary temperature.

The multipole method was introduced by [Hwang]_ and the faster windowed
multipole method by [Josey]_.  In the multipole format, cross section resonances
are represented by poles, :math:`p_j`, and residues, :math:`r_j`, in the complex
plane.  The 0K cross sections in the resolved resonance region can be computed
by summing up a contribution from each pole:

.. math::
   \sigma(E, T=0\text{K}) = \frac{1}{E} \sum_j \text{Re} \left[
   \frac{i r_j}{\sqrt{E} - p_j} \right]

Assuming free-gas thermal motion, cross sections in the multipole form can be
analytically Doppler broadened to give the form:

.. math::
   \sigma(E, T) = \frac{1}{2 E \sqrt{\xi}} \sum_j \text{Re} \left[i r_j
   \sqrt{\pi} W_i(z) - \frac{r_j}{\sqrt{\pi}} C \left(\frac{p_j}{\sqrt{\xi}},
   \frac{u}{2 \sqrt{\xi}}\right)\right]
.. math::
   W_i(z) = \frac{i}{\pi} \int_{-\infty}^\infty dt \frac{e^{-t^2}}{z - t}
.. math::
   C \left(\frac{p_j}{\sqrt{\xi}},\frac{u}{2 \sqrt{\xi}}\right) =
   2p_j \int_0^\infty du' \frac{e^{-(u + u')^2/4\xi}}{p_j^2 - u'^2}
.. math::
   z = \frac{\sqrt{E} - p_j}{2 \sqrt{\xi}}
.. math::
   \xi = \frac{k_B T}{4 A}
.. math::
   u = \sqrt{E}

where :math:`T` is the temperature of the resonant scatterer, :math:`k_B` is the
Boltzmann constant, :math:`A` is the mass of the target nucleus. For
:math:`E \gg k_b T/A`, the :math:`C` integral is approximately zero, simplifying
the cross section to:

.. math::
   \sigma(E, T) = \frac{1}{2 E \sqrt{\xi}} \sum_j \text{Re} \left[i r_j
   \sqrt{\pi} W_i(z)\right]

The :math:`W_i` integral simplifies down to an analytic form.  We define the
Faddeeva function, :math:`W` as:

.. math::
   W(z) = e^{-z^2} \text{Erfc}(-iz)

Through this, the integral transforms as follows:

.. math::
   \text{Im} (z) > 0 : W_i(z) = W(z)
.. math::
   \text{Im} (z) < 0 : W_i(z) = -W(z^*)^*

There are freely available algorithms_ to evaluate the Faddeeva function. For
many nuclides, the Faddeeva function needs to be evaluated thousands of times to
calculate a cross section.  To mitigate that computational cost, the WMP method
only evaluates poles within a certain energy "window" around the incident
neutron energy and accounts for the effect of resonances outside that window
with a polynomial fit.  This polynomial fit is then broadened exactly. This
exact broadening can make up for the removal of the :math:`C` integral, as
typically at low energies, only curve fits are used.

Note that the implementation of WMP in OpenMC currently assumes that inelastic
scattering does not occur in the resolved resonance region.  This is usually,
but not always the case.  Future library versions may eliminate this issue.

The data format used by OpenMC to represent windowed multipole data is specified
in :ref:`io_data_wmp`.

----------------
Multi-Group Data
----------------

The data governing the interaction of particles with various nuclei or materials
are represented using a multi-group library format specific to the OpenMC code.
The format is described in the :ref:`mgxs_lib_spec`.
The data itself can be prepared via traditional paths or directly from a
continuous-energy OpenMC calculation by use of the Python API as is shown in the
:ref:`notebook_mgxs_part_iv` example notebook. This multi-group
library consists of meta-data (such as the energy group structure) and multiple
`xsdata` objects which contains the required microscopic or macroscopic
multi-group data.

At a minimum, the library must contain the absorption cross section
(:math:`\sigma_{a,g}`) and a scattering matrix. If the problem is an eigenvalue
problem then all fissionable materials must also contain either
a fission production matrix cross section
(:math:`\nu\sigma_{f,g\rightarrow g'}`), or
both the fission spectrum data (:math:`\chi_{g'}`) and a fission production
cross section (:math:`\nu\sigma_{f,g}`), or, .  The library must also contain
the fission cross section (:math:`\sigma_{f,g}`) or the fission energy release
cross section (:math:`\kappa\sigma_{f,g}`) if the associated tallies are
required by the model using the library.

After a scattering collision, the outgoing particle experiences a change in both
energy and angle. The probability of a particle resulting in a given outgoing
energy group (`g'`) given a certain incoming energy group (`g`) is provided
by the scattering matrix data.  The angular information can be expressed either
via Legendre expansion of the particle's change-in-angle (:math:`\mu`), a
tabular representation of the probability distribution function of :math:`\mu`,
or a histogram representation of the same PDF. The formats used to
represent these are described in the :ref:`mgxs_lib_spec`.

Unlike the continuous-energy mode, the multi-group mode does not explicitly
track particles produced from scattering multiplication (i.e., :math:`(n,xn)`)
reactions.  These are instead accounted for by adjusting the weight of the
particle after the collision such that the correct total weight is maintained.
The weight adjustment factor is optionally provided by the `multiplicity` data
which is required to be provided in the form of a group-wise matrix.
This data is provided as a group-wise matrix since the probability of producing
multiple particles in a scattering reaction depends on both the incoming energy,
`g`, and the sampled outgoing energy, `g'`. This data represents the average
number of particles emitted from a scattering reaction, given a scattering
reaction has occurred:

.. math::

    multiplicity_{g \rightarrow g'} = \frac{\nu_{scatter}\sigma_{s,g \rightarrow g'}}{
    								   \sigma_{s,g \rightarrow g'}}

If this scattering multiplication information is not provided in the library
then no weight adjustment will be performed. This is equivalent to neglecting
any additional particles produced in scattering multiplication reactions.
However, this assumption will result in a loss of accuracy since the total
particle population would not be conserved. This reduction in accuracy due to
the loss in particle conservation can be mitigated by reducing the absorption
cross section as needed to maintain particle conservation. This adjustment can
be done when generating the library, or by OpenMC. To have OpenMC perform the
adjustment, the total cross section (:math:`\sigma_{t,g}`) must be provided.
With this information, OpenMC will then adjust the absorption cross section as
follows:

.. math::

    \sigma_{a,g} = \sigma_{t,g} - \sum_{g'}\nu_{scatter}\sigma_{s,g \rightarrow g'}

The above method is the same as is usually done with most deterministic solvers.
Note that this method is less accurate than using the scattering multiplication
weight adjustment since simply reducing the absorption cross section does not
include any information about the outgoing energy of the particles produced in
these reactions.

All of the data discussed in this section can be provided to the code
independent of the particle's direction of motion (i.e., isotropic), or the data
can be provided as a tabular distribution of the polar and azimuthal particle
direction angles. The isotropic representation is the most commonly used,
however inaccuracies are to be expected especially near material interfaces
where a material has a very large cross sections relative to the other material
(as can be expected in the resonance range). The angular representation can be
used to minimize this error.

Finally, the above options for representing the physics do not have to be
consistent across the problem.  The number of groups and the structure, however,
does have to be consistent across the data sets. That is to say that each
microscopic or macroscopic data set does not have to apply the same scattering
expansion, treatment of multiplicity or angular representation of the cross
sections. This allows flexibility for the model to use highly anisotropic
scattering information in the water while the fuel can be simulated with linear
or even isotropic scattering.

.. only:: html

   .. rubric:: References

.. [Brown] Forrest B. Brown, "New Hash-based Energy Lookup Algorithm for Monte
           Carlo codes," LA-UR-14-24530, Los Alamos National Laboratory (2014).

.. [Hwang] R. N. Hwang, "A Rigorous Pole Representation of Multilevel Cross
           Sections and Its Practical Application,"  *Nucl. Sci. Eng.*, **96**,
           192-209 (1987).

.. [Josey] Colin Josey, Pablo Ducru, Benoit Forget, and Kord Smith, "Windowed
           Multipole for Cross Section Doppler Broadening," *J. Comp. Phys*,
           **307**, 715-727 (2016). http://dx.doi.org/10.1016/j.jcp.2015.08.013

.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _NJOY: http://t2.lanl.gov/codes.shtml
.. _ENDF/B data: http://www.nndc.bnl.gov/endf
.. _Leppanen: http://dx.doi.org/10.1016/j.anucene.2009.03.019
.. _algorithms: http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
