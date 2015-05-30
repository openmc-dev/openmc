.. _methods_cross_sections:

============================
Cross Section Representation
============================

The data governing the interaction of neutrons with various nuclei are
represented using the ACE format which is used by MCNP_ and Serpent_. ACE-format
data can be generated with the NJOY_ nuclear data processing system which
converts raw `ENDF/B data`_ into linearly-interpolable data as required by most
Monte Carlo codes. The use of a standard cross section format allows for a
direct comparison of OpenMC with other codes since the same cross section
libraries can be used.

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

-------------------
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
-------------------

To speed up energy grid searches, OpenMC uses logarithmic mapping technique
[Brown]_ to limit the range of energies that must be searched for each
nuclide. The entire energy range is divided up into equal-lethargy segments, and
the bounding energies of each segment are mapped to bounding indices on each of
the nuclide energy grids. By default, OpenMC uses 8000 equal-lethargy segments
as recommended by Brown.

Other Methods
-------------

A good survey of other energy grid techniques, including unionized energy grids,
can be found in a paper by Leppanen_.

.. only:: html

   .. rubric:: References

.. [Brown] Forrest B. Brown, "New Hash-based Energy Lookup Algorithm for Monte
           Carlo codes," LA-UR-14-24530, Los Alamos National Laboratory (2014).

.. _MCNP: http://mcnp.lanl.gov
.. _Serpent: http://montecarlo.vtt.fi
.. _NJOY: http://t2.lanl.gov/codes.shtml
.. _ENDF/B data: http://www.nndc.bnl.gov/endf
.. _Leppanen: http://dx.doi.org/10.1016/j.anucene.2009.03.019
