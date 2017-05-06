.. _io_cmfd:

==============================
CMFD Specification -- cmfd.xml
==============================

Coarse mesh finite difference acceleration method has been implemented in
OpenMC. Currently, it allows users to accelerate fission source convergence
during inactive neutron batches. To run CMFD, the ``<run_cmfd>`` element in
``settings.xml`` should be set to "true".

-------------------
``<begin>`` Element
-------------------

The ``<begin>`` element controls what batch CMFD calculations should begin.

  *Default*: 1

------------------------
``<dhat_reset>`` Element
------------------------

The ``<dhat_reset>`` element controls whether :math:`\widehat{D}` nonlinear
CMFD parameters should be reset to zero before solving CMFD eigenproblem.
It can be turned on with "true" and off with "false".

  *Default*: false

---------------------
``<display>`` Element
---------------------

The ``<display>`` element sets one additional CMFD output column. Options are:

* "balance" - prints the RMS [%] of the resdiual from the neutron balance
  equation on CMFD tallies.
* "dominance" - prints the estimated dominance ratio from the CMFD iterations.
  **This will only work for power iteration eigensolver**.
* "entropy" - prints the *entropy* of the CMFD predicted fission source.
  **Can only be used if OpenMC entropy is active as well**.
* "source" - prints the RMS [%] between the OpenMC fission source and CMFD
  fission source.

  *Default*: balance

-------------------------
``<downscatter>`` Element
-------------------------

The ``<downscatter>`` element controls whether an effective downscatter cross
section should be used when using 2-group CMFD. It can be turned on with "true"
and off with "false".

  *Default*: false

----------------------
``<feedback>`` Element
----------------------

The ``<feedback>`` element controls whether or not the CMFD diffusion result is
used to adjust the weight of fission source neutrons on the next OpenMC batch.
It can be turned on with "true" and off with "false".

  *Default*: false

------------------------------------
``<gauss_seidel_tolerance>`` Element
------------------------------------

The ``<gauss_seidel_tolerance>`` element specifies two parameters. The first is
the absolute inner tolerance for Gauss-Seidel iterations when performing CMFD
and the second is the relative inner tolerance for Gauss-Seidel iterations
for CMFD calculations.

  *Default*: 1.e-10 1.e-5

--------------------
``<ktol>`` Element
--------------------

The ``<ktol>`` element specifies the tolerance on the eigenvalue when performing
CMFD power iteration.

  *Default*: 1.e-8

------------------
``<mesh>`` Element
------------------

The CMFD mesh is a structured Cartesian mesh. This element has the following
attributes/sub-elements:

  :lower_left:
    The lower-left corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :upper_right:
    The upper-right corner of the structrued mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :dimension:
    The number of mesh cells in each direction.

  :width:
    The width of mesh cells in each direction.

  :energy:
    Energy bins [in eV], listed in ascending order (e.g. 0.0 0.625 20.0e6)
    for CMFD tallies and acceleration. If no energy bins are listed, OpenMC
    automatically assumes a one energy group calculation over the entire
    energy range.

  :albedo:
    Surface ratio of incoming to outgoing partial currents on global boundary
    conditions. They are listed in the following order: -x +x -y +y -z +z.

    *Default*: 1.0 1.0 1.0 1.0 1.0 1.0

  :map:
    An optional acceleration map can be specified to overlay on the coarse
    mesh spatial grid. If this option is used, a ``1`` is used for a
    non-accelerated region and a ``2`` is used for an accelerated region.
    For a simple 4x4 coarse mesh with a 2x2 fuel lattice surrounded by
    reflector, the map is:

      ``1 1 1 1``

      ``1 2 2 1``

      ``1 2 2 1``

      ``1 1 1 1``

    Therefore a 2x2 system of equations is solved rather than a 4x4. This
    is extremely important to use in reflectors as neutrons will not
    contribute to any tallies far away from fission source neutron regions.
    A ``2`` must be used to identify any fission source region.

    .. note:: Only two of the following three sub-elements are needed:
              ``lower_left``, ``upper_right`` and ``width``. Any combination
              of two of these will yield the third.

------------------
``<norm>`` Element
------------------

The ``<norm>`` element is used to normalize the CMFD fission source distribution
to a particular value. For example, if a fission source is calculated for a
17 x 17 lattice of pins, the fission source may be normalized to the number of
fission source regions, in this case 289. This is useful when visualizing this
distribution as the average peaking factor will be unity. This parameter will
not impact the calculation.

  *Default*: 1.0

---------------------------
``<power_monitor>`` Element
---------------------------

The ``<power_monitor>`` element is used to view the convergence of power
iteration. This option can be turned on with "true" and turned off with "false".

  *Default*: false

-------------------------
``<run_adjoint>`` Element
-------------------------

The ``<run_adjoint>`` element can be turned on with "true" to have an adjoint
calculation be performed on the last batch when CMFD is active.

  *Default*: false

--------------------
``<shift>`` Element
--------------------

The ``<shift>`` element specifies an optional Wielandt shift parameter for
accelerating power iterations. It is by default very large so the impact of the
shift is effectively zero.

  *Default*: 1e6

----------------------
``<spectral>`` Element
----------------------

The ``<spectral>`` element specifies an optional spectral radius that can be set to
accelerate the convergence of Gauss-Seidel iterations during CMFD power iteration
solve.

  *Default*: 0.0

------------------
``<stol>`` Element
------------------

The ``<stol>`` element specifies the tolerance on the fission source when performing
CMFD power iteration.

  *Default*: 1.e-8

-------------------------
``<tally_reset>`` Element
-------------------------

The ``<tally_reset>`` element contains a list of batch numbers in which CMFD tallies
should be reset.

  *Default*: None

----------------------------
``<write_matrices>`` Element
----------------------------

The ``<write_matrices>`` element is used to write the sparse matrices created
when solving CMFD equations. This option can be turned on with "true" and off
with "false".

  *Default*: false
