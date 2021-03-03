.. _io_tallies:

====================================
Tallies Specification -- tallies.xml
====================================

The tallies.xml file allows the user to tell the code what results he/she is
interested in, e.g. the fission rate in a given cell or the current across a
given surface. There are two pieces of information that determine what
quantities should be scored. First, one needs to specify what region of phase
space should count towards the tally and secondly, the actual quantity to be
scored also needs to be specified. The first set of parameters we call *filters*
since they effectively serve to filter events, allowing some to score and
preventing others from scoring to the tally.

The structure of tallies in OpenMC is flexible in that any combination of
filters can be used for a tally. The following types of filter are available:
cell, universe, material, surface, birth region, pre-collision energy,
post-collision energy, and an arbitrary structured mesh.

The five valid elements in the tallies.xml file are ``<tally>``, ``<filter>``,
``<mesh>``, ``<derivative>``, and ``<assume_separate>``.

.. _tally:

-------------------
``<tally>`` Element
-------------------

The ``<tally>`` element accepts the following sub-elements:

  :name:
    An optional string name to identify the tally in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :filters:
    A space-separated list of the IDs of ``filter`` elements.

  :nuclides:
    If specified, the scores listed will be for particular nuclides, not the
    summation of reactions from all nuclides. The format for nuclides should be
    [Atomic symbol]-[Mass number], e.g. "U-235". The reaction rate for all
    nuclides can be obtained with "total". For example, to obtain the reaction
    rates for U-235, Pu-239, and all nuclides in a material, this element should
    be:

    .. code-block:: xml

        <nuclides>U-235 Pu-239 total</nuclides>

    *Default*: total

  :estimator:
    The estimator element is used to force the use of either ``analog``,
    ``collision``, or ``tracklength`` tally estimation.  ``analog`` is generally
    the least efficient though it can be used with every score type.
    ``tracklength`` is generally the most efficient, but neither ``tracklength``
    nor ``collision`` can be used to score a tally that requires post-collision
    information.  For example, a scattering tally with outgoing energy filters
    cannot be used with ``tracklength`` or ``collision`` because the code will
    not know the outgoing energy distribution.

    *Default*: ``tracklength`` but will revert to ``analog`` if necessary.

  :scores:
    A space-separated list of the desired responses to be accumulated. A full
    list of valid scores can be found in the :ref:`user's guide
    <usersguide_scores>`.

  :trigger:
    Precision trigger applied to all filter bins and nuclides for this tally.
    It must specify the trigger's type, threshold and scores to which it will
    be applied. It has the following attributes/sub-elements:

   :type:
     The type of the trigger. Accepted options are "variance", "std_dev",
     and "rel_err".

     :variance:
       Variance of the batch mean :math:`\sigma^2`

     :std_dev:
       Standard deviation of the batch mean :math:`\sigma`

     :rel_err:
       Relative error of the batch mean :math:`\frac{\sigma}{\mu}`

     *Default*: None

   :threshold:
     The precision trigger's convergence criterion for tallied values.

     *Default*: None

   :scores:
     The score(s) in this tally to which the trigger should be applied.

     .. note:: The ``scores`` in ``trigger`` must have been defined in
               ``scores`` in ``tally``. An optional "all" may be used to
               select all scores in this tally.

     *Default*: "all"

  :derivative:
    The id of a ``derivative`` element. This derivative will be applied to all
    scores in the tally. Differential tallies are currently only implemented
    for collision and analog estimators.

     *Default*: None


--------------------
``<filter>`` Element
--------------------

Filters can be used to modify tally behavior. Most tallies (e.g. ``cell``,
``energy``, and ``material``) restrict the tally so that only particles
within certain regions of phase space contribute to the tally.  Others
(e.g. ``delayedgroup`` and ``energyfunction``) can apply some other function
to the scored values. The ``filter`` element has the following
attributes/sub-elements:

  :type:
    The type of the filter. Accepted options are "cell", "cellfrom",
    "cellborn", "surface", "material", "universe", "energy", "energyout", "mu",
    "polar", "azimuthal", "mesh", "distribcell", "delayedgroup",
    "energyfunction", and "particle".

  :bins:
     A description of the bins for each type of filter can be found in
     :ref:`filter_types`.

  :energy:
    ``energyfunction`` filters multiply tally scores by an arbitrary
    function. The function is described by a piecewise linear-linear set of
    (energy, y) values. This entry specifies the energy values. The function
    will be evaluated as zero outside of the bounds of this energy grid.
    (Only used for ``energyfunction`` filters)

  :y:
    ``energyfunction`` filters multiply tally scores by an arbitrary
    function. The function is described by a piecewise linear-linear set of
    (energy, y) values. This entry specifies the y values. (Only used
    for ``energyfunction`` filters)

.. _filter_types:

Filter Types
++++++++++++

For each filter type, the following table describes what the ``bins`` attribute
should be set to:

:cell:
  A list of unique IDs for cells in which the tally should be
  accumulated.

:surface:
  This filter allows the tally to be scored when crossing a surface. A list of
  surface IDs should be given. By default, net currents are tallied, and to
  tally a partial current from one cell to another, this should be used in
  combination with a cell or cell_from filter that defines the other cell.
  This filter should not be used in combination with a meshfilter.

:cellfrom:
  This filter allows the tally to be scored when crossing a surface and the
  particle came from a specified cell. A list of cell IDs should be
  given.
  To tally a partial current from a cell to another, this filter should be
  used in combination with a cell filter, to define the other cell.
  This filter should not be used in combination with a meshfilter.

:cellborn:
  This filter allows the tally to be scored to only when particles were
  originally born in a specified cell. A list of cell IDs should be
  given.

:material:
  A list of unique IDs for materials in which the tally should be accumulated.

:universe:
  A list of unique IDs for universes in which the tally should be accumulated.

:energy:
  In continuous-energy mode, this filter should be provided as a
  monotonically increasing list of bounding **pre-collision** energies
  for a number of groups. For example, if this filter is specified as

  .. code-block:: xml

      <filter type="energy" bins="0.0 1.0e6 20.0e6" />

  then two energy bins will be created, one with energies between 0 and
  1 MeV and the other with energies between 1 and 20 MeV.

  In multi-group mode the bins provided must match group edges
  defined in the multi-group library.

:energyout:
  In continuous-energy mode, this filter should be provided as a
  monotonically increasing list of bounding **post-collision** energies
  for a number of groups. For example, if this filter is specified as

  .. code-block:: xml

      <filter type="energyout" bins="0.0 1.0e6 20.0e6" />

  then two post-collision energy bins will be created, one with
  energies between 0 and 1 MeV and the other with energies between
  1 and 20 MeV.

  In multi-group mode the bins provided must match group edges
  defined in the multi-group library.

:mu:
  A monotonically increasing list of bounding **post-collision** cosines
  of the change in a particle's angle (i.e., :math:`\mu = \hat{\Omega}
  \cdot \hat{\Omega}'`), which represents a portion of the possible
  values of :math:`[-1,1]`.  For example, spanning all of :math:`[-1,1]`
  with five equi-width bins can be specified as:

  .. code-block:: xml

      <filter type="mu" bins="-1.0 -0.6 -0.2 0.2 0.6 1.0" />

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[-1,1]` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[-1,1]` can be instead written as:

  .. code-block:: xml

      <filter type="mu" bins="5" />

:polar:
  A monotonically increasing list of bounding particle polar angles
  which represents a portion of the possible values of :math:`[0,\pi]`.
  For example, spanning all of :math:`[0,\pi]` with five equi-width
  bins can be specified as:

  .. code-block:: xml

      <filter type="polar" bins="0.0 0.6283 1.2566 1.8850 2.5132 3.1416"/>

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[0,\pi]` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[0,\pi]` can be instead written as:

  .. code-block:: xml

      <filter type="polar" bins="5" />

:azimuthal:
  A monotonically increasing list of bounding particle azimuthal angles
  which represents a portion of the possible values of :math:`[-\pi,\pi)`.
  For example, spanning all of :math:`[-\pi,\pi)` with two equi-width
  bins can be specified as:

  .. code-block:: xml

      <filter type="azimuthal" bins="0.0 3.1416 6.2832" />

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[-\pi,\pi)` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[-\pi,\pi)` can be instead written as:

  .. code-block:: xml

      <filter type="azimuthal" bins="2" />

:mesh:
  The unique ID of a mesh to be tallied over.

:distribcell:
  The single cell which should be tallied uniquely for all instances.

  .. note:: The distribcell filter will take a single cell ID and will tally
            each unique occurrence of that cell separately. This filter will not
            accept more than one cell ID. It is not recommended to combine this
            filter with a cell or mesh filter.

:delayedgroup:
  A list of delayed neutron precursor groups for which the tally should
  be accumulated. For instance, to tally to all 6 delayed groups in the
  ENDF/B-VII.1 library the filter is specified as:

  .. code-block:: xml

      <filter type="delayedgroup" bins="1 2 3 4 5 6" />

:energyfunction:
  ``energyfunction`` filters do not use the ``bins`` entry.  Instead
  they use ``energy`` and ``y``.

:particle:
  A list of integers indicating the type of particles to tally ('neutron' = 1,
  'photon' = 2, 'electron' = 3, 'positron' = 4).

------------------
``<mesh>`` Element
------------------

If a mesh is desired as a filter for a tally, it must be specified in a separate
element with the tag name ``<mesh>``. This element has the following
attributes/sub-elements:

  :type:
    The type of mesh. This can be either "regular", "rectilinear", or
    "unstructured".

  :dimension:
    The number of mesh cells in each direction. (For regular mesh only.)

  :lower_left:
    The lower-left corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh. (For regular mesh only.)

  :upper_right:
    The upper-right corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh. (For regular mesh only.)

  :width:
    The width of mesh cells in each direction. (For regular mesh only.)

  :x_grid:
    The mesh divisions along the x-axis. (For rectilinear mesh only.)

  :y_grid:
    The mesh divisions along the y-axis. (For rectilinear mesh only.)

  :z_grid:
    The mesh divisions along the z-axis. (For rectilinear mesh only.)

  :mesh_file:
    The name of the mesh file to be loaded at runtime. (For unstructured mesh
    only.)

  .. note::
      One of ``<upper_right>`` or ``<width>`` must be specified, but not both
      (even if they are consistent with one another).

------------------------
``<derivative>`` Element
------------------------

OpenMC can take the first-order derivative of many tallies with respect to
material perturbations. It works by propagating a derivative through the
transport equation. Essentially, OpenMC keeps track of how each particle's
weight would change as materials are perturbed, and then accounts for that
weight change in the tallies. Note that this assumes material perturbations are
small enough not to change the distribution of fission sites. This element has
the following attributes/sub-elements:

  :id:
    A unique integer that can be used to identify the derivative.

  :variable:
    The independent variable of the derivative. Accepted options are "density",
    "nuclide_density", and "temperature". A "density" derivative will give the
    derivative with respect to the density of the material in [g / cm^3]. A
    "nuclide_density" derivative will give the derivative with respect to the
    density of a particular nuclide in units of [atom / b / cm].  A
    "temperature" derivative is with respect to a material temperature in units
    of [K].  The temperature derivative requires windowed multipole to be
    turned on.  Note also that the temperature derivative only accounts for
    resolved resonance Doppler broadening.  It does not account for thermal
    expansion, S(a, b) scattering, resonance scattering, or unresolved Doppler
    broadening.

  :material:
    The perturbed material. (Necessary for all derivative types)

  :nuclide:
    The perturbed nuclide. (Necessary only for "nuclide_density")

-----------------------------
``<assume_separate>`` Element
-----------------------------

In cases where the user needs to specify many different tallies each of which
are spatially separate, this tag can be used to cut down on some of the tally
overhead. The effect of assuming all tallies are spatially separate is that once
one tally is scored to, the same event is assumed not to score to any other
tallies. This element should be followed by "true" or "false".

  .. warning:: If used incorrectly, the assumption that all tallies are
               spatially separate can lead to incorrect results.

  *Default*: false
