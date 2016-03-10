.. _usersguide_mgxs_library:

========================================
Multi-Group Cross Section Library Format
========================================

OpenMC can be run in continuous-energy mode or multi-group mode, provided the
nuclear data is available.  In continuous-energy mode, the
``cross_sections.xml`` file contains necessary meta-data for each data set,
including the name and a file system location where the complete library
can be found.  In multi-group mode, this ``cross_sections.xml`` file contains
this same meta-data describing the nuclide or material, but also contains the
group-wise nuclear data.  This portion of the manual describes the format of
the multi-group data library required to be used in the ``cross_sections.xml``
file.

Similar to the other input file types, the multi-group library is provided in
the XML_ format.  This library must provide some meta-data about the library
itself (such as the number of groups and the group structure, etc.) as well as
the actual cross section data itself for each of the necessary nuclides or
materials.

.. _XML: http://www.w3.org/XML/

------------------------------------------------
MGXS Library Specification -- cross_sections.xml
------------------------------------------------

The multi-group library meta-data is contained within the groups_,
group_structure_, and inverse_velocities_ elements.
The actual multi-group data itself is contained within the xsdata_ element.

.. _groups:

``<groups>`` Element
----------------------------------

The ``<groups>`` element has no attributes and simply provides the number of
energy groups contained within the library.

  *Default*: None, this must be provided.

.. _group_structure:

``<group_structure>`` Element
-----------------------------

The ``<group_structure>`` element has no attributes and should be provided as a
monotonically increasing list of bounding energies, in MeV, for a number of
groups. To provide proper energy boundaries, the length of the data within the
``<group_structure>`` element should be one more than the number of groups in
the problem.  For example, a two-group problem could be specified as:

.. code-block:: xml

    <group_structure> 0.0 0.625E-6 20.0 </group_structure>

*Default*: None, this must be provided.

.. _inverse_velocities:

``<inverse_velocities>`` Element
--------------------------------

The ``<inverse_velocities>`` element optionally indicates the average
inverse velocity corresponding to each of the groups in the problem.
This element should therefore be an array with a length which matches the
number of groups set in the groups_ element.

*Default*: Should this be needed by the presence of an ``inverse-velocity``
score in the ``tallies.xml`` file and not provided in this element, OpenMC
will simply convert the group mid-point energy to an inverse of the velocity
and use this information for tallying.

.. _xsdata:

``<xsdata>`` Element
--------------------

The ``<xsdata>`` element contains the nuclide or material-specific meta-data as
well as the actual cross section data. The following are the
attributes/sub-elements required to describe the meta-data:

  :name:
    The name of the microscopic or macroscopic data set.  An extension to the
    name must be provided (e.g., the ``.300K`` in ``UO2.300K``).  The name and
    extension together must be twelve or less characters in length.  This
    extension must follow a period and be five characters or less in length.
    similar to the equivalent in the continuous-energy ``cross_sections.xml``
    file, is used to denote variants of the particular nuclide or material of
    interest (i.e. the ``UO2`` data in this example could have been generated
    at a temperature of 300K).

    *Default*: None, this must be provided.

  :alias:
    An alternative name to use for the microscopic or macroscopic data set.

    *Default*: If no alias is provided, it will adopt the value of ``name``.

  :kT:
    The temperature times Boltzmann's constant (in units of MeV) at which the
    data was generated.

    *Default*: Room temperature, 2.53E-8 MeV

  :fissionable:
    This element states whether or not the data in question is fissionable.
    Accepted values are "true" or "false".

    *Default*: None, this element must be provided.

  :representation:
    This element provides the method used to generate and represent the
    multi-group cross sections.  That is, whether they were generated with
    scalar flux weighting (or reduced to an equivalent representation)
    and thus are angle-independent, or if the data was generated with angular
    dependent fluxes and thus the data is angle-dependent.  The options are
    either "isotropic" or "angle".

    *Default*: "isotropic"

  :num_azimuthal:
    This element provides the number of equal width angular bins that the
    azimuthal angular domain is subdivided in the case of angle-dependent
    cross sections (i.e., "angle" is passed to the ``representation`` element).
    Note that these bins are equal in azimuthal angle widths, not equal in the
    cosine of the azimuthal angle widths.

    *Default*: If ``representation`` is "angle", this must be provided.  This
    parameter is not used for other ``representation`` types.

  :num_polar:
    This element provides the number of equal width angular bins that the
    polar angular domain is subdivided in the case of angle-dependent
    cross sections (i.e., "angle" is passed to the ``representation`` element).
    Note that these bins are equal in polar angle widths, not equal in the
    cosine of the polar angle widths.


    *Default*: If ``representation`` is "angle", this must be provided.  This
    parameter is not used for other ``representation`` types.

  :scatt_type:
    This element provides the representation of the angular distribution
    associated with each group-to-group transfer probability. The options are
    either "legendre", "histogram", or "tabular".
    The "legendre" option means the angular distribution has been
    expanded via Legendre polynomials of the order provided in the "order"
    element.
    The "histogram" option means the angular distribution is provided in
    an equi-width histogram format with a number of bins as provided in the
    "order" element.  This is useful when the angular distribution was
    obtained from a Monte Carlo tally and thus is natively in the histogram
    format.
    The "tabular" option means the angular distribution is provided in an
    equi-spaced point-wise representation.

    *Default*: "legendre"

  :order:
    This element provides either the Legendre order, number of bins, or number
    of points used to describe the angular distribution associated with each
    group-to-group transfer probability.  The specific meaning of this bin
    depends upon the value of ``scatt_type`` as discussed above.

    *Default*: None, this element must be provided.

  :tabular_legendre:
    This optional element is used to set how the Legendre scattering kernel, if
    provided via the ``scatt_type`` element above, is represented and thus used
    during the scattering process.  Specifically, the options are to either
    convert the Legendre expansion to a tabular representation or leave it as
    a set of Legendre coefficients.  Converting to a tabular representation will
    cost memory but is likely to decrease runtime compared to leaving as a
    set of Legendre coefficients.  This element has the following
    attributes/sub-elements:

    :enable:
      This attribute/sub-element denotes whether or not the conversion to the
      tabular format should be performed or not.  A value of "true" means
      the conversion should be performed, "false" means it should not.

      *Default*: "true"

    :num_points:
      If the conversion is to take place the number of tabular points is
      required.  This attribute/sub-element allows the user to set the desired
      number of points.

      *Default*: 33

  The following attributes/sub-elements are the cross section values to
  be used during the transport process.

  :total:
    This element requires the group-wise total cross section ordered by
    increasing group index (i.e., fast to thermal).  If ``representation`` is
    "isotropic", then the length of this list should equal the number of
    groups described in the ``groups`` element.  If ``representation`` is
    "angle", then the length of this list should equal the number of groups
    times the number of azimuthal angles times the number of polar angles,
    with the inner-dimension being groups, intermediate-dimension being
    azimuthal angles and outer-dimension being the polar angles.

    *Default*: If not provided, it will be determined by summing the
    absorption and scattering cross sections.

  :absorption:
    This element requires the group-wise absorption cross section ordered by
    increasing group index (i.e., fast to thermal).  If ``representation`` is
    "isotropic", then the length of this list should equal the number of
    groups described in the ``groups`` element.  If ``representation`` is
    "angle", then the length of this list should equal the number of groups
    times the number of azimuthal angles times the number of polar angles,
    with the inner-dimension being groups, intermediate-dimension being
    azimuthal angles and outer-dimension being the polar angles.

    *Default*: None, this must be provided.

  :scatter:
    This element requires the scattering moment matrices presented with the
    columns representing incoming group and rows representing the outgoing
    group.  That is, down-scatter will be above the diagonal of the resultant
    matrix.  This matrix is repeated for every Legendre order (in order of
    increasing orders) if ``scatt_type`` is "legendre"; otherwise, this
    matrix is repeated for every bin of the histogram or tabular
    representation.  Finally, if ``representation`` is "angle", the above
    is repeated for every azimuthal angle and every polar angle, in that
    order.

    *Default*: None, this must be provided.

  :multiplicity:
    This element provides the ratio of neutrons produced in scattering
    collisions to the neutrons which undergo scattering collisions; that is,
    the multiplicity provides the code with a scaling factor to account for
    neutrons being produced in (n,xn) reactions.  This information is assumed
    isotropic and therefore does not need to be repeated for every Legendre
    moment or histogram/tabular bin.  This matrix follows the same arrangement
    as described for the ``scatter`` element, with the exception of the
    data needed to provide the scattering type information.

    *Default*: Multiplicities of 1.0 are assumed (i.e., (n,xn) reactions are
    neglected).

  The following fission-specific data are only needed should ``fissionable``
  be "true".

  :fission:
    This element requires the group-wise fission cross section ordered by
    increasing group index (i.e., fast to thermal).  If ``representation`` is
    "isotropic", then the length of this list should equal the number of
    groups described in the ``groups`` element.  If ``representation`` is
    "angle", then the length of this list should equal the number of groups
    times the number of azimuthal angles times the number of polar angles,
    with the inner-dimension being groups, intermediate-dimension being
    azimuthal angles and outer-dimension being the polar angles.

    *Default*: None, this is required only if fission tallies are
    requested and the material is fissionable.

  :kappa_fission:
    This element requires the group-wise kappa-fission cross section ordered by
    increasing group index (i.e., fast to thermal).  If ``representation`` is
    "isotropic", then the length of this list should equal the number of
    groups described in the ``groups`` element.  If ``representation`` is
    "angle", then the length of this list should equal the number of groups
    times the number of azimuthal angles times the number of polar angles,
    with the inner-dimension being groups, intermediate-dimension being
    azimuthal angles and outer-dimension being the polar angles.

    *Default*: None, this is required only if kappa_fission tallies are
    requested and the material is fissionable.

  :chi:
    This element requires the group-wise fission spectra ordered by
    increasing group index (i.e., fast to thermal).  This element should be
    used if making the common approximation that the fission spectra does
    not depend on incoming energy.  If the user does not wish to make this
    approximation, then this should not be provided and this information
    included in the ``nu_fission`` element instead.  If ``representation`` is
    "isotropic", then the length of this list should equal the number of
    groups described in the ``groups`` element.  If ``representation`` is
    "angle", then the length of this list should equal the number of groups
    times the number of azimuthal angles times the number of polar angles,
    with the inner-dimension being groups, intermediate-dimension being
    azimuthal angles and outer-dimension being the polar angles.

    *Default*: None, either this element is provided or ``nu_fission`` is
    provided in fission matrix form, or the material is not fissionable.

  :nu_fission:
    This element provides either the group-wise fission production cross
    section vector (i.e., if ``chi`` is provided), or is the group-wise fission
    production matrix.  If providing the vector, it should be ordered the same
    as the ``fission`` data.  If providing the matrix, it should be ordered
    the same as the ``multiplicity`` matrix.

    *Default*: None, either this element must be provided if the material
    is fissionable.

