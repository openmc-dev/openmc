.. _io_materials:

========================================
Materials Specification -- materials.xml
========================================

.. _cross_sections:

----------------------------
``<cross_sections>`` Element
----------------------------

The ``<cross_sections>`` element has no attributes and simply indicates the path
to an XML cross section listing file (usually named cross_sections.xml). If this
element is absent from the settings.xml file, the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable will be used to find the
path to the XML cross section listing when in continuous-energy mode, and the
:envvar:`OPENMC_MG_CROSS_SECTIONS` environment variable will be used in
multi-group mode.

.. _multipole_library:

-------------------------------
``<multipole_library>`` Element
-------------------------------

The ``<multipole_library>`` element indicates the directory containing a
windowed multipole library. If a windowed multipole library is available,
OpenMC can use it for on-the-fly Doppler-broadening of resolved resonance range
cross sections. If this element is absent from the settings.xml file, the
:envvar:`OPENMC_MULTIPOLE_LIBRARY` environment variable will be used.

  .. note:: The <temperature_multipole> element must also be set to "true" for
    windowed multipole functionality.

.. _material:

----------------------
``<material>`` Element
----------------------

Each ``material`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the material.

  :name:
    An optional string name to identify the material in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :temperature:
    An element with no attributes which is used to set the default temperature
    of the material in Kelvin.

    *Default*: If a material default temperature is not given and a cell
    temperature is not specified, the :ref:`global default temperature
    <temperature_default>` is used.

  :density:
    An element with attributes/sub-elements called ``value`` and ``units``. The
    ``value`` attribute is the numeric value of the density while the ``units``
    can be "g/cm3", "kg/m3", "atom/b-cm", "atom/cm3", or "sum". The "sum" unit
    indicates that values appearing in ``ao`` or ``wo`` attributes for ``<nuclide>``
    and ``<element>`` sub-elements are to be interpreted as absolute nuclide/element
    densities in atom/b-cm or g/cm3, and the total density of the material is
    taken as the sum of all nuclides/elements. The "macro" unit is used with
    a ``macroscopic`` quantity to indicate that the density is already included
    in the library and thus not needed here.  However, if a value is provided
    for the ``value``, then this is treated as a number density multiplier on
    the macroscopic cross sections in the multi-group data.  This can be used,
    for example, when perturbing the density slightly.

    *Default*: None

    .. note:: A ``macroscopic`` quantity can not be used in conjunction with a
              ``nuclide``, ``element``, or ``sab`` quantity.

  :nuclide:
    An element with attributes/sub-elements called ``name``, and ``ao``
    or ``wo``. The ``name`` attribute is the name of the cross-section for a
    desired nuclide. Finally, the ``ao`` and ``wo`` attributes specify the atom or
    weight percent of that nuclide within the material, respectively. One
    example would be as follows:

    .. code-block:: xml

        <nuclide name="H1" ao="2.0" />
        <nuclide name="O16" ao="1.0" />

    .. note:: If one nuclide is specified in atom percent, all others must also
              be given in atom percent. The same applies for weight percentages.

    An optional attribute/sub-element for each nuclide is ``scattering``. This
    attribute may be set to "data" to use the scattering laws specified by the
    cross section library (default). Alternatively, when set to "iso-in-lab",
    the scattering laws are used to sample the outgoing energy but an
    isotropic-in-lab  distribution is used to sample the outgoing angle at each
    scattering interaction. The ``scattering`` attribute may be most useful
    when using OpenMC to compute multi-group cross-sections for deterministic
    transport codes and to quantify the effects of anisotropic scattering.

    *Default*: None

    .. note:: The ``scattering`` attribute/sub-element is not used in the
              multi-group :ref:`energy_mode`.

  :sab:
    Associates an S(a,b) table with the material. This element has an
    attribute/sub-element called ``name``. The ``name`` attribute
    is the name of the S(a,b) table that should be associated with the material.
    There is also an optional ``fraction`` element which indicates what fraction
    of the relevant nuclides will be affected by the S(a,b) table (e.g. which
    fraction of a material is crystalline versus amorphous).  ``fraction``
    defaults to unity.

    *Default*: None

    .. note:: This element is not used in the multi-group :ref:`energy_mode`.

  :macroscopic:
    The ``macroscopic`` element is similar to the ``nuclide`` element, but,
    recognizes that some multi-group libraries may be providing material
    specific macroscopic cross sections instead of always providing nuclide
    specific data like in the continuous-energy case.  To that end, the
    macroscopic element has one attribute/sub-element called ``name``.
    The ``name`` attribute is the name of the cross-section for a
    desired nuclide. One example would be as follows:

    .. code-block:: xml

        <macroscopic name="UO2" />

    .. note:: This element is only used in the multi-group :ref:`energy_mode`.

    *Default*: None
