.. _usersguide_input:

=======================
Writing XML Input Files
=======================

Unlike many other Monte Carlo codes which use an arbitrary-format ASCII file
with "cards" to specify a particular geometry, materials, and associated run
settings, the input files for OpenMC are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.

Anyone who has ever seen webpages written in HTML will be familiar with the
structure of XML whereby "tags" enclosed in angle brackets denote that a
particular piece of data will follow. Let us examine the follow example:

.. code-block:: xml

    <person>
      <firstname>John</firstname>
      <lastname>Smith</lastname>
      <age>27</age>
      <occupation>Health Physicist</occupation>
    </person>

Here we see that the first tag indicates that the following data will describe a
person. The nested tags *firstname*, *lastname*, *age*, and *occupation*
indicate characteristics about the person being described.

In much the same way, OpenMC input uses XML tags to describe the geometry, the
materials, and settings for a Monte Carlo simulation.

.. _XML: http://www.w3.org/XML/

-----------------
Overview of Files
-----------------

To assemble a complete model for OpenMC, one needs to create separate XML files
for the geometry, materials, and settings. Additionally, there are two optional
input files. The first is a tallies XML file that specifies physical quantities
to be tallied. The second is a plots XML file that specifies regions of geometry
which should be plotted. OpenMC expects that these files are called:

* ``geometry.xml``
* ``materials.xml``
* ``setings.xml``
* ``tallies.xml``
* ``plots.xml``

--------------------------------------
Settings Specification -- settings.xml
--------------------------------------

All simulation parameters and miscellaneous options are specified in the
settings.xml file.

``<criticality>`` Element
-------------------------

The ``<criticality>`` element indicates that a criticality calculation should be
performed. It has the following attributes/sub-elements:

  :batches: 
    The total number of batches, where each batch corresponds to multiple
    fission source iterations. Batching is done to eliminate correlation between
    realizations of random variables.

    *Default*: None

  :generations_per_batch:
    The number of total fission source iterations per batch.

    *Default*: 1

  :inactive:
    The number of inactive batches. In general, the starting cycles in a
    criticality calculation can not be used to contribute to tallies since the
    fission source distribution and eigenvalue are generally not converged
    immediately.

    *Default*: None

  :particles:
    The number of neutrons to simulate per fission source iteration.

    *Default*: None

.. _cross_sections:

``<cross_sections>`` Element
----------------------------

The ``<cross_sections>`` element has no attributes and simply indicates the path
to an XML cross section listing file (usually named cross_sections.xml). If this
element is absent from the settings.xml file, the :envvar:`CROSS_SECTIONS`
environment variable will be used to find the path to the XML cross section
listing.

``<cutoff>`` Element
--------------------

The ``<cutoff>`` element indicates the weight cutoff used below which particles
undergo Russian roulette. Surviving particles are assigned a user-determined
weight. Note that weight cutoffs and Russian rouletting are not turned on by
default. This element has the following attributes/sub-elements:

  :weight:
    The weight below which particles undergo Russian roulette.

    *Default*: 0.25

  :weight_avg:
    The weight that is assigned to particles that are not killed after Russian
    roulette.

    *Default*: 1.0

``<energy_grid>`` Element
-------------------------

The ``<energy_grid>`` element determines the treatment of the energy grid during
a simulation. Setting this element to "nuclide" will cause OpenMC to use a
nuclide's energy grid when determining what points to interpolate between for
determining cross sections (i.e. non-unionized energy grid). To use a unionized
energy grid, set this element to "union". Note that the unionized energy grid
treatment is slightly different than that employed in Serpent.

  *Default*: union

``<entropy>`` Element
---------------------

The ``<entropy>`` element describes a mesh that is used for calculting Shannon
entropy. This mesh should cover all possible fissionable materials in the
problem. It has the following attributes/sub-elements:

  :dimension:
    The number of mesh cells in the x, y, and z directions, respectively.

    *Default*: If this tag is not present, the number of mesh cells is
     automatically determined by the code.

  :lower_left:
    The Cartersian coordinates of the lower-left corner of the mesh.

    *Default*: None

  :upper_right:
    The Cartersian coordinates of the upper-right corner of the mesh.

    *Default*: None

``<ptables>`` Element
---------------------

The ``<ptables>`` element determines whether probability tables should be used
in the unresolved resonance range if available. This element has no attributes
or sub-elements and can be set to either "off" or "on".

  *Default*: on

``<seed>`` Element
------------------

The ``seed`` element is used to set the seed used for the linear congruential
pseudo-random number generator.

  *Default*: 1

``<source>`` Element
--------------------

The ``source`` element gives information on an initial source guess for
criticality calculations. It takes the following attributes:

  :type:
    The type of source distribution. Setting this to "box" indicates that the
    starting source should be sampled uniformly in a parallelepiped. Setting
    this to "point" indicates that the starting source should be sampled from an
    isotropic point source. Setting this to "file" indicates that the starting
    source should be sampled from a ``source.binary`` file.

  :coeffs:
    For a "box" source distribution, ``coeffs`` should be given as six real
    numbers, the first three of which specify the lower-left corner of a
    parallelepiped and the last three of which specify the upper-right
    corner. Source sites are sampled uniformly through that parallelepiped.

    For a "point" source distribution, ``coeffs`` should be given as three real
    numbers which specify the (x,y,z) location of an isotropic point source

    For a "file" source distribution, ``coeffs`` should not be specified.

``<survival_biasing>`` Element
------------------------------

The ``<survival_biasing>`` element has no attributes and assumes wither the
value ``on`` or ``off``. If turned on, this option will enable the use of
survival biasing, otherwise known as implicit capture or absorption.

  *Default*: off

``<trace>`` Element
-------------------

The ``<trace>`` element can be used to print out detailed information about a
single particle during a simulation. This element should be followed by three
integers: the batch number, generation number, and particle number.

  *Default*: None

``<verbosity>`` Element
-----------------------

The ``<verbosity>`` element tells the code how much information to display to
the standard output. A higher verbosity corresponds to more information being
displayed. This element takes the following attributes:

  :value:
    The specified verbosity between 1 and 10.

    *Default*: 5

--------------------------------------
Geometry Specification -- geometry.xml
--------------------------------------

The geometry in OpenMC is described using `constructive solid geometry`_ (CSG),
also sometimes referred to as combinatorial geometry. CSG allows a user to
create complex objects using Boolean operators on a set of simpler surfaces. In
the geometry model, each unique closed volume in defined by its bounding
surfaces. In OpenMC, most `quadratic surfaces`_ can be modeled and used as
bounding surfaces.

Every geometry.xml must have an XML declaration at the beginning of the file and
a root element named geometry. Within the root element the user can define any
number of cells, surfaces, and lattices. Let us look at the following example:

.. code-block:: xml

    <?xml version="1.0"?>
    <geometry>
      <!-- This is a comment -->

      <surface>
        <id>1</id>
        <type>sphere</type>
        <coeffs>0.0 0.0 0.0 5.0</coeffs>
        <boundary>vacuum</boundary>
      <surface>

      <cell>
        <id>1</id>
        <universe>0</universe>
        <material>1</material>
        <surfaces>-1</surfaces>
      </cell>
    </geometry>

At the beginning of this file is a comment, denoted by a tag starting with
``<!--`` and ending with ``-->``. Comments, as well as any other type of input,
may span multiple lines. One convenient feature of the XML input format is that
sub-elements of the ``cell`` and ``surface`` elements can also be equivalently
expressed of attributes of the original element, e.g. the geometry file above
could be written as:

.. code-block:: xml

    <?xml version="1.0">
    <geometry>
      <!-- This is a comment -->

      <surface id="1" type="sphere" coeffs="0.0 0.0 0.0 5.0" boundary="vacuum" />
      <cell id="1" universe="0" material="1" surfaces="-1" />

    </geometry>

``<surface>`` Element
---------------------

Each ``<surface>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the surface.

    *Default*: None

  :type:
    The type of the surfaces. This can be "x-plane", "y-plane", "z-plane",
    "plane", "x-cylinder", "y-cylinder", "z-cylinder", or "sphere".

    *Default*: None

  :coeffs:
    The corresponding coefficients for the given type of surface. See below for
    a list a what coefficients to specify for a given surface

    *Default*: None

  :boundary:
    The boundary condition for the surface. This can be "transmission",
    "vacuum", or "reflective".

    *Default*: "transmission"

The following quadratic surfaces can be modeled:

  :x-plane:
    A plane perpendicular to the x axis, i.e. a surface of the form :math:`x -
    x_0 = 0`. The coefficients specified are ":math:`x_0`".

  :y-plane:
    A plane perpendicular to the y axis, i.e. a surface of the form :math:`y -
    y_0 = 0`. The coefficients specified are ":math:`y_0`".

  :z-plane:
    A plane perpendicular to the z axis, i.e. a surface of the form :math:`z -
    z_0 = 0`. The coefficients specified are ":math:`z_0`".

  :plane:
    An arbitrary plane of the form :math:`Ax + By + Cz = D`. The coefficients
    specified are ":math:`A \: B \: C \: D`".

  :x-cylinder:
    An infinite cylinder whose length is paralle to the x-axis. This is a
    quadratic surface of the form :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`y_0 \: z_0 \: R`".

  :y-cylinder:
    An infinite cylinder whose length is paralle to the y-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: z_0 \: R`".

  :z-cylinder:
    An infinite cylinder whose length is paralle to the z-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: y_0 \: R`".

  :sphere:
    A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 =
    R^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0 \: R`".

``<cell>`` Element
------------------

Each ``<cell>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the surface.

    *Default*: None

  :universe:
    The ``id`` of the universe that this cell is contained in.

    *Default*: 0

  :fill:
    The ``id`` of the universe that fills this cell.

    .. note:: If a fill is specified, no material should be given.

    *Default*: None

  :material:
    The ``id`` of the material that this cell contains.

    .. note:: If a material is specified, no fill should be given.

    *Default*: None

  :surfaces:
    A list of the ``ids`` for surfaces that bound this cell, e.g. if the cell
    is on the negative side of surface 3 and the positive side of surface 5, the
    bounding surfaces would be given as "-3 5".

    *Default*: None

``<lattice>`` Element
---------------------

The ``<lattice>`` can be used to represent repeating structures (e.g. fuel pins
in an assembly) or other geometry which naturally fits into a two-dimensional
structured mesh. Each cell within the lattice is filled with a specified
universe. A ``<lattice>`` accepts the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the surface.

  :type: A string indicating the arrangement of lattice cells. Currently, the
    only accepted option is "rectangular".

    *Default*: rectangular

  :dimension:
    Two integers representing the number of lattice cells in the x- and y-
    directions, respectively.

    *Default*: None

  :lower_left:
    The coordinates of the lower-left corner of the lattice.

    *Default*: None

  :width:
    The width of the lattice cell in the x- and y- directions.

    *Default*: None

  :universes:
    A list of the universe numbers that fill each cell of the lattice.

    *Default*: None

.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry

.. _quadratic surfaces: http://en.wikipedia.org/wiki/Quadric

----------------------------------------
Materials Specification -- materials.xml
----------------------------------------

``<material>`` Element
----------------------

Each ``material`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the material.

  :density:

    An element with attributes/sub-elements called ``value`` and ``units``. The
    ``value`` attribute is the numeric value of the density while the ``units``
    can be "g/cm3", "kg/m3", "atom/b-cm", "atom/cm3", or "sum". The "sum" unit
    indicates that the density should be calculated as the sum of the atom
    fractions for each nuclide in the material. This should not be used in
    conjunction with weight percents.

    *Default*: None

  :nuclide:
    An element with attributes/sub-elements called ``name``, ``xs``, and ``ao``
    or ``wo``. The ``name`` attribute is the name of the cross-section for a
    desired nuclide while the ``xs`` attribute is the cross-section
    identifier. Finally, the ``ao`` and ``wo`` attributes specify the atom or
    weight percent of that nuclide within the material, respectively. One
    example would be as follows::

      <nuclide name="H-1" xs="70c" ao="2.0" />
      <nuclide name="O-16" xs="70c" ao="1.0" />

    .. note:: If one nuclide is specified in atom percent, all others must also
              be given in atom percent. The same applies for weight percentages.

    *Default*: None

  :sab:
    Associates an S(a,b) table with the material. This element has
    attributes/sub-elements called ``name`` and ``xs``. The ``name`` attribute
    is the name of the S(a,b) table that should be associated with the material,
    and ``xs`` is the cross-section identifier for the table.

    *Default*: None

``<default_xs>`` Element
------------------------

In some circumstances, the cross-section identifier may be the same for many or
all nuclides in a given problem. In this case, rather than specifying the
``xs=...`` attribute on every nuclide, a ``<default_xs>`` element can be used to
set the default cross-section identifier for any nuclide without an identifier
explicitly listed. This element has no attributes and accepts a 3-letter string
that indicates the default cross-section identifier, e.g. "70c".

  *Default*: None

------------------------------------
Tallies Specification -- tallies.xml
------------------------------------

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

The two valid elements in the tallies.xml file are ``<tally>`` and ``<mesh>``.

``<tally>`` Element
-------------------

The ``<tally>`` element accepts the following sub-elements:

  :filters:
    A list of filters to specify what region of phase space should contribute to
    the tally. See below for full details on what filters are available.

  :scores:
    The desired responses to be accumulated. See below for full details on what
    responses can be tallied.

The following filters can be specified for a tally:

  :cell:
    A list of cells in which the tally should be accumulated.

  :cellborn:
    This filter allows the tally to be scored to only when particles were
    originally born in a specified cell.

  :surface:
    A list of surfaces for which the tally should be accumulated.

  :material:
    A list of materials for which the tally should be accumulated.

  :universe:
    A list of universes for which the tally should be accumulated.

  :energy:
    A monotonically increasing list of bounding **pre-collision** energies for a
    number of groups. For example, if the following energy filter is specified
    as ``<energy>0.0 1.0 20.0</energy>``, then two energy bins will be created,
    one with energies between 0 and 1 MeV and the other with energies between 1
    and 20 MeV.

  :energyout:
    A monotonically increasing list of bounding **post-collision** energies for
    a number of groups. For example, if the following energy filter is specified
    as ``<energy>0.0 1.0 20.0</energy>``, then two energy bins will be created,
    one with energies between 0 and 1 MeV and the other with energies between 1
    and 20 MeV.

  :mesh:
    The ``id`` of a structured mesh to be tallied over.

The following responses can be tallied.

  :flux:
    Total flux

  :total:
    Total reaction rate

  :scatter:
    Total scattering rate

  :nu-scatter:
    Total production of neutrons due to scattering. This accounts for
    multiplicity from (n,2n), (n,3n), and (n,4n) reactions and should be
    slightly higher than the scattering rate.

  :scatter-1:
    First scattering moment

  :scatter-2:
    Second scattering moment

  :scatter-3:
    Third scattering moment

  :absorption:
    Total absorption rate. This accounts for all reactions which do not produce
    secondary neutrons.

  :fission:
    Total fission rate

  :nu-fission:
    Total production of neutrons due to fission

``<mesh>`` Element
------------------

If a structured mesh is desired as a filter for a tally, it must be specified in
a separate element with the tag name ``<mesh>``. This element has the following
attributes/sub-elements:

  :type:
    The type of structured mesh. Valid options include "rectangular" and
    "hexagonal".

  :lower_left:
    The lower-left corner of the structured mesh. If only two coordinate are
    given, it is assumed that the mesh is an x-y mesh.

  :dimension:
    The number of mesh cells in each direction.

  :width:
    The width of mesh cells in each direction.

``<assume_separate>`` Element
-----------------------------

In cases where the user needs to specify many different tallies each of which
are spatially separate, this tag can be used to cut down on some of the tally
overhead. The effect of assuming all tallies are spatially separate is that once
one tally is scored to, the same event is assumed not to score to any other
tallies. This element should be followed by "yes" or "no"

  .. warning:: If used incorrectly, the assumption that all tallies are spatially
    separate can lead to incorrect results.

  *Default*: no

--------------------------------------------
Geometry Plotting Specification -- plots.xml
--------------------------------------------

A basic 2D plotting capability is available in OpenMC by creating a plots.xml
file and subsequently running with the command-line flag ``-plot``. The root
element of the plots.xml is simply ``<plots>`` and any number output figures can
be defined with ``<plot>`` sub-elements.

``<plot>`` Element
------------------

Each plot must contain a combination of the following attributes or sub-elements:

  :id:
    The unique ``id`` of the plot.

    *Default*: None - Required entry

  :filename:
    Filename for the output plot file.

    *Default*: "plot"

  :color:
    Keyword for plot coloring.  This can only be either ``cell`` or ``mat``,
    which colors regions by cells and materials, respectively.

    *Default*: ``cell``

  :origin:
    Specifies the (x,y,z) coordinate of the center of the plot.  Should be three
    floats separated by spaces.

    *Default*: None - Required entry

  :width:
    Specifies the width of the plot along each of the basis directions.  Should
    be two or three floats separated by spaces for 2D plots and 3D plots,
    respectively.

    *Default*: None - Required entry

  :type:
    Keyword for type of plot to be produced.  Currently only ``slice`` plots are
    implemented, which create 2D pixel maps saved in the PPM file format.  PPM
    files can be displayed in most viewers (e.g. the default Gnome viewer,
    IrfanView, etc.).

    .. note:: Since the PPM format is saved without any kind of compression,
              the resulting file sizes can be quite large.  Saving the image in
              the PNG format can often times reduce the file size by orders of
              magnitude without any loss of image quality.

    *Default*: "slice"

``<plot>`` elements of ``type`` "slice" also contain the following attributes or
sub-elements:

  :basis:
    Keyword specifying the plane of the plot for ``slice`` type plots.  Can be
    one of: "xy", "xz", "yz".

    *Default*: "xy"

  :pixels:
    Specifies the number of pixes to be used along each of the basis directions
    for "slice" plots. Should be two integers separated by spaces.

    .. warning:: The ``pixels`` input determines the output file size.  For the PPM
                 format, 10 million pixels will result in a file just under 30 MB in
                 size.

    .. warning:: If the aspect ratio defined in ``pixels`` does not match the aspect
              ratio defined in ``width`` the plot may appear stretched or squeezed.

    .. warning:: Geometry features along a basis direction smaller than ``width``/``pixels``
                 along that basis direction may not appear in the plot.

    *Default*: None - Required entry for "slice" plots

  :background:
    Specifies the RGB color of the regions where no OpenMC cell can be found. Should
    be three integers separated by spaces.

    *Default*: 0 0 0 (white)

  :col_spec:
    Any number of this optional tag may be included in each ``<plot>`` element, which can
    override the default random colors for cells or materials.  Each ``col_spec``
    element must contain ``id`` and ``rgb`` sub-elements.
  
    :id:
      Specifies the cell or material unique id for the color specification.

    :rgb:
      Specifies the custom color for the cell or material.  Should be 3 intergers separated
      by spaces.

    *Default*: None

  :mask:
    The special ``mask`` sub-element allows for the selective plotting of *only*
    user-specified cells or materials.  Only one ``mask`` element is allowed per ``plot``
    element, and it must contain as atributes or sub-elements a background masking color and
    a list of cells or materials to plot:

    :components:
      List of unique ``id`` numbers of the cells or materials to plot.  Should be any number
      of integers separated by spaces.

    :background:
      Color to apply to all cells or materials not in the ``components`` list of cells or
      materials to plot.  This overrides any ``col_spec`` color specifications.

    *Default*: None
