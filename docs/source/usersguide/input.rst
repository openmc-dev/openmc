.. _usersguide_input:

========================
Creating XML Input Files
========================

Unlike many other Monte Carlo codes which use an arbitrary-format ASCII file
with "cards" to specify a particular geometry, materials, and associated run
settings, the input files for OpenMC are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.

Anyone who has ever seen webpages written in HTML will be familiar with the
structure of XML whereby "tags" enclosed in angle brackets denote that a
particular piece of data will follow. Let us examine the follow example::

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
for the geometry, materials, and settings. Additionally, an optional tallies XML
file specifies physical quantities to be tallied. OpenMC expects that these
files are called:

* ``geometry.xml``
* ``materials.xml``
* ``setings.xml``
* ``tallies.xml``

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
number of cells, surfaces, and lattices. Let us look at the following example::

    <?xml version="1.0">
    <geometry>
      <!-- This is a comment -->

      <surface>
        <uid>1</uid>
        <type>sphere</type>
        <coeffs>0.0 0.0 0.0 5.0</coeffs>
        <boundary>vacuum</boundary>
      <surface>

      <cell>
        <uid>1</uid>
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
could be written as::

    <?xml version="1.0">
    <geometry>
      <!-- This is a comment -->

      <surface uid="1" type="sphere" coeffs="0.0 0.0 0.0 5.0" boundary="vacuum" />
      <cell uid="1" universe="0" material="1" surfaces="-1" />

    </geometry>

.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry

.. _quadratic surfaces: http://en.wikipedia.org/wiki/Quadric

Types of surfaces:

``x-plane``
  A plane perpendicular to the x axis, i.e. a surface of the form :math:`x - x_0
  = 0`. The coefficients specified are ":math:`x_0`".

``y-plane``
  A plane perpendicular to the y axis, i.e. a surface of the form :math:`y - y_0
  = 0`. The coefficients specified are ":math:`y_0`".

``z-plane``
  A plane perpendicular to the z axis, i.e. a surface of the form :math:`z - z_0
  = 0`. The coefficients specified are ":math:`z_0`".

``plane``
  An arbitrary plane of the form :math:`Ax + By + Cz = D`. The coefficients
  specified are ":math:`A \: B \: C \: D`".

``x-cylinder``
  An infinite cylinder whose length is paralle to the x-axis. This is a
  quadratic surface of the form :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`. The
  coefficients specified are ":math:`y_0 \: z_0 \: R`".

``y-cylinder``
  An infinite cylinder whose length is paralle to the y-axis. This is a
  quadratic surface of the form :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`. The
  coefficients specified are ":math:`x_0 \: z_0 \: R`".

``z-cylinder``
  An infinite cylinder whose length is paralle to the z-axis. This is a
  quadratic surface of the form :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`. The
  coefficients specified are ":math:`x_0 \: y_0 \: R`".

``sphere``
  A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 =
  R^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0 \: R`".

----------------------------------------
Materials Specification -- materials.xml
----------------------------------------

--------------------------------------
Settings Specification -- settings.xml
--------------------------------------

------------------------------------
Tallies Specification -- tallies.xml
------------------------------------

