.. _usersguide_geometry:

=================
Defining Geometry
=================

--------------------
Surfaces and Regions
--------------------

The geometry of a model in OpenMC is defined using `constructive solid
geometry`_ (CSG), also sometimes referred to as combinatorial geometry. CSG
allows a user to create complex regions using Boolean operators (intersection,
union, and complement) on simpler regions. In order to define a region that we
can assign to a cell, we must first define surfaces which bound the region. A
surface is a locus of zeros of a function of Cartesian coordinates
:math:`x,y,z`, e.g.

- A plane perpendicular to the :math:`x` axis: :math:`x − x_0 = 0`
- A cylinder perpendicular to the :math:`z` axis: :math:`(x − x_0)^2 + (y −
  y_0)^2 − R^2 = 0`
- A sphere: :math:`(x − x_0)^2 + (y − y_0)^2 + (z − z_0)^2 − R^2 = 0`

Defining a surface alone is not sufficient to specify a volume -- in order to
define an actual volume, one must reference the *half-space* of a surface. A
surface half-space is the region whose points satisfy a positive of negative
inequality of the surface equation. For example, for a sphere of radius one
centered at the origin, the surface equation is :math:`f(x,y,z) = x^2 + y^2 +
z^2 − 1 = 0`. Thus, we say that the negative half-space of the sphere, is
defined as the collection of points satisfying :math:`f(x,y,z) < 0`, which one
can reason is the inside of the sphere. Conversely, the positive half-space of
the sphere would correspond to all points outside of the sphere, satisfying
:math:`f(x,y,z) > 0`.

In the Python API, surfaces are created via subclasses of
:class:`openmc.Surface`. The available surface types and their corresponding
classes are listed in the following table.

.. table:: Surface types available in OpenMC.

    +----------------------+------------------------------+---------------------------+
    | Surface              | Equation                     | Class                     |
    +======================+==============================+===========================+
    | Plane perpendicular  | :math:`x - x_0 = 0`          | :class:`openmc.XPlane`    |
    | to :math:`x`-axis    |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Plane perpendicular  | :math:`y - y_0 = 0`          | :class:`openmc.YPlane`    |
    | to :math:`y`-axis    |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Plane perpendicular  | :math:`z - z_0 = 0`          | :class:`openmc.ZPlane`    |
    | to :math:`z`-axis    |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Arbitrary plane      | :math:`Ax + By + Cz = D`     | :class:`openmc.Plane`     |
    +----------------------+------------------------------+---------------------------+
    | Infinite cylinder    | :math:`(y-y_0)^2 + (z-z_0)^2 | :class:`openmc.XCylinder` |
    | parallel to          | - R^2 = 0`                   |                           |
    | :math:`x`-axis       |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Infinite cylinder    | :math:`(x-x_0)^2 + (z-z_0)^2 | :class:`openmc.YCylinder` |
    | parallel to          | - R^2 = 0`                   |                           |
    | :math:`y`-axis       |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Infinite cylinder    | :math:`(x-x_0)^2 + (y-y_0)^2 | :class:`openmc.ZCylinder` |
    | parallel to          | - R^2 = 0`                   |                           |
    | :math:`z`-axis       |                              |                           |
    +----------------------+------------------------------+---------------------------+
    | Sphere               | :math:`(x-x_0)^2 + (y-y_0)^2 | :class:`openmc.Sphere`    |
    |                      | + (z-z_0)^2 - R^2 = 0`       |                           |
    +----------------------+------------------------------+---------------------------+
    | Cone parallel to the | :math:`(y-y_0)^2 + (z-z_0)^2 | :class:`openmc.XCone`     |
    | :math:`x`-axis       | - R^2(x-x_0)^2 = 0`          |                           |
    +----------------------+------------------------------+---------------------------+
    | Cone parallel to the | :math:`(x-x_0)^2 + (z-z_0)^2 | :class:`openmc.YCone`     |
    | :math:`y`-axis       | - R^2(y-y_0)^2 = 0`          |                           |
    +----------------------+------------------------------+---------------------------+
    | Cone parallel to the | :math:`(x-x_0)^2 + (y-y_0)^2 | :class:`openmc.ZCone`     |
    | :math:`z`-axis       | - R^2(z-z_0)^2 = 0`          |                           |
    +----------------------+------------------------------+---------------------------+
    | General quadric      | :math:`Ax^2 + By^2 + Cz^2 +  |  :class:`openmc.Quadric`  |
    | surface              | Dxy + Eyz + Fxz \\+Gx + Hy + |                           |
    |                      | Jz + K = 0`                  |                           |
    +----------------------+------------------------------+---------------------------+

Each surface is characterized by several parameters. As one example, the
parameters for a sphere are the :math:`x,y,z` coordinates of the center of the
sphere and the radius of the sphere. All of these parameters can be set either
as optional keyword arguments to the class constructor or via attributes::

  sphere = openmc.Sphere(R=10.0)

  # ..or..
  sphere = openmc.Sphere()
  sphere.r = 10.0

Once a surface has been created, half-spaces can be obtained by applying the
unary ``-`` or ``+`` operators, corresponding to the negative and positive
half-spaces, respectively. For example::

   >>> sphere = openmc.Sphere(R=10.0)
   >>> inside_sphere = -sphere
   >>> outside_sphere = +sphere
   >>> type(inside_sphere)
   <class 'openmc.surface.Halfspace'>

Instances of :class:`openmc.Halfspace` can be combined together using the
Boolean operators ``&`` (intersection), ``|`` (union), and ``~`` (complement)::

  >>> inside_sphere = -openmc.Sphere()
  >>> above_plane = +openmc.ZPlane()
  >>> northern_hemisphere = inside_sphere & above_plane
  >>> type(northern_hemisphere)
  <class 'openmc.region.Intersection'>

For many regions, a bounding-box can be determined automatically::

  >>> northern_hemisphere.bounding_box
  (array([-1., -1., 0.]), array([1., 1., 1.]))

Boundary Conditions
-------------------

When a surface is created, by default particles that pass through the surface
will consider it to be transmissive, i.e., they pass through the surface
freely. If your model does not extend to infinity in all spatial dimensions, you
may want to specify different behavior for particles passing through a
surface. To specify a vacuum boundary condition, simply change the
:attr:`Surface.boundary_type` attribute to 'vacuum'::

   outer_surface = openmc.Sphere(R=100.0, boundary_type='vacuum')

   # ..or..
   outer_surface = openmc.Sphere(R=100.0)
   outer_surface.boundary_type = 'vacuum'

Reflective and periodic boundary conditions can be set with the strings
'reflective' and 'periodic'. Vacuum and reflective boundary conditions can be
applied to any type of surface. Periodic boundary conditions can only be applied
to pairs of axis-aligned planar surfaces.

-----
Cells
-----

Once you have a material created and a region of space defined, you need to
define a *cell* that assigns the material to the region. Cells are created using
the :class:`openmc.Cell` class::

  fuel = openmc.Cell(fill=uo2, region=pellet)

  # ..or..
  fuel = openmc.Cell()
  fuel.fill = uo2
  fuel.region = pellet

The classes :class:`Halfspace`, :class:`Intersection`, :class:`Union`, and
:class:`Complement` and all instances of :class:`openmc.Region` and can be
assigned to the :attr:`Cell.region` attribute.

---------
Universes
---------

Similar to MCNP and Serpent, OpenMC is capable of using *universes*, collections
of cells that can be used as repeatable units of geometry. At a minimum, there
must be one "root" universe present in the model. To create a universe, the
:class:`openmc.Universe` is used::

   universe = openmc.Universe(cells=[cell1, cell2, cell3])

   # ..or..
   universe = openmc.Universe()
   universe.add_cells([cell1, cell2])
   universe.add_cell(cell3)

Universes are generally used in three ways:

1. To be assigned to a :class:`Geometry` object (see
   :ref:`usersguide_geom_export`),
2. To be assigned as the fill for a cell via the :attr:`Cell.fill` attribute,
   and
3. To be used in a regular arrangement of universes in a :ref:`lattice
   <usersguide_lattices>`.

.. _usersguide_lattices:

--------
Lattices
--------


------------------
Hexagonal Lattices
------------------


.. _usersguide_geom_export:

--------------------------
Exporting a Geometry Model
--------------------------

.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _quadratic surfaces: http://en.wikipedia.org/wiki/Quadric
