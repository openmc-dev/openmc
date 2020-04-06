.. _usersguide_geometry:

=================
Defining Geometry
=================

.. currentmodule:: openmc

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

- A plane perpendicular to the :math:`x` axis: :math:`x - x_0 = 0`
- A cylinder parallel to the :math:`z` axis: :math:`(x - x_0)^2 + (y -
  y_0)^2 - R^2 = 0`
- A sphere: :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0`

Defining a surface alone is not sufficient to specify a volume -- in order to
define an actual volume, one must reference the *half-space* of a surface. A
surface half-space is the region whose points satisfy a positive or negative
inequality of the surface equation. For example, for a sphere of radius one
centered at the origin, the surface equation is :math:`f(x,y,z) = x^2 + y^2 +
z^2 - 1 = 0`. Thus, we say that the negative half-space of the sphere, is
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
    | surface              | Dxy + Eyz + Fxz + Gx + Hy +  |                           |
    |                      | Jz + K = 0`                  |                           |
    +----------------------+------------------------------+---------------------------+

Each surface is characterized by several parameters. As one example, the
parameters for a sphere are the :math:`x,y,z` coordinates of the center of the
sphere and the radius of the sphere. All of these parameters can be set either
as optional keyword arguments to the class constructor or via attributes::

  sphere = openmc.Sphere(r=10.0)

  # This is equivalent
  sphere = openmc.Sphere()
  sphere.r = 10.0

Once a surface has been created, half-spaces can be obtained by applying the
unary ``-`` or ``+`` operators, corresponding to the negative and positive
half-spaces, respectively. For example::

   >>> sphere = openmc.Sphere(r=10.0)
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

While a bounding box can be determined for regions involving half-spaces of
spheres, cylinders, and axis-aligned planes, it generally cannot be determined
if the region involves cones, non-axis-aligned planes, or other exotic
second-order surfaces. For example, the :func:`openmc.model.hexagonal_prism`
function returns the interior region of a hexagonal prism; because it is bounded
by a :class:`openmc.Plane`, trying to get its bounding box won't work::

  >>> hex = openmc.model.hexagonal_prism()
  >>> hex.bounding_box
  (array([-0.8660254,       -inf,       -inf]),
   array([ 0.8660254,        inf,        inf]))

Boundary Conditions
-------------------

When a surface is created, by default particles that pass through the surface
will consider it to be transmissive, i.e., they pass through the surface
freely. If your model does not extend to infinity in all spatial dimensions, you
may want to specify different behavior for particles passing through a
surface. To specify a vacuum boundary condition, simply change the
:attr:`Surface.boundary_type` attribute to 'vacuum'::

   outer_surface = openmc.Sphere(r=100.0, boundary_type='vacuum')

   # This is equivalent
   outer_surface = openmc.Sphere(r=100.0)
   outer_surface.boundary_type = 'vacuum'

Reflective and periodic boundary conditions can be set with the strings
'reflective' and 'periodic'. Vacuum and reflective boundary conditions can be
applied to any type of surface. Periodic boundary conditions can be applied to
pairs of planar surfaces. For axis-aligned planes, matching periodic surfaces
can be determined automatically. For non-axis-aligned planes, it is necessary to
specify pairs explicitly using the :attr:`Surface.periodic_surface` attribute as
in the following example::

  p1 = openmc.Plane(a=0.3, b=5.0, d=1.0, boundary_type='periodic')
  p2 = openmc.Plane(a=0.3, b=5.0, d=-1.0, boundary_type='periodic')
  p1.periodic_surface = p2

Rotationally-periodic boundary conditions can be specified for a pair of
:class:`XPlane` and :class:`YPlane`; in that case, the
:attr:`Surface.periodic_surface` attribute must be specified manually as well.

.. caution:: When using rotationally-periodic boundary conditions, your geometry
             must be defined in the first quadrant, i.e., above the y-plane and
             to the right of the x-plane.

.. _usersguide_cells:

-----
Cells
-----

Once you have a material created and a region of space defined, you need to
define a *cell* that assigns the material to the region. Cells are created using
the :class:`openmc.Cell` class::

  fuel = openmc.Cell(fill=uo2, region=pellet)

  # This is equivalent
  fuel = openmc.Cell()
  fuel.fill = uo2
  fuel.region = pellet

In this example, an instance of :class:`openmc.Material` is assigned to the
:attr:`Cell.fill` attribute. One can also fill a cell with a :ref:`universe
<usersguide_universes>` or :ref:`lattice <usersguide_lattices>`. If you provide
no fill to a cell or assign a value of `None`, it will be treated as a "void"
cell with no material within. Particles are allowed to stream through the cell but
will undergo no collisions::

  # This cell will be filled with void on export to XML
  gap = openmc.Cell(region=pellet_gap)

The classes :class:`Halfspace`, :class:`Intersection`, :class:`Union`, and
:class:`Complement` and all instances of :class:`openmc.Region` and can be
assigned to the :attr:`Cell.region` attribute.

.. _usersguide_universes:

---------
Universes
---------

Similar to MCNP and Serpent, OpenMC is capable of using *universes*, collections
of cells that can be used as repeatable units of geometry. At a minimum, there
must be one "root" universe present in the model. To define a universe, an
instance of :class:`openmc.Universe` is created and then cells can be added
using the :meth:`Universe.add_cells` or :meth:`Universe.add_cell`
methods. Alternatively, a list of cells can be specified in the constructor::

   universe = openmc.Universe(cells=[cell1, cell2, cell3])

   # This is equivalent
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

Once a universe is constructed, it can actually be used to determine what cell
or material is found at a given location by using the :meth:`Universe.find`
method, which returns a list of universes, cells, and lattices which are
traversed to find a given point. The last element of that list would contain the
lowest-level cell at that location::

  >>> universe.find((0., 0., 0.))[-1]
  Cell
          ID             =    10000
          Name           =    cell 1
          Fill           =    Material 10000
          Region         =    -10000
          Rotation       =    None
          Temperature    =    None
          Translation    =    None

As you are building a geometry, it is also possible to display a plot of single
universe using the :meth:`Universe.plot` method. This method requires that you
have `matplotlib <http://matplotlib.org/>`_ installed.

.. _usersguide_lattices:

--------
Lattices
--------

Many particle transport models involve repeated structures that occur in a
regular pattern such as a rectangular or hexagonal lattice. In such a case, it
would be cumbersome to have to define the boundaries of each of the cells to be
filled with a universe. OpenMC provides a means to define lattice structures
through the :class:`openmc.RectLattice` and :class:`openmc.HexLattice` classes.

Rectangular Lattices
--------------------

A rectangular lattice defines a two-dimensional or three-dimensional array of
universes that are filled into rectangular prisms (lattice elements) each of
which has the same width, length, and height. To completely define a rectangular
lattice, one needs to specify

- The coordinates of the lower-left corner of the lattice
  (:attr:`RectLattice.lower_left`),
- The pitch of the lattice, i.e., the distance between the center of adjacent
  lattice elements (:attr:`RectLattice.pitch`),
- What universes should fill each lattice element
  (:attr:`RectLattice.universes`), and
- A universe that is used to fill any lattice position outside the well-defined
  portion of the lattice (:attr:`RectLattice.outer`).

For example, to create a 3x3 lattice centered at the origin in which each
lattice element is 5cm by 5cm and is filled by a universe ``u``, one could run::

  lattice = openmc.RectLattice()
  lattice.lower_left = (-7.5, -7.5)
  lattice.pitch = (5.0, 5.0)
  lattice.universes = [[u, u, u],
                       [u, u, u],
                       [u, u, u]]

Note that because this is a two-dimensional lattice, the lower-left coordinates
and pitch only need to specify the :math:`x,y` values. The order that the
universes appear is such that the first row corresponds to lattice elements with
the highest :math:`y` -value. Note that the :attr:`RectLattice.universes`
attribute expects a doubly-nested iterable of type :class:`openmc.Universe` ---
this can be normal Python lists, as shown above, or a NumPy array can be used as
well::

  lattice.universes = np.tile(u, (3, 3))

For a three-dimensional lattice, the :math:`x,y,z` coordinates of the lower-left
coordinate need to be given and the pitch should also give dimensions for all
three axes. For example, to make a 3x3x3 lattice where the bottom layer is
universe ``u``, the middle layer is universe ``q`` and the top layer is universe
``z`` would look like::

  lat3d = openmc.RectLattice()
  lat3d.lower_left = (-7.5, -7.5, -7.5)
  lat3d.pitch = (5.0, 5.0, 5.0)
  lat3d.universes = [
      [[u, u, u],
       [u, u, u],
       [u, u, u]],
      [[q, q, q],
       [q, q, q],
       [q, q, q]],
      [[z, z, z],
       [z, z, z]
       [z, z, z]]]

Again, using NumPy can make things easier::

  lat3d.universes = np.empty((3, 3, 3), dtype=openmc.Universe)
  lat3d.universes[0, ...] = u
  lat3d.universes[1, ...] = q
  lat3d.universes[2, ...] = z

Finally, it's possible to specify that lattice positions that aren't normally
without the bounds of the lattice be filled with an "outer" universe. This
allows one to create a truly infinite lattice if desired. An outer universe is
set with the :attr:`RectLattice.outer` attribute.

Hexagonal Lattices
------------------

OpenMC also allows creation of 2D and 3D hexagonal lattices. Creating a
hexagonal lattice is similar to creating a rectangular lattice with a few
differences:

- The center of the lattice must be specified (:attr:`HexLattice.center`).
- For a 2D hexagonal lattice, a single value for the pitch should be specified,
  although it still needs to appear in a list. For a 3D hexagonal lattice, the
  pitch in the radial and axial directions should be given.
- For a hexagonal lattice, the :attr:`HexLattice.universes` attribute cannot be
  given as a NumPy array for reasons explained below.
- As with rectangular lattices, the :attr:`HexLattice.outer` attribute will
  specify an outer universe.

For a 2D hexagonal lattice, the :attr:`HexLattice.universes` attribute should be
set to a two-dimensional list of universes filling each lattice element. Each
sub-list corresponds to one ring of universes and is ordered from the outermost
ring to the innermost ring. The universes within each sub-list are ordered from
the "top" (position with greatest y value) and proceed in a clockwise fashion
around the ring. The :meth:`HexLattice.show_indices` static method can be used
to help figure out how to place universes::

  >>> print(openmc.HexLattice.show_indices(3))
              (0, 0)
        (0,11)      (0, 1)
  (0,10)      (1, 0)      (0, 2)
        (1, 5)      (1, 1)
  (0, 9)      (2, 0)      (0, 3)
        (1, 4)      (1, 2)
  (0, 8)      (1, 3)      (0, 4)
        (0, 7)      (0, 5)
              (0, 6)


Note that by default, hexagonal lattices are positioned such that each lattice
element has two faces that are parallel to the :math:`y` axis. As one example,
to create a three-ring lattice centered at the origin with a pitch of 10 cm
where all the lattice elements centered along the :math:`y` axis are filled with
universe ``u`` and the remainder are filled with universe ``q``, the following
code would work::

  hexlat = openmc.HexLattice()
  hexlat.center = (0, 0)
  hexlat.pitch = [10]

  outer_ring = [u, q, q, q, q, q, u, q, q, q, q, q]
  middle_ring = [u, q, q, u, q, q]
  inner_ring = [u]
  hexlat.universes = [outer_ring, middle_ring, inner_ring]

If you need to create a hexagonal boundary (composed of six planar surfaces) for
a hexagonal lattice, :func:`openmc.model.hexagonal_prism` can be used.

.. _usersguide_geom_export:

--------------------------
Exporting a Geometry Model
--------------------------

Once you have finished building your geometry by creating surfaces, cell, and,
if needed, lattices, the last step is to create an instance of
:class:`openmc.Geometry` and export it to an XML file that the
:ref:`scripts_openmc` executable can read using the
:meth:`Geometry.export_to_xml` method. This can be done as follows::

   geom = openmc.Geometry(root_univ)
   geom.export_to_xml()

   # This is equivalent
   geom = openmc.Geometry()
   geom.root_universe = root_univ
   geom.export_to_xml()

Note that it's not strictly required to manually create a root universe. You can
also pass a list of cells to the :class:`openmc.Geometry` constructor and it
will handle creating the unverse::

   geom = openmc.Geometry([cell1, cell2, cell3])
   geom.export_to_xml()

.. _constructive solid geometry: https://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _quadratic surfaces: https://en.wikipedia.org/wiki/Quadric

--------------------------
Using CAD-based Geometry
--------------------------

OpenMC relies on the Direct Accelerated Geometry Monte Carlo toolkit (`DAGMC
<https://svalinn.github.io/DAGMC/>`_) to represent CAD-based geometry in a
surface mesh format. A DAGMC run can be enabled in OpenMC by setting the
``dagmc`` property to ``True`` in the model Settings either via the Python
:class:`openmc.settings` Python class::

  settings = openmc.Settings()
  settings.dagmc = True

or in the :ref:`settings.xml <io_settings>` file::

  <dagmc>true</dagmc>

With ``dagmc`` set to true, OpenMC will load the DAGMC model (from a local file
named ``dagmc.h5m``) when initializing a simulation. If a `geometry.xml
<../io_formats/geometry.html>`_ is present as well, it will be ignored.

  **Note:** DAGMC geometries used in OpenMC are currently required to be clean,
  meaning that all surfaces have been `imprinted and merged
  <https://svalinn.github.io/DAGMC/usersguide/trelis_workflow.html>`_
  successfully and that the model is `watertight
  <https://svalinn.github.io/DAGMC/usersguide/tools.html#make-watertight>`_. Future
  implementations of DAGMC geometry will support small volume overlaps and
  un-merged surfaces.

-------------------------
Calculating Atoms Content
-------------------------

If the total volume occupied by all instances of a cell in the geometry is known
by the user, it is possible to assign this volume to a cell without performing a
:ref:`stochastic volume <usersguide_volume>` calculation::

  from uncertainties import ufloat

  # Set known total volume in [cc]
  cell = openmc.Cell()
  cell.volume = 17.0

  # Set volume if it is known with some uncertainty
  cell.volume = ufloat(17.0, 0.1)

Once a volume is set, and a cell is filled with a material or distributed
materials, it is possible to use the :func:`~openmc.Cell.atoms` method to obtain
a dictionary of nuclides and their total number of atoms in all instances
of a cell (e.g. ``{'H1': 1.0e22, 'O16': 0.5e22, ...}``)::

  cell = openmc.Cell(fill = u02)
  cell.volume = 17.0

  O16_atoms = cell.atoms['O16']
