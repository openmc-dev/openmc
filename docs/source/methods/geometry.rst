.. _methods_geometry:

========
Geometry
========

---------------------------
Constructive Solid Geometry
---------------------------

OpenMC uses a technique known as `constructive solid geometry`_ (CSG) to build
arbitrarily complex three-dimensional models in Euclidean space. In a CSG model,
every unique object is described as the union, intersection, or difference of
half-spaces created by bounding `surfaces`_. Every surface divides all of space
into exactly two half-spaces. We can mathematically define a surface as a
collection of points that satisfy an equation of the form :math:`f(x,y,z) = 0`
where :math:`f(x,y,z)` is a given function. The region for which :math:`f(x,y,z)
< 0` can be called the negative half-space (or simply the "negative side") and
the region for which :math:`f(x,y,z) > 0` can be called the positive half-space.

Let us take the example of a sphere centered at the point :math:`(x_0,y_0,z_0)`
with radius :math:`R`. One would normally write the equation of the sphere as

.. math::

    (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 = R^2

By subtracting the right-hand term from both sides of the equation, we can then
write the surface equation:

.. math::

    f(x,y,z) = (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0

One can confirm that any point inside this sphere will correspond to
:math:`f(x,y,z) < 0` and any point outside the sphere will correspond to
:math:`f(x,y,z) > 0`.

In OpenMC, every surface defined by the user is assigned an integer to uniquely
identify it. We can then refer to either of the two half-spaces created by a
surface by a combination of the unique ID of the surface and a positive/negative
sign. For example, to refer to the negative half-space of a sphere (the volume
inside the sphere) with unique ID 35, the reference would be -35. These
references to half-spaces are used in created regions in space of homogeneous
material, known as "cells".


.. figure:: ../../img/halfspace.svg
   :align: center
   :figclass: align-center
   
.. figure:: ../../img/union.svg
   :align: center
   :figclass: align-center

In OpenMC, any second-order surface of the form

.. math::

    f(x,y,z) = Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0

can be modeled in OpenMC. For example, the equation for a sphere centered at
:math:`(\bar{x},\bar{y},\bar{z})` and of radius :math:`R` can be written as

Universes
---------

Lattices
--------

------------------------------------------
Computing the Distance to Nearest Boundary
------------------------------------------

One of the most basic algorithms in any Monte Carlo code is determining the
distance to the nearest surface within a cell. Since each cell is defined by
the surfaces that bound it, if we compute the distance to all surfaces bounding
a cell, we can determine the nearest one.

With the possibility of a particle having coordinates on multiple levels
(universes) in a geometry, we must exercise care when calculating the distance
to the nearest surface. Each different level of geometry has a set of boundaries
with which the particle's direction of travel may intersect. Thus, it is
necessary to check the distance to the surfaces bounding the cell in each
level. This should be done starting the highest (most global) level going down
to the lowest (most local) level. That ensures that if two surfaces on different
levels are coincident, by default the one on the higher level will be selected
as the nearest surface.

The following procedure is used to calculate the distance to each bounding
surface. Suppose we have a particle at :math:`(x,y,z)` traveling in the
direction :math:`u,v,w`. To find the distance :math:`d` to a surface
:math:`f(x,y,z) = 0`, we need to solve the equation:

.. math::
    :label: dist-to-boundary-1

    f(x + du, y + dv, z + dw) = 0

If no solutions to equation :eq:`dist-to-boundary-1` exists or the only
solutions are complex, then the particle's direction of travel will not
intersect the surface. If the solution to equation :eq:`dist-to-boundary-1` is
negative, this means that the surface is "behind" the particle, i.e. if the
particle continues traveling in its current direction, it will not hit the
surface. The complete derivation for different types of surfaces used in OpenMC
will be presented in the following sections.

Once a distance has been computed to a boundary, we need to check if it is
closer than previously-computed distances to surfaces. Unfortunately, we cannot
just use the minimum function because some distances may be almost identical but
still different due to the use of floating-point arithmetic. Consequently, we
should first check for floating-point equality of the current distance
calculated and the minimum found thus far. This is done by checking if

.. math::
    :label: fp-distance

    \frac{| d - d_{min} |}{d_{min}} < \epsilon

where :math:`d` is the distance to a surface just calculated, :math:`d_{min}` is
the minimum distance found thus far, and :math:`\epsilon` is a small number. In
OpenMC, this parameter is set to :math:`\epsilon = 10^{-14}` since all floating
calculations are done on 8-byte floating point numbers.

Although they are not explicitly defined, it is also necessary to check the
distance to surfaces representing lattice boundaries if a lattice exists on a
given level.

Plane Perpendicular to an Axis
------------------------------

The equation for a plane perpendicular to, for example, the x-axis is simply
:math:`x - x_0 = 0`. As such, we need to solve :math:`x + du - x_0 = 0`. The
solution for the distance is

.. math::
    :label: dist-xplane

    d = \frac{x_0 - x}{u}

Note that if the particle's direction of flight is parallel to the x-axis,
i.e. :math:`u = 0`, the distance to the surface will be infinity. While the
example here was for a plane perpendicular to the x-axis, the same formula can
be applied for the surfaces :math:`y = y_0` and :math:`z = z_0`.

Generic Plane
-------------

The equation for a generic plane is :math:`Ax + By + Cz = D`. Thus, we need to
solve the equation :math:`A(x + du) + B(y + dv) + C(z + dw) = D`. The solution
to this equation for the distance is

.. math::
    :label: dist-plane

    d = \frac{D - Ax - By - Cz}{Au + Bv + Cw}

Again, we need to check whether the denominator is zero. If so, this means that
the particle's direction of flight is parallel to the plane and it will
therefore never hit the plane.

Cylinder Parallel to an Axis
----------------------------

The equation for a cylinder parallel to, for example, the x-axis is :math:`(y -
y_0)^2 + (z - z_0)^2 = R^2`. Thus, we need to solve :math:`(y + dv - y_0)^2 +
(z + dw - z_0)^2 = R^2`. Let us define :math:`\bar{y} = y - y_0` and
:math:`\bar{z} = z - z_0`. We then have

.. math::
    :label: dist-xcylinder-1

    (\bar{y} + dv)^2 + (\bar{z} + dw)^2 = R^2

Expanding equation :eq:`dist-xcylinder-1` and rearranging terms, we obtain

.. math::
    :label: dist-xcylinder-2

    (v^2 + w^2) d^2 + 2 (\bar{y}v + \bar{z}w) d + (\bar{y}^2 + \bar{z}^2 - R^2)
    = 0

This is a quadratic equation for :math:`d`. To simplify notation, let us define
:math:`a = v^2 + w^2`, :math:`k = \bar{y}v + \bar{z}w`, and :math:`c =
\bar{y}^2 + \bar{z}^2 - R^2`. Thus, the distance is just the solution to
:math:`ad^2 + 2kd + c = 0`:

.. math::
    :label: dist-xcylinder-3

    d = \frac{-k \pm \sqrt{k^2 - ac}}{a}

A few conditions must be checked for. If :math:`a = 0`, this means the particle
is parallel to the cylinder and will thus never intersect it. Also, if
:math:`k^2 - ac < 0`, this means that both solutions to the quadratic are
complex. In physical terms, this means that the ray along which the particle is
traveling does not make any intersections with the cylinder.

If we do have intersections and :math:`c < 0`, this means that the particle is
inside the cylinder. Thus, one solution should be positive and one should be
negative. Clearly, the positive distance will occur when the sign on the
square root of the discriminant is positive since :math:`a > 0`.

If we have intersections and :math:`c > 0` this means that the particle is
outside the cylinder. Thus, the solutions to the quadratic are either both
positive or both negative. If they are both positive, the smaller (closer) one
will be the solution with a negative sign on the square root of the
discriminant.

The same equations and logic here can be used for cylinders that are parallel to
the y- or z-axis with appropriate substitution of constants.

Sphere
------

The equation for a sphere is :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 =
R^2`. Thus, we need to solve the equation

.. math::
    :label: dist-sphere-1

    (x + du - x_0)^2 + (y + dv - y_0)^2 + (z + dw - z_0)^2 = R^2

Let us define :math:`\bar{x} = x - x_0`, :math:`\bar{y} = y - y_0`, and
:math:`\bar{z} = z - z_0`. We then have

.. math::
    :label: dist-sphere-2

    (\bar{x} + du)^2 + (\bar{y} + dv)^2 + (\bar{z} - dw)^2 = R^2

Expanding equation :eq:`dist-sphere-2` and rearranging terms, we obtain

.. math::
    :label: dist-sphere-3

    d^2 + 2 (\bar{x}u + \bar{y}v + \bar{z}w) d + (\bar{x}^2 + \bar{y}^2 +
    \bar{z}^2 - R^2) = 0

This is a quadratic equation for :math:`d`. To simplify notation, let us define
:math:`k = \bar{x}u + \bar{y}v + \bar{z}w` and :math:`c = \bar{x}^2 +
\bar{y}^2 + \bar{z}^2 - R^2`. Thus, the distance is just the solution to
:math:`d^2 + 2kd + c = 0`:

.. math::
    :label: dist-sphere-4

    d = -k \pm \sqrt{k^2 - c}

If the discriminant :math:`k^2 - c < 0`, this means that both solutions to the
quadratic are complex. In physical terms, this means that the ray along which
the particle is traveling does not make any intersections with the sphere.

If we do have intersections and :math:`c < 0`, this means that the particle is
inside the sphere. Thus, one solution should be positive and one should be
negative. The positive distance will occur when the sign on the square root of
the discriminant is positive. If we have intersections but :math:`c > 0` this
means that the particle is outside the sphere. The solutions to the quadratic
will then be either both positive or both negative. If they are both positive,
the smaller (closer) one will be the solution with a negative sign on the square
root of the discriminant.

.. _find-cell:

----------------------------
Finding a Cell Given a Point
----------------------------

Another basic algorithm is to determine which cell contains a given point in the
global coordinate system, i.e. if the particle's position is :math:`(x,y,z)`,
what cell is it currently in. This is done in the following manner in
OpenMC. With the possibility of multiple levels of coordinates, we must perform
a recursive search for the cell. First, we start in the highest (most global)
universe which we call the base universe and do a loop over each cell within
that universe. For each cell, we check whether the specified point is inside the
cell using the algorithm described in :ref:`cell-contains`. If the cell is
filled with a normal material, the search is done and we have identified the
cell containing the point. If the cell is filled with another universe, we then
search all cells within that universe to see if any of them contain the
specified point. If the cell is filled with a lattice, the position within the
lattice is determined, and then whatever universe fills that lattice position is
recursively searched. The search ends once a cell containing a normal material
is found that contains the specified point.

.. _cell-contains:

----------------------------------------
Determining if a Coordinate is in a Cell
----------------------------------------

One aspect of being able to determine what cell a particle is in is determining
if a particle's coordinates lie within a given cell. The current geometry
implementation in OpenMC limits all cells to being simple cells, i.e. they are
defined only with intersection of half-spaces and not unions, differences,
etc. This makes the job of determining if a point is in a cell quite simple.

The algorithm for determining if a cell contains a point is as follows. For each
surface that bounds a cell, we determine the particle's sense with respect to
the surface. As explained earlier, if we have a point :math:`(x_0,y_0,z_0)` and
a surface :math:`f(x,y,z) = 0`, the point is said to have negative sense if
:math:`f(x_0,y_0,z_0) < 0` and positive sense if :math:`f(x_0,y_0,z_0) > 0`. If
for all surfaces, the sense of the particle with respect to the surface matches
the specified sense that defines the half-space within the cell, then the point
is inside the cell.

Let us illustrate this idea with a concept. Let's say we have a cell defined as

.. code-block:: xml

    <cell id="1" surfaces="-1 2 -3" />

    <surface id="1" type="sphere"  coeffs="0 0 0 10" />
    <surface id="2" type="x-plane" coeffs="-3" />
    <surface id="3" type="y-plane" coeffs="2" />

This means that the cell is defined as the intersection of the negative half
space of a sphere, the positive half-space of an x-plane, and the negative
half-space of a y-plane. Said another way, any point inside this cell must
satisfy the following equations

.. math::
    :label: cell-contains-example

    x^2 + y^2 + z^2 - 10^2 < 0 \\
    x - (-3) > 0 \\
    x - 2 < 0

So in order to determine if a point is inside the cell, we would plug its
coordinates into equation :eq:`cell-contains-example` and if the inequalities
are satisfied, than the point is indeed inside the cell.

--------------------------
Handling Surface Crossings
--------------------------

A particle will cross a surface if the distance to the nearest surface is closer
than the distance sampled to the next collision. A number of things happen when
a particle hits a surface. First, we need to check if a non-transmissive
boundary condition has been applied to the surface. If a vacuum boundary
condition has been applied, the particle is killed and any surface current
tallies are scored to as needed. If a reflective boundary condition has been
applied to the surface, surface current tallies are scored to and then the
particle's direction is changed according to the procedure in :ref:`reflection`.

Next, we need to determine what cell is beyond the surface in the direction of
travel of the particle so that we can evaluate cross sections based on its
material properties. At initialization, a list of neighboring cells is created
for each surface in the problem as described in :ref:`neighbor-lists`. The
algorithm outlined in :ref:`find-cell` is used to find a cell containing the
particle except rather than searching all cells in the base universe, only the
list of neighboring cells is searched. If this search is unsuccessful, then a
search is done over every cell in the base universe.

.. _neighbor-lists:

-----------------------
Building Neighbor Lists
-----------------------

After the geometry has been loaded and stored in memory from an input file,
OpenMC builts a list for each surface containing any cells that contain the
surface in their specification in order to speed up processing of surface
crossings. The algorithm to build these lists is as follows. First, we loop over
all cells in the geometry and count up how many times each surface appears in a
specification as bounding a negative half-space and bounding a positive
half-space. Two arrays are then allocated for each surface, one that lists each
cell that contains the negative half-space of the surface and one that lists
each cell that contains the positive half-space of the surface. Another loop is
performed over all cells and the neighbor lists are populated for each surface.

.. _reflection:

------------------------------
Reflective Boundary Conditions
------------------------------

If the velocity of a particle is :math:`\mathbf{v}` and it crosses a surface of
the form :math:`f(x,y,z) = 0` with a reflective boundary condition, it can be
shown based on geometric arguments that the velocity vector will then become

.. math::
    :label: reflection-v

    \mathbf{v'} = \mathbf{v} - 2 (\mathbf{v} \cdot \hat{\mathbf{n}})
    \hat{\mathbf{n}}

where :math:`\hat{\mathbf{n}}` is a unit vector normal to the surface at the
point of the surface crossing. The rationale for this can be understood by
noting that :math:`(\mathbf{v} \cdot \hat{\mathbf{n}}) \hat{\mathbf{n}}` is the
projection of the velocity vector onto the normal vector. By subtracting two
times this projection, the velocity is reflected with respect to the surface
normal. Since the velocity of the particle will not change as it undergoes
reflection, we can work with the direction of the particle instead, simplifying
equation :eq:`reflection-v` to

.. math::
    :label: reflection-omega

    \mathbf{\Omega'} = \mathbf{\Omega} - 2 (\mathbf{\Omega} \cdot
    \hat{\mathbf{n}}) \hat{\mathbf{n}}


The direction of the surface normal will be the gradient to the surface at the
point of crossing, i.e. :math:`\mathbf{n} = \nabla f(x,y,z)`. Substituting this
into equation :eq:`reflection-omega`, we get

.. math::
    :label: reflection-omega-2

    \mathbf{\Omega'} = \mathbf{\Omega} - \frac{2 ( \mathbf{\Omega} \cdot \nabla
    f )}{|| \nabla f ||^2} \nabla f


If we write the initial and final directions in terms of their vector
components, :math:`\mathbf{\Omega} = (u,v,w)` and :math:`\mathbf{\Omega'} = (u',
v', w')`, this allows us to represent equation :eq:`reflection-omega` as a
series of equations:

.. math::
    :label: reflection-system

    u' = u - \frac{2 ( \mathbf{\Omega} \cdot \nabla f )}{|| \nabla f ||^2}
    \frac{\partial f}{\partial x} \\

    v' = v - \frac{2 ( \mathbf{\Omega} \cdot \nabla f )}{|| \nabla f ||^2}
    \frac{\partial f}{\partial y} \\

    w' = w - \frac{2 ( \mathbf{\Omega} \cdot \nabla f )}{|| \nabla f ||^2}
    \frac{\partial f}{\partial z}

We can now use this form to develop rules for how to transform a particle's
direction for different types of surfaces.

Plane Perpendicular to an Axis
------------------------------

For a plane that is perpendicular to an axis, the rule for reflection is almost
so simple that no derivation is needed at all. Nevertheless, we will proceed
with the derivation to confirm that the rules of geometry agree with our
intuition. The gradient of the surface :math:`f(x,y,z) = x - x_0 = 0` is simply
:math:`\nabla f = (1, 0, 0)`. Note that this vector is already normalized,
i.e. :math:`|| \nabla f || = 1`. The second two equations in
:eq:`reflection-system` tell us that :math:`v` and :math:`w` do not change and
the first tell us that

.. math::
    :label: reflection-xplane

    u' = u - 2u = -u

We see that reflection for a plane perpendicular to an axis only entails
negating the directional cosine for that axis.

Generic Plane
-------------

A generic plane has the form :math:`f(x,y,z) = Ax + By + Cz - D = 0`. Thus, the
gradient to the surface is simply :math:`\nabla f = (A,B,C)` whose norm squared
is :math:`A^2 + B^2 + C^2`. This implies that

.. math::
    :label: reflection-plane-constant

    \frac{2 (\mathbf{\Omega} \cdot \nabla f)}{|| \nabla f ||^2} = \frac{2(Au +
    Bv + Cw)}{A^2 + B^2 + C^2}

Substituting equation :eq:`reflection-plane-constant` into equation
:eq:`reflection-system` gives us the form of the solution. For example, the
x-component of the reflected direction will be

.. math::
    :label: reflection-plane

    u' = u - \frac{2A(Au + Bv + Cw)}{A^2 + B^2 + C^2}


Cylinder Parallel to an Axis
----------------------------

A cylinder parallel to, for example, the x-axis has the form :math:`f(x,y,z) =
(y - y_0)^2 + (z - z_0)^2 - R^2 = 0`. Thus, the gradient to the surface is

.. math:: :label: reflection-cylinder-grad

    \nabla f = 2 \left ( \begin{array}{c} 0 \\ y - y_0 \\ z - z_0 \end{array}
    \right ) = 2 \left ( \begin{array}{c} 0 \\ \bar{y} \\ \bar{z} \end{array}
    \right )

where we have introduced the constants :math:`\bar{y}` and
:math:`\bar{z}`. Taking the square of the norm of the gradient, we find that

.. math:: :label: reflection-cylinder-norm

    || \nabla f ||^2 = 4 \bar{y}^2 + 4 \bar{z}^2 = 4 R^2

This imples that

.. math:: :label: reflection-cylinder-constant

    \frac{2 (\mathbf{\Omega} \cdot \nabla f)}{|| \nabla f ||^2} =
    \frac{\bar{y}v + \bar{z}w}{R^2}

Substituting equations :eq:`reflection-cylinder-constant` and
:eq:`reflection-cylinder-grad` into equation :eq:`reflection-system` gives us
the form of the solution. In this case, the x-component will not change. The y-
and z-components of the reflected direction will be

.. math:: :label: reflection-cylinder

    v' = v - \frac{2 ( \bar{y}v + \bar{z}w ) \bar{y}}{R^2} \\

    w' = w - \frac{2 ( \bar{y}v + \bar{z}w ) \bar{z}}{R^2}


Sphere
------

The surface equation for a sphere has the form :math:`f(x,y,z) = (x - x_0)^2 +
(y - y_0)^2 + (z - z_0)^2 - R^2 = 0`. Thus, the gradient to the surface is

.. math:: :label: reflection-sphere-grad

    \nabla f = 2 \left ( \begin{array}{c} x - x_0 \\ y - y_0 \\ z - z_0
    \end{array} \right ) = 2 \left ( \begin{array}{c} \bar{x} \\ \bar{y} \\
    \bar{z} \end{array} \right )

where we have introduced the constants :math:`\bar{x}, \bar{y}, \bar{z}`. Taking
the square of the norm of the gradient, we find that

.. math:: :label: reflection-sphere-norm

    || \nabla f ||^2 = 4 \bar{x}^2 + 4 \bar{y}^2 + 4 \bar{z}^2 = 4 R^2

This imples that

.. math:: :label: reflection-sphere-constant

    \frac{2 (\mathbf{\Omega} \cdot \nabla f)}{|| \nabla f ||^2} =
    \frac{\bar{x}u + \bar{y}v + \bar{z}w}{R^2}

Substituting equations :eq:`reflection-sphere-constant` and
:eq:`reflection-sphere-grad` into equation :eq:`reflection-system` gives us the
form of the solution:

.. math:: :label: reflection-sphere

    u' = u - \frac{2 ( \bar{x}u + \bar{y}v + \bar{z}w ) \bar{x} }{R^2} \\

    v' = v - \frac{2 ( \bar{x}u + \bar{y}v + \bar{z}w ) \bar{y} }{R^2} \\

    w' = w - \frac{2 ( \bar{x}u + \bar{y}v + \bar{z}w ) \bar{z} }{R^2} \\


.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _surfaces: http://en.wikipedia.org/wiki/Surface
