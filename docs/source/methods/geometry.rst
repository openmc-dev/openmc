.. _methods_geometry:

========
Geometry
========

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


-------------------
Reflective Surfaces
-------------------

In general, a surface can be written in the form :math:`f(x,y,z) = 0`. If a
neutron is traveling in direction :math:`\mathbf{v}` and crosses a reflective
surface of the above form, it can be shown that the velocity vector will then
become

.. math::

    \mathbf{v'} = \mathbf{v} - 2 (\mathbf{v} \cdot \hat{\mathbf{n}})
    \hat{\mathbf{n}}

where :math:`\hat{\mathbf{n}}` is a unit vector normal to the surface at the
point of the surface crossing. The direction of the surface normal will be the
gradient to the surface at the point of crossing, i.e. :math:`\mathbf{n} =
\nabla f(x,y,z)`.

.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry
.. _surfaces: http://en.wikipedia.org/wiki/Surface
