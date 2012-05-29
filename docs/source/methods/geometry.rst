.. _methods_geometry:

========
Geometry
========

OpenMC uses a technique known as `constructive solid geometry`_ (CSG) to build
arbitrarily complex 3D models. In a CSG model, every simply-connected region is
described as the union, intersection, or difference of half-spaces created by
bounding surfaces. In OpenMC, any second-order surface of the form

.. math::

    f(x,y,z) = Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0

can be modeled in OpenMC. For example, the equation for a sphere centered at
:math:`(\bar{x},\bar{y},\bar{z})` and of radius :math:`R` can be written as

.. math::

    f(x,y,z) = (x - \bar{x})^2 + (y - \bar{y})^2 + (z - \bar{z})^2 + R^2 = 0

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
