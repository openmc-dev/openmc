from random import uniform

import numpy as np
import openmc
import pytest


def test_rectangular_parallelepiped():
    xmin = uniform(-5., 5.)
    xmax = xmin + uniform(0., 5.)
    ymin = uniform(-5., 5.)
    ymax = ymin + uniform(0., 5.)
    zmin = uniform(-5., 5.)
    zmax = zmin + uniform(0., 5.)
    s = openmc.model.RectangularParallelepiped(xmin, xmax, ymin, ymax, zmin, zmax)
    assert isinstance(s.xmin, openmc.XPlane)
    assert isinstance(s.xmax, openmc.XPlane)
    assert isinstance(s.ymin, openmc.YPlane)
    assert isinstance(s.ymax, openmc.YPlane)
    assert isinstance(s.zmin, openmc.ZPlane)
    assert isinstance(s.zmax, openmc.ZPlane)

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    for axis in 'xyz':
        assert getattr(s, '{}min'.format(axis)).boundary_type == 'reflective'
        assert getattr(s, '{}max'.format(axis)).boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ur == pytest.approx((xmax, ymax, zmax))
    assert ll == pytest.approx((xmin, ymin, zmin))

    # __contains__ on associated half-spaces
    assert (xmin - 0.1, 0., 0.) in +s
    assert (xmin - 0.1, 0., 0.) not in -s
    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    assert (xmin + dx/2, ymin + dy/2, zmin + dz/2) in -s
    assert (xmin + dx/2, ymin + dy/2, zmin + dz/2) not in +s

    # translate method
    t = uniform(-5.0, 5.0)
    s_t = s.translate((t, t, t))
    ll_t, ur_t = (-s_t).bounding_box
    assert ur_t == pytest.approx(ur + t)
    assert ll_t == pytest.approx(ll + t)

    # Make sure repr works
    repr(s)


@pytest.mark.parametrize(
    "axis, indices", [
        ("X", [0, 1, 2]),
        ("Y", [1, 2, 0]),
        ("Z", [2, 0, 1]),
    ]
)
def test_right_circular_cylinder(axis, indices):
    x, y, z = 1.0, -2.5, 3.0
    h, r = 5.0, 3.0
    s = openmc.model.RightCircularCylinder((x, y, z), h, r, axis=axis.lower())
    assert isinstance(s.cyl, getattr(openmc, axis + "Cylinder"))
    assert isinstance(s.top, getattr(openmc, axis + "Plane"))
    assert isinstance(s.bottom, getattr(openmc, axis + "Plane"))

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.cyl.boundary_type == 'reflective'
    assert s.bottom.boundary_type == 'reflective'
    assert s.top.boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx((x, y, z) + np.roll([0, -r, -r], indices[0]))
    assert ur == pytest.approx((x, y, z) + np.roll([h, r, r], indices[0]))

    # __contains__ on associated half-spaces
    point_pos = (x, y, z) + np.roll([h/2, r+1, r+1], indices[0])
    assert point_pos in +s
    assert point_pos not in -s
    point_neg = (x, y, z) + np.roll([h/2, 0, 0], indices[0])
    assert point_neg in -s
    assert point_neg not in +s

    # translate method
    t = uniform(-5.0, 5.0)
    s_t = s.translate((t, t, t))
    ll_t, ur_t = (-s_t).bounding_box
    assert ur_t == pytest.approx(ur + t)
    assert ll_t == pytest.approx(ll + t)

    # Make sure repr works
    repr(s)


@pytest.mark.parametrize(
    "axis, point_pos, point_neg, ll_true", [
        ("X", (8., 0., 0.), (12., 0., 0.), (10., -np.inf, -np.inf)),
        ("Y", (10., -2., 0.), (10., 2., 0.), (-np.inf, 0., -np.inf)),
        ("Z", (10., 0., -3.), (10., 0., 3.), (-np.inf, -np.inf, 0.))
    ]
)
def test_cone_one_sided(axis, point_pos, point_neg, ll_true):
    cone_oneside = getattr(openmc.model, axis + "ConeOneSided")
    cone_twoside = getattr(openmc, axis + "Cone")
    plane = getattr(openmc, axis + "Plane")

    x, y, z = 10., 0., 0.
    r2 = 4.
    s = cone_oneside(x, y, z, r2, True)
    assert isinstance(s.cone, cone_twoside)
    assert isinstance(s.plane, plane)
    assert s.up

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.cone.boundary_type == 'reflective'
    assert s.plane.boundary_type == 'transmission'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert np.all(np.isinf(ur))
    assert ll == pytest.approx(ll_true)

    # __contains__ on associated half-spaces
    assert point_pos in +s
    assert point_pos not in -s
    assert point_neg in -s
    assert point_neg not in +s

    # translate method
    t = uniform(-5.0, 5.0)
    s_t = s.translate((t, t, t))
    ll_t, ur_t = (-s_t).bounding_box
    assert ur_t == pytest.approx(ur + t)
    assert ll_t == pytest.approx(ll + t)

    # Make sure repr works
    repr(s)
