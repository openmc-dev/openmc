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
    s = openmc.model.RectangularParallelepiped(
        xmin, xmax, ymin, ymax, zmin, zmax)
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


@pytest.mark.parametrize(
    "axis, indices", [
        ("X", [0, 1, 2]),
        ("Y", [1, 2, 0]),
        ("Z", [2, 0, 1]),
    ]
)
def test_cylinder_sector(axis, indices):
    r1, r2 = 0.5, 1.5
    d = (r2 - r1) / 2
    phi1 = -60.
    phi2 = 60
    s = openmc.model.CylinderSector(r1, r2, phi1, phi2,
                                    axis=axis.lower())
    assert isinstance(s.outer_cyl, getattr(openmc, axis + "Cylinder"))
    assert isinstance(s.inner_cyl, getattr(openmc, axis + "Cylinder"))
    assert isinstance(s.plane1, openmc.Plane)
    assert isinstance(s.plane2, openmc.Plane)

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.outer_cyl.boundary_type == 'reflective'
    assert s.inner_cyl.boundary_type == 'reflective'
    assert s.plane1.boundary_type == 'reflective'
    assert s.plane2.boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx(np.roll([-np.inf, -r2, -r2], indices[0]))
    assert ur == pytest.approx(np.roll([np.inf, r2, r2], indices[0]))

    # __contains__ on associated half-spaces
    point_pos = np.roll([0, r2 + 1, 0], indices[0])
    assert point_pos in +s
    assert point_pos not in -s
    point_neg = np.roll([0, r1 + d, r1 + d], indices[0])
    assert point_neg in -s
    assert point_neg not in +s

    # translate method
    t = uniform(-5.0, 5.0)
    s_t = s.translate((t, t, t))
    ll_t, ur_t = (-s_t).bounding_box
    assert ur_t == pytest.approx(ur + t)
    assert ll_t == pytest.approx(ll + t)

    # Check invalid r1, r2 combinations
    with pytest.raises(ValueError):
        openmc.model.CylinderSector(r2, r1, phi1, phi2)

    # Check invalid angles
    with pytest.raises(ValueError):
        openmc.model.CylinderSector(r1, r2, phi2, phi1)

    # Make sure repr works
    repr(s)


def test_cylinder_sector_from_theta_alpha():
    r1, r2 = 0.5, 1.5
    d = (r2 - r1) / 2
    theta = 120.
    alpha = -60.
    theta1 = alpha
    theta2 = alpha + theta
    s = openmc.model.CylinderSector(r1, r2, theta1, theta2)
    s_alt = openmc.model.CylinderSector.from_theta_alpha(r1,
                                                     r2,
                                                     theta,
                                                     alpha)

    # Check that the angles are correct
    assert s.plane1.coefficients == s_alt.plane1.coefficients
    assert s.plane2.coefficients == s_alt.plane2.coefficients
    assert s.inner_cyl.coefficients == s_alt.inner_cyl.coefficients
    assert s.outer_cyl.coefficients == s_alt.outer_cyl.coefficients

    # Check invalid sector width
    with pytest.raises(ValueError):
        openmc.model.CylinderSector.from_theta_alpha(r1, r2, 360, alpha)
    with pytest.raises(ValueError):
        openmc.model.CylinderSector.from_theta_alpha(r1, r2, -1, alpha)


@pytest.mark.parametrize(
    "axis, plane_tb, plane_lr, axis_idx", [
        ("x", "Z", "Y", 0),
        ("y", "X", "Z", 1),
        ("z", "Y", "X", 2),
    ]
)
def test_isogonal_octagon(axis, plane_tb, plane_lr, axis_idx):
    center = np.array([0., 0.])
    point_pos = np.array([0.8, 0.8])
    point_neg = np.array([0.7, 0.7])
    r1 = 1.
    r2 = 1.
    plane_top_bottom = getattr(openmc, plane_tb + "Plane")
    plane_left_right = getattr(openmc, plane_lr + "Plane")
    s = openmc.model.IsogonalOctagon(center, r1, r2, axis=axis)
    assert isinstance(s.top, plane_top_bottom)
    assert isinstance(s.bottom, plane_top_bottom)
    assert isinstance(s.right, plane_left_right)
    assert isinstance(s.left, plane_left_right)
    assert isinstance(s.upper_right, openmc.Plane)
    assert isinstance(s.lower_right, openmc.Plane)
    assert isinstance(s.upper_left, openmc.Plane)
    assert isinstance(s.lower_left, openmc.Plane)

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.top.boundary_type == 'reflective'
    assert s.bottom.boundary_type == 'reflective'
    assert s.right.boundary_type == 'reflective'
    assert s.left.boundary_type == 'reflective'
    assert s.upper_right.boundary_type == 'reflective'
    assert s.lower_right.boundary_type == 'reflective'
    assert s.lower_left.boundary_type == 'reflective'
    assert s.upper_left.boundary_type == 'reflective'

    # Check bounding box
    center = np.insert(center, axis_idx, np.inf)
    xmax, ymax, zmax = center + r1
    coord_min = center - r1
    coord_min[axis_idx] *= -1
    xmin, ymin, zmin = coord_min
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ur == pytest.approx((xmax, ymax, zmax))
    assert ll == pytest.approx((xmin, ymin, zmin))

    # __contains__ on associated half-spaces
    point_pos = np.insert(point_pos, axis_idx, 0)
    point_neg = np.insert(point_neg, axis_idx, 0)
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

    # Check invalid r1, r2 combinations
    with pytest.raises(ValueError):
        openmc.model.IsogonalOctagon(center, r1=1.0, r2=10.)
    with pytest.raises(ValueError):
        openmc.model.IsogonalOctagon(center, r1=10., r2=1.)

    # Make sure repr works
    repr(s)
