from functools import partial
from random import uniform, seed

import numpy as np
import openmc
import pytest


def assert_infinite_bb(s):
    ll, ur = (-s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))


def test_plane():
    s = openmc.Plane(a=1, b=2, c=-1, d=3, name='my plane')
    assert s.a == 1
    assert s.b == 2
    assert s.c == -1
    assert s.d == 3
    assert s.boundary_type == 'transmission'
    assert s.name == 'my plane'
    assert s.type == 'plane'

    # Generic planes don't have well-defined bounding boxes
    assert_infinite_bb(s)

    # evaluate method
    x, y, z = (4, 3, 6)
    assert s.evaluate((x, y, z)) == pytest.approx(s.a*x + s.b*y + s.c*z - s.d)

    # translate method
    st = s.translate((1.0, 0.0, 0.0))
    assert (st.a, st.b, st.c, st.d) == (s.a, s.b, s.c, 4)

    # Make sure repr works
    repr(s)


def test_plane_from_points():
    # Generate the plane x - y = 1 given three points
    p1 = (0, -1, 0)
    p2 = (1, 0, 0)
    p3 = (1, 0, 1)
    s = openmc.Plane.from_points(p1, p2, p3)

    # Confirm correct coefficients
    assert s.a == 1.0
    assert s.b == -1.0
    assert s.c == 0.0
    assert s.d == 1.0


def test_xplane():
    s = openmc.XPlane(x0=3., boundary_type='reflective')
    assert s.x0 == 3.
    assert s.boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert ll == pytest.approx((3., -np.inf, -np.inf))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ur == pytest.approx((3., np.inf, np.inf))
    assert np.all(np.isinf(ll))

    # __contains__ on associated half-spaces
    assert (5, 0, 0) in +s
    assert (5, 0, 0) not in -s
    assert (-2, 1, 10) in -s
    assert (-2, 1, 10) not in +s

    # evaluate method
    assert s.evaluate((5., 0., 0.)) == pytest.approx(2.)

    # translate method
    st = s.translate((1.0, 0.0, 0.0))
    assert st.x0 == s.x0 + 1

    # Make sure repr works
    repr(s)


def test_yplane():
    s = openmc.YPlane(y0=3.)
    assert s.y0 == 3.

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert ll == pytest.approx((-np.inf, 3., -np.inf))
    assert np.all(np.isinf(ur))
    ll, ur = s.bounding_box('-')
    assert ur == pytest.approx((np.inf, 3., np.inf))
    assert np.all(np.isinf(ll))

    # __contains__ on associated half-spaces
    assert (0, 5, 0) in +s
    assert (0, 5, 0) not in -s
    assert (-2, 1, 10) in -s
    assert (-2, 1, 10) not in +s

    # evaluate method
    assert s.evaluate((0., 0., 0.)) == pytest.approx(-3.)

    # translate method
    st = s.translate((0.0, 1.0, 0.0))
    assert st.y0 == s.y0 + 1


def test_zplane():
    s = openmc.ZPlane(z0=3.)
    assert s.z0 == 3.

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert ll == pytest.approx((-np.inf, -np.inf, 3.))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ur == pytest.approx((np.inf, np.inf, 3.))
    assert np.all(np.isinf(ll))

    # __contains__ on associated half-spaces
    assert (0, 0, 5) in +s
    assert (0, 0, 5) not in -s
    assert (-2, 1, -10) in -s
    assert (-2, 1, -10) not in +s

    # evaluate method
    assert s.evaluate((0., 0., 10.)) == pytest.approx(7.)

    # translate method
    st = s.translate((0.0, 0.0, 1.0))
    assert st.z0 == s.z0 + 1

    # Make sure repr works
    repr(s)


def test_xcylinder():
    y, z, r = 3, 5, 2
    s = openmc.XCylinder(y0=y, z0=z, r=r)
    assert s.y0 == y
    assert s.z0 == z
    assert s.r == r

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx((-np.inf, y-r, z-r))
    assert ur == pytest.approx((np.inf, y+r, z+r))

    # evaluate method
    assert s.evaluate((0, y, z)) == pytest.approx(-r**2)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.y0 == s.y0 + 1
    assert st.z0 == s.z0 + 1
    assert st.r == s.r

    # Make sure repr works
    repr(s)


def test_periodic():
    x = openmc.XPlane(boundary_type='periodic')
    y = openmc.YPlane(boundary_type='periodic')
    x.periodic_surface = y
    assert y.periodic_surface == x
    with pytest.raises(TypeError):
        x.periodic_surface = openmc.Sphere()


def test_ycylinder():
    x, z, r = 3, 5, 2
    s = openmc.YCylinder(x0=x, z0=z, r=r)
    assert s.x0 == x
    assert s.z0 == z
    assert s.r == r

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx((x-r, -np.inf, z-r))
    assert ur == pytest.approx((x+r, np.inf, z+r))

    # evaluate method
    assert s.evaluate((x, 0, z)) == pytest.approx(-r**2)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.z0 == s.z0 + 1
    assert st.r == s.r


def test_zcylinder():
    x, y, r = 3, 5, 2
    s = openmc.ZCylinder(x0=x, y0=y, r=r)
    assert s.x0 == x
    assert s.y0 == y
    assert s.r == r

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx((x-r, y-r, -np.inf))
    assert ur == pytest.approx((x+r, y+r, np.inf))

    # evaluate method
    assert s.evaluate((x, y, 0)) == pytest.approx(-r**2)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.y0 == s.y0 + 1
    assert st.r == s.r

    # Make sure repr works
    repr(s)


def test_sphere():
    x, y, z, r = -3, 5, 6, 2
    s = openmc.Sphere(x0=x, y0=y, z0=z, r=r)
    assert s.x0 == x
    assert s.y0 == y
    assert s.z0 == z
    assert s.r == r

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll == pytest.approx((x-r, y-r, z-r))
    assert ur == pytest.approx((x+r, y+r, z+r))

    # evaluate method
    assert s.evaluate((x, y, z)) == pytest.approx(-r**2)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.y0 == s.y0 + 1
    assert st.z0 == s.z0 + 1
    assert st.r == s.r

    # Make sure repr works
    repr(s)


def cone_common(apex, r2, cls):
    x, y, z = apex
    s = cls(x0=x, y0=y, z0=z, r2=r2)
    assert s.x0 == x
    assert s.y0 == y
    assert s.z0 == z
    assert s.r2 == r2

    # Check bounding box
    assert_infinite_bb(s)

    # evaluate method -- should be zero at apex
    assert s.evaluate((x, y, z)) == pytest.approx(0.0)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.y0 == s.y0 + 1
    assert st.z0 == s.z0 + 1
    assert st.r2 == s.r2

    # Make sure repr works
    repr(s)


def test_xcone():
    apex = (10, 0, 0)
    r2 = 4
    cone_common(apex, r2, openmc.XCone)


def test_ycone():
    apex = (10, 0, 0)
    r2 = 4
    cone_common(apex, r2, openmc.YCone)


def test_zcone():
    apex = (10, 0, 0)
    r2 = 4
    cone_common(apex, r2, openmc.ZCone)


def test_quadric():
    # Make a sphere from a quadric
    r = 10.0
    coeffs = {'a': 1, 'b': 1, 'c': 1, 'k': -r**2}
    s = openmc.Quadric(**coeffs)
    assert s.a == coeffs['a']
    assert s.b == coeffs['b']
    assert s.c == coeffs['c']
    assert s.k == coeffs['k']

    # All other coeffs should be zero
    for coeff in ('d', 'e', 'f', 'g', 'h', 'j'):
        assert getattr(s, coeff) == 0.0

    # Check bounding box
    assert_infinite_bb(s)

    # evaluate method
    assert s.evaluate((0., 0., 0.)) == pytest.approx(coeffs['k'])
    assert s.evaluate((1., 1., 1.)) == pytest.approx(3 + coeffs['k'])

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    for coeff in 'abcdef':
        assert getattr(s, coeff) == getattr(st, coeff)
    assert (st.g, st.h, st.j) == (-2, -2, -2)
    assert st.k == s.k + 3


def test_cylinder_from_points():
    seed(1)  # Make random numbers reproducible
    for _ in range(100):
        # Generate cylinder in random direction
        xi = partial(uniform, -10.0, 10.0)
        p1 = np.array([xi(), xi(), xi()])
        p2 = np.array([xi(), xi(), xi()])
        r = uniform(1.0, 100.0)
        s = openmc.model.cylinder_from_points(p1, p2, r)

        # Points p1 and p2 need to be inside cylinder
        assert p1 in -s
        assert p2 in -s

        # Points further along the line should be inside cylinder as well
        t = uniform(-100.0, 100.0)
        p = p1 + t*(p2 - p1)
        assert p in -s

        # Check that points outside cylinder are in positive half-space and
        # inside are in negative half-space. We do this by constructing a plane
        # that includes the cylinder's axis, finding the normal to the plane,
        # and using it to find a point slightly more/less than one radius away
        # from the axis.
        plane = openmc.Plane.from_points(p1, p2, (0., 0., 0.))
        n = np.array([plane.a, plane.b, plane.c])
        n /= np.linalg.norm(n)
        assert p1 + 1.1*r*n in +s
        assert p2 + 1.1*r*n in +s
        assert p1 + 0.9*r*n in -s
        assert p2 + 0.9*r*n in -s


def test_cylinder_from_points_axis():
    # Create axis-aligned cylinders and confirm the coefficients are as expected

    # (x - 3)^2 + (y - 4)^2 = 2^2
    # x^2 + y^2 - 6x - 8y + 21 = 0
    s = openmc.model.cylinder_from_points((3., 4., 0.), (3., 4., 1.), 2.)
    assert (s.a, s.b, s.c) == pytest.approx((1., 1., 0.))
    assert (s.d, s.e, s.f) == pytest.approx((0., 0., 0.))
    assert (s.g, s.h, s.j) == pytest.approx((-6., -8., 0.))
    assert s.k == pytest.approx(21.)

    # (y + 7)^2 + (z - 1)^2 = 3^2
    # y^2 + z^2 + 14y - 2z + 41 = 0
    s = openmc.model.cylinder_from_points((0., -7, 1.), (1., -7., 1.), 3.)
    assert (s.a, s.b, s.c) == pytest.approx((0., 1., 1.))
    assert (s.d, s.e, s.f) == pytest.approx((0., 0., 0.))
    assert (s.g, s.h, s.j) == pytest.approx((0., 14., -2.))
    assert s.k == 41.

    # (x - 2)^2 + (z - 5)^2 = 4^2
    # x^2 + z^2 - 4x - 10z + 13 = 0
    s = openmc.model.cylinder_from_points((2., 0., 5.), (2., 1., 5.), 4.)
    assert (s.a, s.b, s.c) == pytest.approx((1., 0., 1.))
    assert (s.d, s.e, s.f) == pytest.approx((0., 0., 0.))
    assert (s.g, s.h, s.j) == pytest.approx((-4., 0., -10.))
    assert s.k == pytest.approx(13.)
