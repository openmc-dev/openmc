from functools import partial
from random import uniform, seed

import numpy as np
import math
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

    # rotate method
    yp = openmc.YPlane(abs(s.d)/math.sqrt(s.a**2 + s.b**2 + s.c**2))
    psi = math.degrees(math.atan2(1, 2))
    phi = math.degrees(math.atan2(1, math.sqrt(5)))
    sr = s.rotate((phi, 0., psi), order='zyx')
    assert yp.normalize() == pytest.approx(sr.normalize())
    # test rotation ordering
    phi = math.degrees(math.atan2(1, math.sqrt(2)))
    sr = s.rotate((0., -45., phi), order='xyz')
    assert yp.normalize() == pytest.approx(sr.normalize())

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
    s = openmc.XPlane(3., boundary_type='reflective')
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

    # rotate method
    # make sure rotating around x axis does nothing to coefficients
    sr = s.rotate((37.4, 0., 0.))
    assert s._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())
    # rotating around z by 90 deg then x by -90 deg should give negative z-plane
    sr = s.rotate((-90., 0., 90), order='zyx')
    assert (0., 0., -1., 3.) == pytest.approx(sr._get_base_coeffs())

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

    # rotate method
    # make sure rotating around y axis does nothing to coefficients
    sr = s.rotate((0., -12.4, 0.), order='yxz')
    assert s._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())
    # rotate around x by -90 deg and y by 90 deg should give negative x-plane
    sr = s.rotate((-90, 90, 0.))
    assert (-1, 0., 0., 3.) == pytest.approx(sr._get_base_coeffs())

    # Make sure repr works
    repr(s)


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

    # rotate method
    # make sure rotating around z axis does nothing to coefficients
    sr = s.rotate((0., 0., 123), order='zxy')
    assert s._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())
    # rotate around x by -90 deg and y by 90 deg should give negative x-plane
    sr = s.rotate((-90, 0., 90.))
    assert (-1., 0., 0., 3.) == pytest.approx(sr._get_base_coeffs())

    # Make sure repr works
    repr(s)


def test_cylinder():
    x0, y0, z0, r = 2, 3, 4, 2
    dx, dy, dz = 1, -1, 1
    s = openmc.Cylinder(x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, r=r)
    assert s.x0 == 2
    assert s.y0 == 3
    assert s.z0 == 4
    assert s.dx == 1
    assert s.dy == -1
    assert s.dz == 1
    assert s.r == 2

    # Check bounding box
    assert_infinite_bb(s)

    # evaluate method
    # |(p - p1) ⨯ (p - p2)|^2 / |p2 - p1|^2 - r^2
    p1 = s._origin
    p2 = p1 + s._axis
    perp = np.array((1, -2, 1))*(1 / s._axis)
    divisor = np.linalg.norm(p2 - p1)
    pin = p1 + 5*s._axis # point inside cylinder
    pout = np.array((4., 0., 2.5)) # point outside the cylinder
    pon = p1 + s.r*perp / np.linalg.norm(perp) # point on cylinder
    for p, fn in zip((pin, pout, pon), (np.less, np.greater, np.isclose)):
        c1 = np.linalg.norm(np.cross(p - p1, p - p2)) / divisor
        val = c1*c1 - s.r*s.r
        p_eval = s.evaluate(p)
        assert fn(p_eval, 0.)
        assert p_eval == pytest.approx(val)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.y0 == s.y0 + 1
    assert st.z0 == s.z0 + 1
    assert st.dx == s.dx
    assert st.dy == s.dy
    assert st.dz == s.dz
    assert st.r == s.r

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2.is_equal(sr)

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

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

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

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

    # Make sure repr works
    repr(s)


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

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

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

    # rotate method
    pivot = np.array([1, -2, 3])
    sr = s.rotate((90, 90, 90), pivot=pivot)
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx((R @ (s._origin - pivot)) + pivot)
    # test passing in rotation matrix
    sr2 = s.rotate(R, pivot=pivot)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

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

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

    # Make sure repr works
    repr(s)


def test_cone():
    x0, y0, z0, r2 = 2, 3, 4, 4
    dx, dy, dz = 1, -1, 1
    s = openmc.Cone(x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, r2=r2)
    assert s.x0 == 2
    assert s.y0 == 3
    assert s.z0 == 4
    assert s.dx == 1
    assert s.dy == -1
    assert s.dz == 1
    assert s.r2 == 4

    # Check bounding box
    assert_infinite_bb(s)

    # evaluate method
    # cos^2(theta) * ((p - p1))**2 - (d @ (p - p1))^2
    # The argument r2 for cones is actually tan^2(theta) so that
    # cos^2(theta) = 1 / (1 + r2)
    #
    # This makes the evaluation equation shown below where p is the evaluation
    # point (x, y, z) p1 is the apex (origin) of the cone and r2 is related to
    # the aperature of the cone as described above
    # (p - p1) @ (p - p1) / (1 + r2) - (d @ (p - p1))^2
    # point inside
    p1 = s._origin
    d = s._axis
    perp = np.array((1, -2, 1))*(1 / d)
    perp /= np.linalg.norm(perp)
    pin = p1 + 5*d # point inside cone
    pout = p1 + 3.2*perp # point outside cone
    pon = p1 + 3.2*d + 3.2*math.sqrt(s.r2)*perp # point on cone
    for p, fn in zip((pin, pout, pon), (np.less, np.greater, np.isclose)):
        val = np.sum((p - p1)**2) / (1 + s.r2) - np.sum((d @ (p - p1))**2)
        p_eval = s.evaluate(p)
        assert fn(p_eval, 0.)
        assert p_eval == pytest.approx(val)

    # translate method
    st = s.translate((1.0, 1.0, 1.0))
    assert st.x0 == s.x0 + 1
    assert st.y0 == s.y0 + 1
    assert st.z0 == s.z0 + 1
    assert st.dx == s.dx
    assert st.dy == s.dy
    assert st.dz == s.dz
    assert st.r2 == s.r2

    # rotate method
    sr = s.rotate((90, 90, 90))
    R = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    assert sr._origin == pytest.approx(R @ s._origin)
    assert sr._axis == pytest.approx(R @ s._axis)
    # test passing in rotation matrix
    sr2 = s.rotate(R)
    assert sr2._get_base_coeffs() == pytest.approx(sr._get_base_coeffs())

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
    assert openmc.Sphere(r=10).is_equal(s)

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

    # rotate method
    x0, y0, z0, r2 = 2, 3, 4, 4
    dx, dy, dz = 1, -1, 1
    s = openmc.Cone(x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, r2=r2)
    q = openmc.Quadric(*s._get_base_coeffs())
    qr = q.rotate((45, 60, 30))
    sr = s.rotate((45, 60, 30))
    assert qr.is_equal(sr)


def test_cylinder_from_points():
    seed(1)  # Make random numbers reproducible
    for _ in range(100):
        # Generate cylinder in random direction
        xi = partial(uniform, -10.0, 10.0)
        p1 = np.array([xi(), xi(), xi()])
        p2 = np.array([xi(), xi(), xi()])
        r = uniform(1.0, 100.0)
        s = openmc.Cylinder.from_points(p1, p2, r)

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
    s = openmc.Cylinder.from_points((3., 4., 0.), (3., 4., 1.), 2.)
    a, b, c, d, e, f, g, h, j, k = s._get_base_coeffs()
    assert (a, b, c) == pytest.approx((1., 1., 0.))
    assert (d, e, f) == pytest.approx((0., 0., 0.))
    assert (g, h, j) == pytest.approx((-6., -8., 0.))
    assert k == pytest.approx(21.)

    # (y + 7)^2 + (z - 1)^2 = 3^2
    # y^2 + z^2 + 14y - 2z + 41 = 0
    s = openmc.Cylinder.from_points((0., -7, 1.), (1., -7., 1.), 3.)
    a, b, c, d, e, f, g, h, j, k = s._get_base_coeffs()
    assert (a, b, c) == pytest.approx((0., 1., 1.))
    assert (d, e, f) == pytest.approx((0., 0., 0.))
    assert (g, h, j) == pytest.approx((0., 14., -2.))
    assert k == 41.

    # (x - 2)^2 + (z - 5)^2 = 4^2
    # x^2 + z^2 - 4x - 10z + 13 = 0
    s = openmc.Cylinder.from_points((2., 0., 5.), (2., 1., 5.), 4.)
    a, b, c, d, e, f, g, h, j, k = s._get_base_coeffs()
    assert (a, b, c) == pytest.approx((1., 0., 1.))
    assert (d, e, f) == pytest.approx((0., 0., 0.))
    assert (g, h, j) == pytest.approx((-4., 0., -10.))
    assert k == pytest.approx(13.)
