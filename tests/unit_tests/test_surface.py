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
    s = openmc.Plane(A=1, B=2, C=-1, D=3, name='my plane')
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

    # Make sure repr works
    repr(s)


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

    # Make sure repr works
    repr(s)


def test_xcylinder():
    y, z, r = 3, 5, 2
    s = openmc.XCylinder(y0=y, z0=z, R=r)
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
    s = openmc.YCylinder(x0=x, z0=z, R=r)
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


def test_zcylinder():
    x, y, r = 3, 5, 2
    s = openmc.ZCylinder(x0=x, y0=y, R=r)
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

    # Make sure repr works
    repr(s)


def test_sphere():
    x, y, z, r = -3, 5, 6, 2
    s = openmc.Sphere(x0=x, y0=y, z0=z, R=r)
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

    # Make sure repr works
    repr(s)


def cone_common(apex, r2, cls):
    x, y, z = apex
    s = cls(x0=x, y0=y, z0=z, R2=r2)
    assert s.x0 == x
    assert s.y0 == y
    assert s.z0 == z
    assert s.r2 == r2

    # Check bounding box
    assert_infinite_bb(s)

    # evaluate method -- should be zero at apex
    assert s.evaluate((x, y, z)) == pytest.approx(0.0)

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
