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
    s_r = openmc.model.RightCircularCylinder((x, y, z), h, r, axis=axis.lower(),
                                             upper_fillet_radius=1.6,
                                             lower_fillet_radius=1.6)
    for s in (s, s_r):
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
    "axis, indices, center", [
        ("X", [2, 0, 1], (0., 0.)),
        ("Y", [0, 2, 1], (0., 0.)),
        ("Z", [0, 1, 2], (0., 0.)),
        ("X", [2, 0, 1], (10., 5.)),
        ("Y", [0, 2, 1], (10., 5.)),
        ("Z", [0, 1, 2], (10., 5.)),

    ]
)
def test_cylinder_sector(axis, indices, center):
    c1, c2 = center
    r1, r2 = 0.5, 1.5
    d = (r2 - r1) / 2
    phi1 = -60.
    phi2 = 60
    s = openmc.model.CylinderSector(r1, r2, phi1, phi2, center=center,
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
    test_point_ll = np.array([-r2 + c1, -r2 + c2, -np.inf])
    assert ll == pytest.approx(test_point_ll[indices])
    test_point_ur = np.array([r2 + c1, r2 + c2, np.inf])
    assert ur == pytest.approx(test_point_ur[indices])

    # __contains__ on associated half-spaces
    point_pos = np.array([0 + c1, r2 + 1 + c2, 0])
    assert point_pos[indices] in +s
    assert point_pos[indices] not in -s
    point_neg = np.array([r1 + d + c1, r1 + d + c2, 0])
    assert point_neg[indices] in -s
    assert point_neg[indices] not in +s

    # Check __contains__ for sector with reflex angle
    s_reflex = openmc.model.CylinderSector(
        r1, r2, 0., 270., center=center, axis=axis.lower())
    points = [
        np.array([c1 + r1 + d, c2 + 0.01, 0.]),
        np.array([c1, c2 + r1 + d, 0.]),
        np.array([c1 - r1 - d, c2, 0.]),
        np.array([c1 - 0.01, c2 - r1 - d, 0.])
    ]
    for point_neg in points:
        assert point_neg[indices] in -s_reflex
        assert point_neg[indices] not in +s_reflex

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
        ("y", "Z", "X", 1),
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


def test_polygon():
    # define a 5 pointed star centered on 1, 1
    star = np.array([[1.        , 2.        ],
                     [0.70610737, 1.4045085 ],
                     [0.04894348, 1.30901699],
                     [0.52447174, 0.8454915 ],
                     [0.41221475, 0.19098301],
                     [1.        , 0.5       ],
                     [1.58778525, 0.19098301],
                     [1.47552826, 0.8454915 ],
                     [1.95105652, 1.30901699],
                     [1.29389263, 1.4045085 ],
                     [1.        , 2.        ]])
    points_in = [(1, 1, 0), (0, 1, 1), (1, 0, 1), (.707, .707, 1)]
    for i, basis in enumerate(('xy', 'yz', 'xz', 'rz')):
        star_poly = openmc.model.Polygon(star, basis=basis)
        assert points_in[i] in -star_poly
        assert any([points_in[i] in reg for reg in star_poly.regions])
        assert points_in[i] not in +star_poly
        assert (0, 0, 0) not in -star_poly
        if basis != "rz":
            for offsets in [0.6, np.array([0.6] * 10), [0.6] * 10]:
                offset_star = star_poly.offset(offsets)
                assert (0, 0, 0) in -offset_star
                assert any([(0, 0, 0) in reg for reg in offset_star.regions])
            with pytest.raises(ValueError):
                star_poly.offset([0.6, 0.6])

    # check invalid Polygon input points
    # duplicate points not just at start and end
    rz_points = np.array([[6.88, 3.02],
                          [6.88, 2.72],
                          [6.88, 3.02],
                          [7.63, 0.0],
                          [5.75, 0.0],
                          [5.75, 1.22],
                          [7.63, 0.0],
                          [6.30, 1.22],
                          [6.30, 3.02],
                          [6.88, 3.02]])
    with pytest.raises(ValueError):
        openmc.model.Polygon(rz_points)

    # segment traces back on previous segment
    rz_points = np.array([[6.88, 3.02],
                          [6.88, 2.72],
                          [6.88, 2.32],
                          [6.88, 2.52],
                          [7.63, 0.0],
                          [5.75, 0.0],
                          [6.75, 0.0],
                          [5.75, 1.22],
                          [6.30, 1.22],
                          [6.30, 3.02],
                          [6.88, 3.02]])
    with pytest.raises(ValueError):
        openmc.model.Polygon(rz_points)

    # segments intersect (line-line)
    rz_points = np.array([[6.88, 3.02],
                          [5.88, 2.32],
                          [7.63, 0.0],
                          [5.75, 0.0],
                          [5.75, 1.22],
                          [6.30, 1.22],
                          [6.30, 3.02],
                          [6.88, 3.02]])
    with pytest.raises(ValueError):
        openmc.model.Polygon(rz_points)

    # segments intersect (line-point)
    rz_points = np.array([[6.88, 3.02],
                          [6.3, 2.32],
                          [7.63, 0.0],
                          [5.75, 0.0],
                          [5.75, 1.22],
                          [6.30, 1.22],
                          [6.30, 3.02],
                          [6.88, 3.02]])
    with pytest.raises(ValueError):
        openmc.model.Polygon(rz_points)

    # Test "M" shaped polygon
    points = np.array([[8.5151581, -17.988337],
                       [10.381711000000001, -17.988337],
                       [12.744357, -24.288728000000003],
                       [15.119406000000001, -17.988337],
                       [16.985959, -17.988337],
                       [16.985959, -27.246687],
                       [15.764328, -27.246687],
                       [15.764328, -19.116951],
                       [13.376877, -25.466951],
                       [12.118039, -25.466951],
                       [9.7305877, -19.116951],
                       [9.7305877, -27.246687],
                       [8.5151581, -27.246687]])

    # Test points inside and outside by using offset method
    m_polygon = openmc.model.Polygon(points, basis='xz')
    inner_pts = m_polygon.offset(-0.1).points
    assert all([(pt[0], 0, pt[1]) in -m_polygon for pt in inner_pts])
    outer_pts = m_polygon.offset(0.1).points
    assert all([(pt[0], 0, pt[1]) in +m_polygon for pt in outer_pts])

    # Offset of -0.2 will cause self-intersection
    with pytest.raises(ValueError):
        m_polygon.offset(-0.2)


@pytest.mark.parametrize("axis", ["x", "y", "z"])
def test_cruciform_prism(axis):
    center = x0, y0 = (3., 4.)
    distances = [2., 3., 5.]
    s = openmc.model.CruciformPrism(distances, center, axis=axis)

    if axis == 'x':
        i1, i2 = 1, 2
    elif axis == 'y':
        i1, i2 = 0, 2
    elif axis == 'z':
        i1, i2 = 0, 1
    plane_cls = (openmc.XPlane, openmc.YPlane, openmc.ZPlane)

    # Check type of surfaces
    for i in range(3):
        assert isinstance(getattr(s, f'hmin{i}'), plane_cls[i1])
        assert isinstance(getattr(s, f'hmax{i}'), plane_cls[i1])
        assert isinstance(getattr(s, f'vmin{i}'), plane_cls[i2])
        assert isinstance(getattr(s, f'vmax{i}'), plane_cls[i2])

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    for i in range(3):
        assert getattr(s, f'hmin{i}').boundary_type == 'reflective'
        assert getattr(s, f'hmax{i}').boundary_type == 'reflective'
        assert getattr(s, f'vmin{i}').boundary_type == 'reflective'
        assert getattr(s, f'vmax{i}').boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ur[i1] == pytest.approx(x0 + distances[-1])
    assert ur[i2] == pytest.approx(y0 + distances[-1])
    assert ll[i1] == pytest.approx(x0 - distances[-1])
    assert ll[i2] == pytest.approx(y0 - distances[-1])

    # __contains__ on associated half-spaces
    point_pos, point_neg = np.zeros(3), np.zeros(3)
    point_pos[i1] = x0 + 3.1
    point_pos[i2] = y0 + 2.05
    point_neg[i1] = x0 + 3.5
    point_neg[i2] = y0 + 1.99
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

    # Check that non-monotonic distances fail
    with pytest.raises(ValueError):
        openmc.model.CruciformPrism([1.0, 0.5, 2.0, 3.0])
    with pytest.raises(ValueError):
        openmc.model.CruciformPrism([3.0, 2.0, 0.5, 1.0])


def test_box():
    v = (-1.0, -1.0, -2.5)
    a1 = (2.0, -1.0, 0.0)
    a2 = (1.0, 2.0, 0.0)
    a3 = (0.0, 0.0, 5.0)
    s = openmc.model.OrthogonalBox(v, a1, a2, a3)
    for num in (1, 2, 3):
        assert isinstance(getattr(s, f'ax{num}_min'), openmc.Plane)
        assert isinstance(getattr(s, f'ax{num}_max'), openmc.Plane)

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    for num in (1, 2, 3):
        assert getattr(s, f'ax{num}_min').boundary_type == 'reflective'
        assert getattr(s, f'ax{num}_max').boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll[2] == pytest.approx(-2.5)
    assert ur[2] == pytest.approx(2.5)

    # __contains__ on associated half-spaces
    assert (0., 0., 0.) in -s
    assert (-2., 0., 0.) not in -s
    assert (0., 0.9, 0.) in -s
    assert (0., 0., -3.) in +s
    assert (0., 0., 3.) in +s

    # translate method
    s_t = s.translate((1., 1., 0.))
    assert (-0.01, 0., 0.) in +s_t
    assert (0.01, 0., 0.) in -s_t

    # Make sure repr works
    repr(s)

    # Version with infinite 3rd dimension
    s = openmc.model.OrthogonalBox(v, a1, a2)
    assert not hasattr(s, 'ax3_min')
    assert not hasattr(s, 'ax3_max')
    ll, ur = (-s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    assert (0., 0., 0.) in -s
    assert (-2., 0., 0.) not in -s
    assert (0., 0.9, 0.) in -s
    assert (0., 0., -3.) not in +s
    assert (0., 0., 3.) not in +s


def test_conical_frustum():
    center_base = (0.0, 0.0, -3)
    axis = (0., 0., 3.)
    r1 = 2.0
    r2 = 0.5
    s = openmc.model.ConicalFrustum(center_base, axis, r1, r2)
    assert isinstance(s.cone, openmc.Cone)
    assert isinstance(s.plane_bottom, openmc.Plane)
    assert isinstance(s.plane_top, openmc.Plane)

    # Make sure boundary condition propagates
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.cone.boundary_type == 'reflective'
    assert s.plane_bottom.boundary_type == 'reflective'
    assert s.plane_top.boundary_type == 'reflective'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert ll[2] == pytest.approx(-3.0)
    assert ur[2] == pytest.approx(0.0)

    # __contains__ on associated half-spaces
    assert (0., 0., -1.) in -s
    assert (0., 0., -4.) not in -s
    assert (0., 0., 1.) not in -s
    assert (1., 1., -2.99) in -s
    assert (1., 1., -0.01) in +s

    # translate method
    s_t = s.translate((1., 1., 0.))
    assert (1., 1., -0.01) in -s_t

    # Make sure repr works
    repr(s)

    # Denegenerate case with r1 = r2
    s = openmc.model.ConicalFrustum(center_base, axis, r1, r1)
    assert (1., 1., -0.01) in -s


def test_vessel():
    center = (3.0, 2.0)
    r = 1.0
    p1, p2 = -5.0, 5.0
    h1 = h2 = 1.0
    s = openmc.model.Vessel(r, p1, p2, h1, h2, center)
    assert isinstance(s.cyl, openmc.Cylinder)
    assert isinstance(s.plane_bottom, openmc.Plane)
    assert isinstance(s.plane_top, openmc.Plane)
    assert isinstance(s.bottom, openmc.Quadric)
    assert isinstance(s.top, openmc.Quadric)

    # Make sure boundary condition propagates (but not for planes)
    s.boundary_type = 'reflective'
    assert s.boundary_type == 'reflective'
    assert s.cyl.boundary_type == 'reflective'
    assert s.bottom.boundary_type == 'reflective'
    assert s.top.boundary_type == 'reflective'
    assert s.plane_bottom.boundary_type == 'transmission'
    assert s.plane_top.boundary_type == 'transmission'

    # Check bounding box
    ll, ur = (+s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))
    ll, ur = (-s).bounding_box
    assert np.all(np.isinf(ll))
    assert np.all(np.isinf(ur))

    # __contains__ on associated half-spaces
    assert (3., 2., 0.) in -s
    assert (3., 2., -5.0) in -s
    assert (3., 2., 5.0) in -s
    assert (3., 2., -5.9) in -s
    assert (3., 2., 5.9) in -s
    assert (3., 2., -6.1) not in -s
    assert (3., 2., 6.1) not in -s
    assert (4.5, 2., 0.) in +s
    assert (3., 3.2, 0.) in +s
    assert (3., 2., 7.) in +s

    # translate method
    s_t = s.translate((0., 0., 1.))
    assert (3., 2., 6.1) in -s_t

    # Make sure repr works
    repr(s)
