from math import cos, sin, radians
import random

import openmc
import pytest


@pytest.mark.parametrize("angle", [30., 45., 60., 90., 120.])
def test_rotational_periodic_bc(angle):
    # Pick random starting angle
    start = random.uniform(0., 360.)
    degrees = angle
    ang1 = radians(start)
    ang2 = radians(start + degrees)

    # Define three points on each plane and then randomly shuffle them
    p1_points = [(0., 0., 0.), (cos(ang1), sin(ang1), 0.), (0., 0., 1.)]
    p2_points = [(0., 0., 0.), (cos(ang2), sin(ang2), 0.), (0., 0., 1.)]
    random.shuffle(p1_points)
    random.shuffle(p2_points)

    # Create periodic planes and a cylinder
    p1 = openmc.Plane.from_points(*p1_points, boundary_type='periodic')
    p2 = openmc.Plane.from_points(*p2_points, boundary_type='periodic')
    p1.periodic_surface = p2
    zcyl = openmc.ZCylinder(r=5., boundary_type='vacuum')

    # Figure out which side of planes to use based on a point in the middle
    ang_mid = radians(start + degrees/2.)
    mid_point = (cos(ang_mid), sin(ang_mid), 0.)
    r1 = -p1 if mid_point in -p1 else +p1
    r2 = -p2 if mid_point in -p2 else +p2

    # Create one cell bounded by the two planes and the cylinder
    mat = openmc.Material(density=1.0, density_units='g/cm3', components={'U235': 1.0})
    cell = openmc.Cell(fill=mat, region=r1 & r2 & -zcyl)

    # Make the model complete
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point(mid_point))
    model.settings.particles = 1000
    model.settings.batches = 10
    model.settings.inactive = 5

    # Run the model
    model.run()
