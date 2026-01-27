"""
This test is designed to ensure that we account for potential corner crossings
in floating point precision.

"""

from math import pi, cos, sin

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()

    # Length of lattice on each side
    lat_size = 20.0

    # Angle that we're crossing the corner at
    phi = pi/4.0

    air = openmc.Material()
    air.set_density('g/cm3', 0.001)
    air.add_nuclide('N14', 1.0)
    metal = openmc.Material()
    metal.set_density('g/cm3', 7.0)
    metal.add_nuclide('Fe56', 1.0)

    metal_cell = openmc.Cell(fill=metal)
    metal_uni = openmc.Universe(cells=[metal_cell])

    air_cell = openmc.Cell(fill=air)
    air_uni = openmc.Universe(cells=[air_cell])

    # Define a checkerboard lattice
    lattice = openmc.RectLattice()
    lattice.lower_left = (-lat_size/2.0, -lat_size/2.0)
    lattice.pitch = (lat_size/2, lat_size/2)
    lattice.universes = [
        [metal_uni, air_uni],
        [air_uni, metal_uni]
    ]

    box = openmc.model.RectangularPrism(lat_size, lat_size)
    cyl = openmc.ZCylinder(r=lat_size, boundary_type='vacuum')
    outside_lattice = openmc.Cell(region=-cyl & +box, fill=air)
    inside_lattice = openmc.Cell(region=-box, fill=lattice)

    model.geometry = openmc.Geometry([outside_lattice, inside_lattice])

    # Set all runtime parameters
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 1000

    # Define a source located outside the lattice and pointing straight into its
    # corner at 45 degrees
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((-cos(phi), -sin(phi), 0.0)),
        angle=openmc.stats.Monodirectional((cos(phi), sin(phi), 0.0))
    )

    # Create a mesh tally
    mesh = openmc.RegularMesh()
    mesh.dimension = (10, 10)
    mesh.lower_left = (-lat_size, -lat_size)
    mesh.upper_right = (lat_size, lat_size)
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally(tally_id=1)
    tally.filters = [mesh_filter]
    tally.scores = ['flux']
    tally.estimator = 'tracklength'
    model.tallies = [tally]

    return model


def test_lattice_corner_crossing(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
