"""Regression test for rectangular lattices with large pitch values.

Large pitches used to trigger a segmentation fault in ``RectLattice::distance``
because the boundary-crossing check compared a reconstructed crossing position
against an absolute tolerance that did not scale with the geometry (see #3852).
This test transports particles across a lattice with a large pitch to ensure the
crossing logic remains robust.

"""

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.Model()

    # Large pitch that previously caused floating-point cancellation
    pitch = 100.0
    n = 3

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

    lattice = openmc.RectLattice()
    lattice.lower_left = (-pitch*n/2, -pitch*n/2)
    lattice.pitch = (pitch, pitch)
    lattice.outer = air_uni
    lattice.universes = [
        [metal_uni, air_uni, metal_uni],
        [air_uni, metal_uni, air_uni],
        [metal_uni, air_uni, metal_uni],
    ]

    box = openmc.model.RectangularPrism(pitch*n, pitch*n, boundary_type='vacuum')
    root_cell = openmc.Cell(region=-box, fill=lattice)
    model.geometry = openmc.Geometry([root_cell])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 1000

    mesh = openmc.RegularMesh()
    mesh.dimension = (6, 6)
    mesh.lower_left = (-pitch*n/2, -pitch*n/2)
    mesh.upper_right = (pitch*n/2, pitch*n/2)
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally(tally_id=1)
    tally.filters = [mesh_filter]
    tally.scores = ['flux']
    model.tallies = [tally]

    return model


def test_lattice_large_pitch(model):
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
