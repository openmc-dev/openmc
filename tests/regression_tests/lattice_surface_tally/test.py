"""Test that surface tallies work correctly in lattice geometries.

This is a regression test for a bug introduced in PR #3742 where the surface
reference in event_cross_surface() was dereferenced unconditionally, including
for lattice boundary crossings where surface() = SURFACE_NONE (0). This caused
an out-of-bounds access on model::surfaces when surface tallies were active.
"""

import openmc

from tests.testing_harness import PyAPITestHarness


def model():
    """Simple 3x3 lattice of UO2 with a surface tally."""
    mat = openmc.Material()
    mat.add_nuclide('U235', 0.04)
    mat.add_nuclide('U238', 0.96)
    mat.add_nuclide('O16', 2.0)
    mat.set_density('g/cm3', 10.0)

    # Create a 3x3 lattice of infinite cells filled with UO2
    pitch = 2.0
    nx, ny = 3, 3
    universes = []
    for _ in range(nx * ny):
        cell = openmc.Cell(fill=mat)
        universes.append(openmc.Universe(cells=[cell]))

    lattice = openmc.RectLattice()
    lattice.lower_left = [-nx * pitch / 2, -ny * pitch / 2]
    lattice.pitch = [pitch, pitch]
    lattice.universes = [
        [universes[6], universes[7], universes[8]],
        [universes[3], universes[4], universes[5]],
        [universes[0], universes[1], universes[2]],
    ]

    xmin = openmc.XPlane(-nx * pitch / 2, boundary_type='vacuum')
    xmax = openmc.XPlane(nx * pitch / 2, boundary_type='vacuum')
    ymin = openmc.YPlane(-ny * pitch / 2, boundary_type='vacuum')
    ymax = openmc.YPlane(ny * pitch / 2, boundary_type='vacuum')

    root_cell = openmc.Cell(fill=lattice, region=+xmin & -xmax & +ymin & -ymax)

    model = openmc.Model()
    model.geometry = openmc.Geometry([root_cell])
    model.materials = openmc.Materials([mat])

    model.settings.particles = 1000
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box([-1, -1, -1], [1, 1, 1]),
        constraints={'fissionable': True}
    )

    # Surface tally — this triggers the bug because score_surface_tally is
    # called for all events (including lattice crossings) when
    # active_surface_tallies is non-empty
    surf_filter = openmc.SurfaceFilter([xmin, xmax, ymin, ymax])
    tally = openmc.Tally()
    tally.filters = [surf_filter]
    tally.scores = ['current']
    model.tallies = [tally]

    return model


def test_lattice_surface_tally():
    harness = PyAPITestHarness('statepoint.5.h5', model=model())
    harness.main()
