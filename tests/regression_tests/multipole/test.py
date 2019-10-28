import os

import openmc
import openmc.model
import pytest

from tests.testing_harness import TestHarness, PyAPITestHarness


def make_model():
    model = openmc.model.Model()

    # Materials
    moderator = openmc.Material(material_id=1)
    moderator.set_density('g/cc', 1.0)
    moderator.add_nuclide('H1', 2.0)
    moderator.add_nuclide('O16', 1.0)
    moderator.add_s_alpha_beta('c_H_in_H2O')

    dense_fuel = openmc.Material(material_id=2)
    dense_fuel.set_density('g/cc', 4.5)
    dense_fuel.add_nuclide('U235', 1.0)

    model.materials += [moderator, dense_fuel]

    # Geometry
    c1 = openmc.Cell(cell_id=1, fill=moderator)
    mod_univ = openmc.Universe(universe_id=1, cells=(c1,))

    r0 = openmc.ZCylinder(r=0.3)
    c11 = openmc.Cell(cell_id=11, fill=dense_fuel, region=-r0)
    c11.temperature = [500, 700, 0, 800]
    c12 = openmc.Cell(cell_id=12, fill=moderator, region=+r0)
    fuel_univ = openmc.Universe(universe_id=11, cells=(c11, c12))

    lat = openmc.RectLattice(lattice_id=101)
    lat.dimension = [2, 2]
    lat.lower_left = [-2.0, -2.0]
    lat.pitch = [2.0, 2.0]
    lat.universes = [[fuel_univ]*2]*2
    lat.outer = mod_univ

    x0 = openmc.XPlane(x0=-3.0)
    x1 = openmc.XPlane(x0=3.0)
    y0 = openmc.YPlane(y0=-3.0)
    y1 = openmc.YPlane(y0=3.0)
    for s in [x0, x1, y0, y1]:
        s.boundary_type = 'reflective'
    c101 = openmc.Cell(cell_id=101, fill=lat, region=+x0 & -x1 & +y0 & -y1)
    model.geometry.root_universe = openmc.Universe(universe_id=0, cells=(c101,))

    # Settings
    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000
    model.settings.source = openmc.Source(space=openmc.stats.Box(
        [-1, -1, -1], [1, 1, 1]))
    model.settings.temperature = {'tolerance': 1000, 'multipole': True}

    # Tallies
    tally = openmc.Tally()
    tally.nuclides = ['U235', 'O16', 'total']
    tally.scores = ['total', 'fission', '(n,gamma)', 'elastic', '(n,p)']
    model.tallies.append(tally)

    return model


class MultipoleTestHarness(PyAPITestHarness):
    def _get_results(self):
        outstr = super()._get_results()
        su = openmc.Summary('summary.h5')
        outstr += str(su.geometry.get_all_cells()[11])
        return outstr


def test_multipole():
    model = make_model()
    harness = MultipoleTestHarness('statepoint.5.h5', model)
    harness.main()
