import hashlib

import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material(name='UO2')
    fuel.set_density('g/cm3', 10.29769)
    fuel.add_nuclide("U234", 4.4843e-6)
    fuel.add_nuclide("U235", 5.5815e-4)
    fuel.add_nuclide("U238", 2.2408e-2)
    fuel.add_nuclide("O16", 4.5829e-2)
    water = openmc.Material(name='light water')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    model.materials.extend([fuel, water])

    cyl = openmc.ZCylinder(r=0.4)
    pin = openmc.model.pin([cyl], [fuel, water])
    d = 1.2
    lattice = openmc.RectLattice()
    lattice.lower_left = (-d, -d)
    lattice.pitch = (d, d)
    lattice.outer = pin
    lattice.universes = [
        [pin, pin],
        [pin, pin],
    ]
    box = openmc.model.rectangular_prism(2*d, 2*d, boundary_type='reflective')
    main_cell = openmc.Cell(fill=lattice, region=box)
    model.geometry = openmc.Geometry([main_cell])

    model.settings.batches = 10
    model.settings.inactive = 5
    model.settings.particles = 1000

    energy_filter = openmc.EnergyFilter([0.0, 0.253, 1.0e3, 1.0e6, 20.0e6])
    distrib_filter = openmc.DistribcellFilter(pin.cells[1])
    tally = openmc.Tally(name='distribcell tally')
    tally.filters = [energy_filter, distrib_filter]
    tally.scores = ['nu-fission', 'total']
    tally.nuclides = ['U234', 'U235', 'U238']
    model.tallies.append(tally)

    return model



class TallyAggregationTestHarness(PyAPITestHarness):
    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Extract the tally of interest
        tally = sp.get_tally(name='distribcell tally')

        # Perform tally aggregations across filter bins, nuclides and scores
        outstr = ''

        # Sum across all energy filter bins
        tally_sum = tally.summation(filter_type=openmc.EnergyFilter)
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all distribcell filter bins
        tally_sum = tally.summation(filter_type=openmc.DistribcellFilter)
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all nuclides
        tally_sum = tally.summation(nuclides=['U234', 'U235', 'U238'])
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        # Sum across all scores
        tally_sum = tally.summation(scores=['nu-fission', 'total'])
        outstr += ', '.join(map(str, tally_sum.mean))
        outstr += ', '.join(map(str, tally_sum.std_dev))

        return outstr


def test_tally_aggregation(model):
    harness = TallyAggregationTestHarness('statepoint.10.h5', model)
    harness.main()
