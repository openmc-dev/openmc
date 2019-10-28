import hashlib

import numpy as np
import openmc
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def model():
    model = openmc.model.Model()

    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.0)
    fuel.add_nuclide('U234', 1.0)
    fuel.add_nuclide('U235', 4.0)
    fuel.add_nuclide('U238', 95.0)
    water = openmc.Material(name='light water')
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')
    model.materials.extend([fuel, water])

    cyl1 = openmc.ZCylinder(r=5.0)
    cyl2 = openmc.ZCylinder(r=10.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fuel, region=-cyl1)
    cell2 = openmc.Cell(fill=water, region=+cyl1 & -cyl2)
    model.geometry = openmc.Geometry([cell1, cell2])

    model.settings.batches = 5
    model.settings.inactive = 0
    model.settings.particles = 1000

    mesh = openmc.RegularMesh()
    mesh.dimension = (2, 2)
    mesh.lower_left = (-10.0, -10.0)
    mesh.upper_right = (10.0, 10.0)
    energy_filter = openmc.EnergyFilter((0.0, 10.0, 20.0e6))
    material_filter = openmc.MaterialFilter((fuel, water))
    mesh_filter = openmc.MeshFilter(mesh)

    tally = openmc.Tally(name='tally 1')
    tally.filters = [material_filter, energy_filter]
    tally.scores = ['nu-fission', 'total']
    tally.nuclides = ['U234', 'U235']
    model.tallies.append(tally)
    tally = openmc.Tally(name='tally 2')
    tally.filters = [energy_filter, mesh_filter]
    tally.scores = ['total', 'fission']
    tally.nuclides = ['U238', 'U235']
    model.tallies.append(tally)

    return model


class TallyArithmeticTestHarness(PyAPITestHarness):
    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Load the tallies
        tally_1 = sp.get_tally(name='tally 1')
        tally_2 = sp.get_tally(name='tally 2')

        # Perform all the tally arithmetic operations and output results
        output = []
        with np.printoptions(precision=5, threshold=np.inf):
            mean = (tally_1 * tally_2).mean
            output.append(str(mean[np.nonzero(mean)]))

            mean = tally_1.hybrid_product(
                tally_2, '*', 'entrywise', 'tensor', 'tensor').mean
            output.append(str(mean[np.nonzero(mean)]))

            mean = tally_1.hybrid_product(
                tally_2, '*', 'entrywise', 'entrywise', 'tensor').mean
            output.append(str(mean[np.nonzero(mean)]))

            mean = tally_1.hybrid_product(
                tally_2, '*', 'entrywise', 'tensor', 'entrywise').mean
            output.append(str(mean[np.nonzero(mean)]))

            mean = tally_1.hybrid_product(
                tally_2, '*', 'entrywise', 'entrywise', 'entrywise').mean
            output.append(str(mean[np.nonzero(mean)]))

        # Hash the results if necessary
        outstr = ''.join(output)
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


def test_tally_arithmetic(model):
    harness = TallyArithmeticTestHarness('statepoint.5.h5', model)
    harness.main()
