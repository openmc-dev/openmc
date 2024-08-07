import os

import openmc
from openmc.utility_funcs import change_directory
from openmc.examples import random_ray_three_region_cube
import pytest

from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


@pytest.mark.parametrize("domain_type", ["cell", "material", "universe"])
def test_random_ray_fixed_source(domain_type):
    with change_directory(domain_type):
        openmc.reset_auto_ids()
        model = random_ray_three_region_cube()

        # Based on the parameter, we need to adjust
        # the particle source constraints
        source = model.settings.source[0]
        constraints = source.constraints

        if domain_type == 'cell':
            cells = model.geometry.get_all_cells()
            for key, cell in cells.items():
                print(cell.name)
                if cell.name == 'infinite source region':
                    constraints['domain_type'] = 'cell'
                    constraints['domain_ids'] = [cell.id]
        elif domain_type == 'material':
            materials = model.materials
            for material in materials:
                if material.name == 'source':
                    constraints['domain_type'] = 'material'
                    constraints['domain_ids'] = [material.id]
        elif domain_type == 'universe':
            universes = model.geometry.get_all_universes()
            for key, universe in universes.items():
                if universe.name == 'source universe':
                    constraints['domain_type'] = 'universe'
                    constraints['domain_ids'] = [universe.id]

        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
