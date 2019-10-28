import hashlib

import openmc
import openmc.mgxs
from openmc.examples import pwr_pin_cell

from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        # Generate inputs using parent class routine
        super().__init__(*args, **kwargs)

        # Initialize a two-group structure
        energy_groups = openmc.mgxs.EnergyGroups(group_edges=[0, 0.625, 20.e6])

        # Initialize MGXS Library for a few cross section types
        self.mgxs_lib = openmc.mgxs.Library(self._model.geometry)
        self.mgxs_lib.by_nuclide = False

        # Test all MGXS types
        self.mgxs_lib.mgxs_types = ['scatter matrix', 'nu-scatter matrix',
                                    'consistent scatter matrix',
                                    'consistent nu-scatter matrix']
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.correction = 'P0'
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.build_library()

        # Add tallies
        self.mgxs_lib.add_to_tallies_file(self._model.tallies, merge=False)

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Load the MGXS library from the statepoint
        self.mgxs_lib.load_from_statepoint(sp)

        # Build a string from Pandas Dataframe for each MGXS
        outstr = ''
        for domain in self.mgxs_lib.domains:
            for mgxs_type in self.mgxs_lib.mgxs_types:
                mgxs = self.mgxs_lib.get_mgxs(domain, mgxs_type)
                df = mgxs.get_pandas_dataframe()
                outstr += df.to_string() + '\n'

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr


def test_mgxs_library_correction():
    model = pwr_pin_cell()
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
