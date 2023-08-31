import hashlib

import openmc
import openmc.mgxs
from openmc.examples import pwr_pin_cell

from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Initialize a two-group structure
        energy_groups = openmc.mgxs.EnergyGroups(group_edges=[0, 0.625, 20.e6])

        # Initialize MGXS Library for a few cross section types
        self.mgxs_lib = openmc.mgxs.Library(self._model.geometry)
        self.mgxs_lib.by_nuclide = True

        # Test relevant MGXS types
        relevant_MGXS_TYPES = [item for item in openmc.mgxs.MGXS_TYPES
                               if item != 'current']
        # Add in a subset of openmc.mgxs.ARBITRARY_VECTOR_TYPES and
        # openmc.mgxs.ARBITRARY_MATRIX_TYPES so we can see the code works,
        # but not use too much resources
        relevant_MGXS_TYPES += [
            "(n,elastic)", "(n,level)", "(n,2n)", "(n,na)", "(n,nc)",
            "(n,gamma)", "(n,a)", "(n,Xa)", "heating", "damage-energy",
            "(n,n1)", "(n,a0)", "(n,nc) matrix", "(n,n1) matrix",
            "(n,2n) matrix"]
        self.mgxs_lib.mgxs_types = tuple(relevant_MGXS_TYPES)
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.legendre_order = 3
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.nuclides = ['U235', 'Zr90', 'H1']
        self.mgxs_lib.build_library()

        # Add tallies
        self.mgxs_lib.add_to_tallies_file(self._model.tallies, merge=True)

    def _get_results(self, hash_output=True):
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


def test_mgxs_library_specific_nuclides():
    model = pwr_pin_cell()
    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
