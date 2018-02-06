import hashlib

import openmc
import openmc.mgxs

from tests.testing_harness import PyAPITestHarness


class MGXSTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Initialize a one-group structure
        energy_groups = openmc.mgxs.EnergyGroups(group_edges=[0, 20.e6])

        # Initialize MGXS Library for a few cross section types
        # for one material-filled cell in the geometry
        self.mgxs_lib = openmc.mgxs.Library(self._model.geometry)
        self.mgxs_lib.by_nuclide = False

        # Test all MGXS types
        self.mgxs_lib.mgxs_types = openmc.mgxs.MGXS_TYPES + \
                                   openmc.mgxs.MDGXS_TYPES
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.num_delayed_groups = 6
        self.mgxs_lib.legendre_order = 3
        self.mgxs_lib.domain_type = 'mesh'

        # Instantiate a tally mesh
        mesh = openmc.Mesh(mesh_id=1)
        mesh.type = 'regular'
        mesh.dimension = [2, 2]
        mesh.lower_left = [-100., -100.]
        mesh.width = [100., 100.]

        self.mgxs_lib.domains = [mesh]
        self.mgxs_lib.build_library()

        # Add tallies
        self.mgxs_lib.add_to_tallies_file(self._model.tallies, merge=False)

    def _get_results(self, hash_output=False):
        """Digest info in the statepoint and return as a string."""

        # Read the statepoint file.
        sp = openmc.StatePoint(self._sp_name)

        # Load the MGXS library from the statepoint
        self.mgxs_lib.load_from_statepoint(sp)

        # Build a string from Pandas Dataframe for each 1-group MGXS
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


def test_mgxs_library_mesh():
    harness = MGXSTestHarness('statepoint.10.h5')
    harness.main()
