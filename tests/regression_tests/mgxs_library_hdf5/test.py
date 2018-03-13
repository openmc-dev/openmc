import os
import hashlib

import numpy as np
import h5py
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
        self.mgxs_lib.mgxs_types = openmc.mgxs.MGXS_TYPES + \
                                   openmc.mgxs.MDGXS_TYPES
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.num_delayed_groups = 6
        self.mgxs_lib.legendre_order = 3
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

        # Export the MGXS Library to an HDF5 file
        self.mgxs_lib.build_hdf5_store(directory='.')

        # Open the MGXS HDF5 file
        with h5py.File('mgxs.h5', 'r') as f:

            # Build a string from the datasets in the HDF5 file
            outstr = ''
            for domain in self.mgxs_lib.domains:
                for mgxs_type in self.mgxs_lib.mgxs_types:
                    outstr += 'domain={0} type={1}\n'.format(domain.id, mgxs_type)
                    avg_key = 'material/{0}/{1}/average'.format(domain.id, mgxs_type)
                    std_key = 'material/{0}/{1}/std. dev.'.format(domain.id, mgxs_type)
                    outstr += '{}\n{}\n'.format(f[avg_key][...], f[std_key][...])

        # Hash the results if necessary
        if hash_output:
            sha512 = hashlib.sha512()
            sha512.update(outstr.encode('utf-8'))
            outstr = sha512.hexdigest()

        return outstr

    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mgxs_library_hdf5():
    try:
        np.set_printoptions(formatter={'float_kind': '{:.8e}'.format})
        model = pwr_pin_cell()
        harness = MGXSTestHarness('statepoint.10.h5', model)
        harness.main()
    finally:
        np.set_printoptions(formatter=None)
