import os

import openmc
import openmc.mgxs
from openmc.examples import pwr_pin_cell

from tests.testing_harness import PyAPITestHarness
from tests.regression_tests import config


class MGXSTestHarness(PyAPITestHarness):
    def __init__(self, *args, **kwargs):
        # Generate inputs using parent class routine
        super().__init__(*args, **kwargs)

        # Initialize a two-group structure
        energy_groups = openmc.mgxs.EnergyGroups(group_edges=[0, 0.625, 20.e6])

        # Initialize MGXS Library for a few cross section types
        self.mgxs_lib = openmc.mgxs.Library(self._model.geometry)
        self.mgxs_lib.by_nuclide = True
        self.mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission matrix',
                                    'nu-scatter matrix', 'multiplicity matrix']
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.correction = None
        self.mgxs_lib.legendre_order = 3
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.build_library()

        # Initialize a tallies file
        self.mgxs_lib.add_to_tallies_file(self._model.tallies, merge=False)

    def _run_openmc(self):
        # Initial run
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(openmc_exec=config['exe'], mpi_args=mpi_args)
        else:
            openmc.run(openmc_exec=config['exe'])

        # Build MG Inputs
        # Get data needed to execute Library calculations.
        with openmc.StatePoint(self._sp_name) as sp:
            self.mgxs_lib.load_from_statepoint(sp)
        self._model.mgxs_file, self._model.materials, \
            self._model.geometry = self.mgxs_lib.create_mg_mode()

        # Modify materials and settings so we can run in MG mode
        self._model.materials.cross_sections = './mgxs.h5'
        self._model.settings.energy_mode = 'multi-group'
        # Dont need tallies so clear them from the model
        self._model.tallies = openmc.Tallies()

        # Write modified input files
        self._model.export_to_model_xml()
        self._model.mgxs_file.export_to_hdf5()

        # Re-run MG mode.
        if config['mpi']:
            mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
            openmc.run(openmc_exec=config['exe'], mpi_args=mpi_args)
        else:
            openmc.run(openmc_exec=config['exe'])

    def _cleanup(self):
        super()._cleanup()
        f = 'mgxs.h5'
        if os.path.exists(f):
            os.remove(f)


def test_mgxs_library_ce_to_mg():
    # Set the input set to use the pincell model
    model = pwr_pin_cell()

    harness = MGXSTestHarness('statepoint.10.h5', model)
    harness.main()
