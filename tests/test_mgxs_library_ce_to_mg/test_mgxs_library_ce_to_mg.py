#!/usr/bin/env python

import os
import sys
import glob
import hashlib
sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
from input_set import PinCellInputSet
import openmc
import openmc.mgxs


class MGXSTestHarness(PyAPITestHarness):
    def _build_inputs(self):
        # Set the input set to use the pincell model
        self._input_set = PinCellInputSet()

        # Generate inputs using parent class routine
        super(MGXSTestHarness, self)._build_inputs()

        # Initialize a two-group structure
        energy_groups = openmc.mgxs.EnergyGroups(group_edges=[0, 0.625e-6,
                                                              20.])

        # Initialize MGXS Library for a few cross section types
        self.mgxs_lib = openmc.mgxs.Library(self._input_set.geometry)
        self.mgxs_lib.by_nuclide = False
        self.mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission matrix',
                                    'nu-scatter matrix', 'multiplicity matrix']
        self.mgxs_lib.energy_groups = energy_groups
        self.mgxs_lib.correction = None
        self.mgxs_lib.legendre_order = 3
        self.mgxs_lib.domain_type = 'material'
        self.mgxs_lib.build_library()

        # Initialize a tallies file
        self._input_set.tallies = openmc.Tallies()
        self.mgxs_lib.add_to_tallies_file(self._input_set.tallies, merge=False)
        self._input_set.tallies.export_to_xml()

    def _run_openmc(self):
        # Initial run
        if self._opts.mpi_exec is not None:
            returncode = openmc.run(mpi_procs=self._opts.mpi_np,
                                    openmc_exec=self._opts.exe,
                                    mpi_exec=self._opts.mpi_exec)

        else:
            returncode = openmc.run(openmc_exec=self._opts.exe)

        assert returncode == 0, 'CE OpenMC calculation did not exit' \
                                'successfully.'

        # Build MG Inputs
        # Get data needed to execute Library calculations.
        statepoint = glob.glob(os.path.join(os.getcwd(), self._sp_name))[0]
        sp = openmc.StatePoint(statepoint)
        self.mgxs_lib.load_from_statepoint(sp)
        self._input_set.mgxs_file, self._input_set.materials, \
            self._input_set.geometry = self.mgxs_lib.create_mg_mode()

        # Modify settings so we can run in MG mode
        self._input_set.settings.cross_sections = './mgxs.xml'
        self._input_set.settings.energy_mode = 'multi-group'

        # Write modified input files
        self._input_set.settings.export_to_xml()
        self._input_set.geometry.export_to_xml()
        self._input_set.materials.export_to_xml()
        self._input_set.mgxs_file.export_to_xml()
        # Dont need tallies.xml, so remove the file
        if os.path.exists('./tallies.xml'):
            os.remove('./tallies.xml')

        # Re-run MG mode.
        if self._opts.mpi_exec is not None:
            returncode = openmc.run(mpi_procs=self._opts.mpi_np,
                                    openmc_exec=self._opts.exe,
                                    mpi_exec=self._opts.mpi_exec)

        else:
            returncode = openmc.run(openmc_exec=self._opts.exe)

    def _cleanup(self):
        super(MGXSTestHarness, self)._cleanup()
        f = os.path.join(os.getcwd(), 'mgxs.xml')
        if os.path.exists(f):
            os.remove(f)


if __name__ == '__main__':
    harness = MGXSTestHarness('statepoint.10.*', False)
    harness.main()
