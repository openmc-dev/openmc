import glob
import os
import copy

from tests.testing_harness import CMFDTestHarness
from openmc import cmfd
import numpy as np


class CMFDRestartTestHarness(CMFDTestHarness):
    def __init__(self, final_sp, restart_sp, cmfd_run1, cmfd_run2):
        super().__init__(final_sp, cmfd_run1)
        self._cmfd_restart_run = cmfd_run2
        self._restart_sp = restart_sp

    def execute_test(self):
        try:
            # Compare results from first CMFD run
            self._test_output_created()
            results = self._get_results()
            results += self._cmfdrun_results
            self._write_results(results)
            self._compare_results()

            # Run CMFD from restart file
            statepoint = glob.glob(os.path.join(os.getcwd(), self._restart_sp))
            assert len(statepoint) == 1
            statepoint = statepoint[0]
            self._cmfd_restart_run.run(args=['-r', statepoint])

            # Compare results from second CMFD run
            self._test_output_created()
            self._create_cmfd_result_str(self._cmfd_restart_run)
            results = self._get_results()
            results += self._cmfdrun_results
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()


def test_cmfd_restart():
    """Test 1 group CMFD solver with restart run"""
    # Initialize and set CMFD mesh, create a copy for second run
    cmfd_mesh = cmfd.CMFDMesh()
    cmfd_mesh.lower_left = (-10.0, -1.0, -1.0)
    cmfd_mesh.upper_right = (10.0, 1.0, 1.0)
    cmfd_mesh.dimension = (10, 1, 1)
    cmfd_mesh.albedo = (0.0, 0.0, 1.0, 1.0, 1.0, 1.0)
    cmfd_mesh2 = copy.deepcopy(cmfd_mesh)

    # Initialize and run first CMFDRun object
    cmfd_run = cmfd.CMFDRun()
    cmfd_run.mesh = cmfd_mesh
    cmfd_run.tally_begin = 5
    cmfd_run.solver_begin = 5
    cmfd_run.feedback = True
    cmfd_run.gauss_seidel_tolerance = [1.e-15, 1.e-20]
    cmfd_run.run()

    # Initialize second CMFDRun object which will be run from restart file
    cmfd_run2 = cmfd.CMFDRun()
    cmfd_run2.mesh = cmfd_mesh2
    cmfd_run2.tally_begin = 5
    cmfd_run2.solver_begin = 5
    cmfd_run2.feedback = True
    cmfd_run2.gauss_seidel_tolerance = [1.e-15, 1.e-20]

    # Initialize and run CMFD restart test harness
    harness = CMFDRestartTestHarness('statepoint.20.h5', 'statepoint.15.h5',
                                     cmfd_run, cmfd_run2)
    harness.main()
