""" Full system test suite. """

import shutil
import unittest

import numpy as np

import opendeplete
from opendeplete import results
from opendeplete import utilities
import test.example_geometry as example_geometry


class TestFull(unittest.TestCase):
    """ Full system test suite.

    Runs an entire OpenMC simulation with depletion coupling and verifies
    that the outputs match a reference file.  Sensitive to changes in
    OpenMC.
    """

    def test_full(self):
        """
        This test runs a complete OpenMC simulation and tests the outputs.
        It will take a while.
        """

        n_rings = 2
        n_wedges = 4

        # Load geometry from example
        geometry, lower_left, upper_right = \
            example_geometry.generate_problem(n_rings=n_rings, n_wedges=n_wedges)

        # Create dt vector for 3 steps with 15 day timesteps
        dt1 = 15*24*60*60  # 15 days
        dt2 = 1.5*30*24*60*60  # 1.5 months
        N = np.floor(dt2/dt1)

        dt = np.repeat([dt1], N)

        # Create settings variable
        settings = opendeplete.OpenMCSettings()

        settings.chain_file = "chains/chain_simple.xml"
        settings.openmc_call = "openmc"
        settings.openmc_npernode = 2
        settings.particles = 100
        settings.batches = 100
        settings.inactive = 40
        settings.lower_left = lower_left
        settings.upper_right = upper_right
        settings.entropy_dimension = [10, 10, 1]

        settings.round_number = True
        settings.constant_seed = 1

        joule_per_mev = 1.6021766208e-13
        settings.power = 2.337e15*4*joule_per_mev  # MeV/second cm from CASMO
        settings.dt_vec = dt
        settings.output_dir = "test_full"

        op = opendeplete.OpenMCOperator(geometry, settings)

        # Perform simulation using the predictor algorithm
        opendeplete.integrator.predictor(op)

        # Load the files
        res_test = results.read_results(settings.output_dir + "/results.h5")

        # Load the reference
        res_old = results.read_results("test/test_reference.h5")

        # Assert same mats
        for mat in res_old[0].mat_to_ind:
            self.assertIn(mat, res_test[0].mat_to_ind,
                          msg="Cell " + mat + " not in new results.")
        for nuc in res_old[0].nuc_to_ind:
            self.assertIn(nuc, res_test[0].nuc_to_ind,
                          msg="Nuclide " + nuc + " not in new results.")

        for mat in res_test[0].mat_to_ind:
            self.assertIn(mat, res_old[0].mat_to_ind,
                          msg="Cell " + mat + " not in old results.")
        for nuc in res_test[0].nuc_to_ind:
            self.assertIn(nuc, res_old[0].nuc_to_ind,
                          msg="Nuclide " + nuc + " not in old results.")

        for mat in res_test[0].mat_to_ind:
            for nuc in res_test[0].nuc_to_ind:
                _, y_test = utilities.evaluate_single_nuclide(res_test, mat, nuc)
                _, y_old = utilities.evaluate_single_nuclide(res_old, mat, nuc)

                # Test each point

                tol = 1.0e-6

                correct = True
                for i, ref in enumerate(y_old):
                    if ref != y_test[i]:
                        if ref != 0.0:
                            if np.abs(y_test[i] - ref) / ref > tol:
                                correct = False
                        else:
                            correct = False

                self.assertTrue(correct,
                                msg="Discrepancy in mat " + mat + " and nuc " + nuc
                                + "\n" + str(y_old) + "\n" + str(y_test))

    def tearDown(self):
        """ Clean up files"""
        opendeplete.comm.barrier()
        if opendeplete.comm.rank == 0:
            shutil.rmtree("test_full", ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
