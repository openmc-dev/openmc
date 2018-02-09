""" Regression tests for predictor.py"""

import os
import unittest

import numpy as np

import opendeplete
from opendeplete import results
from opendeplete import utilities
import test.dummy_geometry as dummy_geometry

class TestPredictorRegression(unittest.TestCase):
    """ Regression tests for opendeplete.integrator.predictor algorithm.

    These tests integrate a simple test problem described in dummy_geometry.py.
    """

    @classmethod
    def setUpClass(cls):
        """ Save current directory in case integrator crashes."""
        cls.cwd = os.getcwd()
        cls.results = "test_integrator_regression"

    def test_predictor(self):
        """ Integral regression test of integrator algorithm using CE/CM. """

        settings = opendeplete.Settings()
        settings.dt_vec = [0.75, 0.75]
        settings.output_dir = self.results

        op = dummy_geometry.DummyGeometry(settings)

        # Perform simulation using the predictor algorithm
        opendeplete.predictor(op, print_out=False)

        # Load the files
        res = results.read_results(settings.output_dir + "/results.h5")

        _, y1 = utilities.evaluate_single_nuclide(res, "1", "1")
        _, y2 = utilities.evaluate_single_nuclide(res, "1", "2")

        # Mathematica solution
        s1 = [2.46847546272295, 0.986431226850467]
        s2 = [4.11525874568034, -0.0581692232513460]

        tol = 1.0e-13

        self.assertLess(np.absolute(y1[1] - s1[0]), tol)
        self.assertLess(np.absolute(y2[1] - s1[1]), tol)

        self.assertLess(np.absolute(y1[2] - s2[0]), tol)
        self.assertLess(np.absolute(y2[2] - s2[1]), tol)

    @classmethod
    def tearDownClass(cls):
        """ Clean up files"""

        os.chdir(cls.cwd)

        opendeplete.comm.barrier()
        if opendeplete.comm.rank == 0:
            os.remove(os.path.join(cls.results, "results.h5"))
            os.rmdir(cls.results)


if __name__ == '__main__':
    unittest.main()
