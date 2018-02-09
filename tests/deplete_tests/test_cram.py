""" Tests for cram.py """

import unittest

import numpy as np
import scipy.sparse as sp

from opendeplete.integrator import CRAM16, CRAM48

class TestCram(unittest.TestCase):
    """ Tests for cram.py

    Compares a few Mathematica matrix exponentials to CRAM16/CRAM48.
    """

    def test_CRAM16(self):
        """ Test 16-term CRAM. """
        x = np.array([1.0, 1.0])
        mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
        dt = 0.1

        z = CRAM16(mat, x, dt)

        # Solution from mathematica
        z0 = np.array((0.904837418035960, 0.576799023327476))

        tol = 1.0e-15

        self.assertLess(np.linalg.norm(z - z0), tol)

    def test_CRAM48(self):
        """ Test 48-term CRAM. """
        x = np.array([1.0, 1.0])
        mat = sp.csr_matrix([[-1.0, 0.0], [-2.0, -3.0]])
        dt = 0.1

        z = CRAM48(mat, x, dt)

        # Solution from mathematica
        z0 = np.array((0.904837418035960, 0.576799023327476))

        tol = 1.0e-15

        self.assertLess(np.linalg.norm(z - z0), tol)


if __name__ == '__main__':
    unittest.main()
