""" Tests for integrator.py """

import copy
import os
import unittest
from unittest.mock import MagicMock

import numpy as np

from opendeplete import integrator, ReactionRates, results, comm


class TestIntegrator(unittest.TestCase):
    """ Tests for integrator.py

    It is worth noting that opendeplete.integrate is extremely complex, to
    the point I am unsure if it can be reasonably unit-tested.  For the time
    being, it will be left unimplemented and testing will be done via
    regression (in test_integrator_regression.py)
    """

    def test_save_results(self):
        """ Test data save module """

        stages = 3

        np.random.seed(comm.rank)

        # Mock geometry
        op = MagicMock()

        vol_dict = {}
        full_burn_dict = {}

        j = 0
        for i in range(comm.size):
            vol_dict[str(2*i)] = 1.2
            vol_dict[str(2*i + 1)] = 1.2
            full_burn_dict[str(2*i)] = j
            full_burn_dict[str(2*i + 1)] = j + 1
            j += 2

        burn_list = [str(i) for i in range(2*comm.rank, 2*comm.rank + 2)]
        nuc_list = ["na", "nb"]

        op.get_results_info.return_value = vol_dict, nuc_list, burn_list, full_burn_dict

        # Construct x
        x1 = []
        x2 = []

        for i in range(stages):
            x1.append([np.random.rand(2), np.random.rand(2)])
            x2.append([np.random.rand(2), np.random.rand(2)])

        # Construct r
        cell_dict = {s:i for i, s in enumerate(burn_list)}
        r1 = ReactionRates(cell_dict, {"na":0, "nb":1}, {"ra":0, "rb":1})
        r1.rates = np.random.rand(2, 2, 2)

        rate1 = []
        rate2 = []

        for i in range(stages):
            rate1.append(copy.deepcopy(r1))
            r1.rates = np.random.rand(2, 2, 2)
            rate2.append(copy.deepcopy(r1))
            r1.rates = np.random.rand(2, 2, 2)

        # Create global terms
        eigvl1 = np.random.rand(stages)
        eigvl2 = np.random.rand(stages)
        seed1 = [np.random.randint(100) for i in range(stages)]
        seed2 = [np.random.randint(100) for i in range(stages)]

        eigvl1 = comm.bcast(eigvl1, root=0)
        eigvl2 = comm.bcast(eigvl2, root=0)
        seed1 = comm.bcast(seed1, root=0)
        seed2 = comm.bcast(seed2, root=0)

        t1 = [0.0, 1.0]
        t2 = [1.0, 2.0]

        integrator.save_results(op, x1, rate1, eigvl1, seed1, t1, 0)
        integrator.save_results(op, x2, rate2, eigvl2, seed2, t2, 1)

        # Load the files
        res = results.read_results("results.h5")

        for i in range(stages):
            for mat_i, mat in enumerate(burn_list):

                for nuc_i, nuc in enumerate(nuc_list):
                    self.assertEqual(res[0][i, mat, nuc], x1[i][mat_i][nuc_i])
                    self.assertEqual(res[1][i, mat, nuc], x2[i][mat_i][nuc_i])
                    np.testing.assert_array_equal(res[0].rates[i][mat, nuc, :],
                                                  rate1[i][mat, nuc, :])
                    np.testing.assert_array_equal(res[1].rates[i][mat, nuc, :],
                                                  rate2[i][mat, nuc, :])

        np.testing.assert_array_equal(res[0].k, eigvl1)
        np.testing.assert_array_equal(res[0].seeds, seed1)
        np.testing.assert_array_equal(res[0].time, t1)

        np.testing.assert_array_equal(res[1].k, eigvl2)
        np.testing.assert_array_equal(res[1].seeds, seed2)
        np.testing.assert_array_equal(res[1].time, t2)

        # Delete files
        comm.barrier()
        if comm.rank == 0:
            os.remove("results.h5")


if __name__ == '__main__':
    unittest.main()
