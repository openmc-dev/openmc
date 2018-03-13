"""Tests for saving results

It is worth noting that openmc.deplete.integrate is extremely complex, to the
point I am unsure if it can be reasonably unit-tested.  For the time being, it
will be left unimplemented and testing will be done via regression.

"""

import copy
import os
from unittest.mock import MagicMock

import numpy as np
from openmc.deplete import (ReactionRates, Results, ResultsList, comm,
                            OperatorResult)


def test_results_save(run_in_tmpdir):
    """Test data save module"""

    stages = 3

    np.random.seed(comm.rank)

    # Mock geometry
    op = MagicMock()

    vol_dict = {}
    full_burn_list = []

    for i in range(comm.size):
        vol_dict[str(2*i)] = 1.2
        vol_dict[str(2*i + 1)] = 1.2
        full_burn_list.append(str(2*i))
        full_burn_list.append(str(2*i + 1))

    burn_list = full_burn_list[2*comm.rank : 2*comm.rank + 2]
    nuc_list = ["na", "nb"]

    op.get_results_info.return_value = vol_dict, nuc_list, burn_list, full_burn_list

    # Construct x
    x1 = []
    x2 = []

    for i in range(stages):
        x1.append([np.random.rand(2), np.random.rand(2)])
        x2.append([np.random.rand(2), np.random.rand(2)])

    # Construct r
    r1 = ReactionRates(burn_list, ["na", "nb"], ["ra", "rb"])
    r1[:] = np.random.rand(2, 2, 2)

    rate1 = []
    rate2 = []

    for i in range(stages):
        rate1.append(copy.deepcopy(r1))
        r1[:] = np.random.rand(2, 2, 2)
        rate2.append(copy.deepcopy(r1))
        r1[:] = np.random.rand(2, 2, 2)

    # Create global terms
    eigvl1 = np.random.rand(stages)
    eigvl2 = np.random.rand(stages)

    eigvl1 = comm.bcast(eigvl1, root=0)
    eigvl2 = comm.bcast(eigvl2, root=0)

    t1 = [0.0, 1.0]
    t2 = [1.0, 2.0]

    op_result1 = [OperatorResult(k, rates) for k, rates in zip(eigvl1, rate1)]
    op_result2 = [OperatorResult(k, rates) for k, rates in zip(eigvl2, rate2)]
    Results.save(op, x1, op_result1, t1, 0)
    Results.save(op, x2, op_result2, t2, 1)

    # Load the files
    res = ResultsList("depletion_results.h5")

    for i in range(stages):
        for mat_i, mat in enumerate(burn_list):
            for nuc_i, nuc in enumerate(nuc_list):
                assert res[0][i, mat, nuc] == x1[i][mat_i][nuc_i]
                assert res[1][i, mat, nuc] == x2[i][mat_i][nuc_i]
        np.testing.assert_array_equal(res[0].rates[i], rate1[i])
        np.testing.assert_array_equal(res[1].rates[i], rate2[i])

    np.testing.assert_array_equal(res[0].k, eigvl1)
    np.testing.assert_array_equal(res[0].time, t1)

    np.testing.assert_array_equal(res[1].k, eigvl2)
    np.testing.assert_array_equal(res[1].time, t2)
