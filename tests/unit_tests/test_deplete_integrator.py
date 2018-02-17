"""Tests for integrator.py

It is worth noting that openmc.deplete.integrate is extremely complex, to the
point I am unsure if it can be reasonably unit-tested.  For the time being, it
will be left unimplemented and testing will be done via regression.

"""

import copy
import os
from unittest.mock import MagicMock

import numpy as np
from openmc.deplete import (integrator, ReactionRates, results, comm,
                            OperatorResult)


def test_save_results(run_in_tmpdir):
    """Test data save module"""

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
    cell_dict = {s: i for i, s in enumerate(burn_list)}
    r1 = ReactionRates(cell_dict, {"na": 0, "nb": 1}, {"ra": 0, "rb": 1})
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

    eigvl1 = comm.bcast(eigvl1, root=0)
    eigvl2 = comm.bcast(eigvl2, root=0)

    t1 = [0.0, 1.0]
    t2 = [1.0, 2.0]

    op_result1 = [OperatorResult(k, rates) for k, rates in zip(eigvl1, rate1)]
    op_result2 = [OperatorResult(k, rates) for k, rates in zip(eigvl2, rate2)]
    integrator.save_results(op, x1, op_result1, t1, 0)
    integrator.save_results(op, x2, op_result2, t2, 1)

    # Load the files
    res = results.read_results("depletion_results.h5")

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
