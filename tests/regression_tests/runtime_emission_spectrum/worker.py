"""Run matched control/perturbed transport, including under mpiexec."""
from copy import deepcopy
import sys

import numpy as np

import openmc.lib
from openmc.mpi import comm


def run(modify):
    banks = []
    try:
        openmc.lib.init(args=['-e'] if '--event' in sys.argv else [])
        openmc.lib.simulation_init()
        index = openmc.lib.nuclides['U235']._index
        original = openmc.lib.get_emission_spectrum(index, 18)
        openmc.lib.next_batch()
        banks.append(openmc.lib.source_bank()['E'].copy())
        if modify:
            changed = deepcopy(original)
            for table in changed['tables']:
                table['pdf'] *= np.exp(-table['energy'] / 3.e5)
            openmc.lib.set_emission_spectrum(index, 18, changed)
            readback = openmc.lib.get_emission_spectrum(index, 18)
            assert not np.array_equal(readback['tables'][0]['pdf'],
                                      original['tables'][0]['pdf'])
        for _ in range(3):
            openmc.lib.next_batch()
            banks.append(openmc.lib.source_bank()['E'].copy())
        return np.asarray(banks)
    finally:
        openmc.lib.finalize()


control = run(False)
modified = run(True)
np.testing.assert_array_equal(control[0], modified[0])
control_sum = comm.allreduce(control[1:].sum())
modified_sum = comm.allreduce(modified[1:].sum())
# Same seed, same first batch, only the native spectrum changes. A large shift
# in energies of fission sites confirms that transport sampled the new tables.
assert modified_sum < 0.6 * control_sum
