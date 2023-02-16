import openmc

import pytest

def test_restart(run_in_tmpdir):

    pincell = openmc.examples.pwr_pin_cell()

    # run the pincell
    sp_file = pincell.run()

    # run a restart with the resulting statepoint
    # and the settings unchanged

    with pytest.raises(RuntimeError, match='is smaller than the number of batches'):
        pincell.run(restart_file=sp_file)

    # update the number of batches and run again
    pincell.settings.batches = 15
    sp_file = pincell.run(restart_file=sp_file)

    sp = openmc.StatePoint(sp_file)
    assert sp.n_batches == 15