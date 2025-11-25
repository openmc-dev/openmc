import os
import glob

import openmc
import openmc.examples


def test_statepoint_latest_pyapi():
    """Test that negative batch numbers keep running statepoints for last N batches."""
    # Use the PWR pin cell example model which is known to work
    model = openmc.examples.pwr_pin_cell()
    
    # Set batches: need at least 1 inactive + 1 active
    model.settings.inactive = 1
    model.settings.batches = 5
    model.settings.particles = 1000
    
    # Enable running statepoints: keep last 2 batches
    model.settings.statepoint = {'batches': -2}
    
    # Run the model
    sp_path = model.run(output=False)
    
    # Check that exactly 2 running statepoint files exist
    running = sorted(glob.glob('statepoint.running.*.h5'))
    assert len(running) == 2, f'Expected 2 running statepoint files, found {len(running)}: {running}'
    
    # Verify the batch numbers are 4 and 5 (last 2 of 5 batches)
    batch_nums = sorted([int(f.split('.')[2]) for f in running])
    assert batch_nums == [4, 5], f'Expected batches [4, 5], found {batch_nums}'
    
    # Clean up the running statepoint files
    for f in running:
        if os.path.exists(f):
            os.remove(f)
