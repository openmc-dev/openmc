from random import random

import h5py
import numpy as np
import openmc
import pytest


def test_source_file(run_in_tmpdir):
    # write_source_file shouldn't accept non-SourceParticle items
    with pytest.raises(TypeError):
        openmc.write_source_file([1, 2, 3], 'test_source.h5')

    # Create source particles
    source = []
    n = 1000
    for i in range(n):
        source.append(openmc.SourceParticle(
            r=(random(), i, 0),
            u=(0., 0., 1.),
            E=float(n - i),
        ))

    # Create source file
    openmc.write_source_file(source, 'test_source.h5')

    # Get array of source particles from file
    with h5py.File('test_source.h5', 'r') as fh:
        filetype = fh.attrs['filetype']
        arr = fh['source_bank'][...]

    # Ensure data is consistent
    assert filetype == b'source'
    r = arr['r']
    assert np.all((r['x'] > 0.0) & (r['x'] < 1.0))
    assert np.all(r['y'] == np.arange(1000))
    assert np.all(r['z'] == 0.0)
    u = arr['u']
    assert np.all(u['x'] == 0.0)
    assert np.all(u['y'] == 0.0)
    assert np.all(u['z'] == 1.0)
    assert np.all(arr['E'] == n - np.arange(n))
    assert np.all(arr['wgt'] == 1.0)
    assert np.all(arr['delayed_group'] == 0)
    assert np.all(arr['particle'] == 0)
