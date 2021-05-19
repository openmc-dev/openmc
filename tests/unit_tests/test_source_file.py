from random import random
import subprocess

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


def test_wrong_source_attributes(run_in_tmpdir):
    # Create a source file with animal attributes
    source_dtype = np.dtype([
        ('platypus', '<f8'),
        ('axolotl', '<f8'),
        ('narwhal', '<i4'),
    ])
    arr = np.array([(1.0, 2.0, 3), (4.0, 5.0, 6), (7.0, 8.0, 9)], dtype=source_dtype)
    with h5py.File('animal_source.h5', 'w') as fh:
        fh.attrs['filetype'] = np.string_("source")
        fh.create_dataset('source_bank', data=arr)

    # Create a simple model that uses this lovely animal source
    m = openmc.Material()
    m.add_nuclide('U235', 0.02)
    openmc.Materials([m]).export_to_xml()
    s = openmc.Sphere(r=10.0, boundary_type='vacuum')
    c = openmc.Cell(fill=m, region=-s)
    openmc.Geometry([c]).export_to_xml()
    settings =  openmc.Settings()
    settings.particles = 100
    settings.batches = 10
    settings.source = openmc.Source(filename='animal_source.h5')
    settings.export_to_xml()

    # When we run the model, it should error out with a message that includes
    # the names of the wrong attributes
    with pytest.raises(RuntimeError) as excinfo:
        openmc.run()
    assert 'platypus, axolotl, narwhal' in str(excinfo.value)
