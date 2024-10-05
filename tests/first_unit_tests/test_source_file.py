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

    # Ensure sites read in are consistent
    sites = openmc.ParticleList.from_hdf5('test_source.h5')

    xs = np.array([site.r[0] for site in sites])
    ys = np.array([site.r[1] for site in sites])
    zs = np.array([site.r[2] for site in sites])
    assert np.all((xs > 0.0) & (xs < 1.0))
    assert np.all(ys == np.arange(1000))
    assert np.all(zs == 0.0)
    u = np.array([s.u for s in sites])
    assert np.all(u[..., 0] == 0.0)
    assert np.all(u[..., 1] == 0.0)
    assert np.all(u[..., 2] == 1.0)
    E = np.array([s.E for s in sites])
    assert np.all(E == n - np.arange(n))
    wgt = np.array([s.wgt for s in sites])
    assert np.all(wgt == 1.0)
    dgs = np.array([s.delayed_group for s in sites])
    assert np.all(dgs == 0)
    p_types = np.array([s.particle for s in sites])
    assert np.all(p_types == 0)

    # Ensure a ParticleList item is a SourceParticle
    site = sites[0]
    assert isinstance(site, openmc.SourceParticle)
    assert site.E == pytest.approx(n)

    # Ensure site slice read in and exported are consistent
    sites_slice = sites[:10]
    sites_slice.export_to_hdf5("test_source_slice.h5")
    sites_slice = openmc.ParticleList.from_hdf5('test_source_slice.h5')

    assert isinstance(sites_slice, openmc.ParticleList)
    assert len(sites_slice) == 10
    E = np.array([s.E for s in sites_slice])
    np.testing.assert_allclose(E, n - np.arange(10))

    # Ensure site list read in and exported are consistent
    df = sites.to_dataframe()
    sites_filtered = sites[df[df.E <= 10.0].index.tolist()]
    sites_filtered.export_to_hdf5("test_source_filtered.h5")
    sites_filtered = openmc.read_source_file('test_source_filtered.h5')

    assert isinstance(sites_filtered, openmc.ParticleList)
    assert len(sites_filtered) == 10
    E = np.array([s.E for s in sites_filtered])
    np.testing.assert_allclose(E, np.arange(10, 0, -1))


def test_wrong_source_attributes(run_in_tmpdir):
    # Create a source file with animal attributes
    source_dtype = np.dtype([
        ('platypus', '<f8'),
        ('axolotl', '<f8'),
        ('narwhal', '<i4'),
    ])
    arr = np.array([(1.0, 2.0, 3), (4.0, 5.0, 6), (7.0, 8.0, 9)], dtype=source_dtype)
    with h5py.File('animal_source.h5', 'w') as fh:
        fh.attrs['filetype'] = np.bytes_("source")
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
    settings.source = openmc.FileSource(path='animal_source.h5')
    settings.export_to_xml()

    # When we run the model, it should error out with a message that includes
    # the names of the wrong attributes
    with pytest.raises(RuntimeError) as excinfo:
        openmc.run()
    assert 'platypus, axolotl, narwhal' in str(excinfo.value)


def test_source_file_transport(run_in_tmpdir):
    # Create a source file with a single particle
    particle = openmc.SourceParticle()
    openmc.write_source_file([particle], 'source.h5')

    # Created simple model to use source file
    model = openmc.Model()
    al = openmc.Material()
    al.add_element('Al', 1.0)
    al.set_density('g/cm3', 2.7)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=al, region=-sph)
    model.geometry = openmc.Geometry([cell])
    model.settings.source = openmc.FileSource(path='source.h5')
    model.settings.particles = 10
    model.settings.batches = 3
    model.settings.run_mode = 'fixed source'

    # Try running OpenMC
    model.run()
