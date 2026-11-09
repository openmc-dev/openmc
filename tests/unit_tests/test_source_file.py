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
    assert np.all(arr['particle'] == 2112)  # PDG number for neutron

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
    assert np.all(p_types == 2112)  # PDG number for neutron

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


def test_source_file_photon_transport(run_in_tmpdir):
    # Create a source file containing a photon. Note that photon_transport is
    # not explicitly enabled in the settings -- it should be turned on
    # automatically because the source file contains a photon.
    particle = openmc.SourceParticle(E=1.0e6, particle='photon')
    openmc.write_source_file([particle], 'photon_source.h5')

    # Create simple model to use the photon source file
    model = openmc.Model()
    al = openmc.Material()
    al.add_element('Al', 1.0)
    al.set_density('g/cm3', 2.7)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=al, region=-sph)
    model.geometry = openmc.Geometry([cell])
    model.settings.source = openmc.FileSource(path='photon_source.h5')
    model.settings.particles = 10
    model.settings.batches = 3
    model.settings.run_mode = 'fixed source'

    # Running OpenMC should succeed
    model.run()


def _rotated_universe_model(angle):
    """Plane inside a universe filled into a cell rotated by `angle` about z.

    Particles are born at (-5, 0, 0) travelling along +x and cross the plane
    at the origin. The plane's normal in its own frame is +x, but in the root
    frame it is rotated by `angle`, so for angle > 90 degrees the particle
    enters the negative half-space even though its root-frame direction has a
    positive dot product with the unrotated normal.
    """
    openmc.reset_auto_ids()

    xplane = openmc.XPlane(0.0, surface_id=7)
    inner1 = openmc.Cell(cell_id=1, region=-xplane)
    inner2 = openmc.Cell(cell_id=2, region=+xplane)
    inner_univ = openmc.Universe(cells=[inner1, inner2])

    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')
    root_cell = openmc.Cell(cell_id=3, region=-sph, fill=inner_univ)
    root_cell.rotation = (0.0, 0.0, angle)

    model = openmc.Model()
    model.geometry = openmc.Geometry([root_cell])

    src = openmc.IndependentSource()
    src.space = openmc.stats.Point((-5.0, 0.0, 0.0))
    src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 20
    model.settings.source = src
    return model, inner1, inner2


def test_surface_source_half_space_rotated_universe(run_in_tmpdir):
    """A surface source on a surface in a rotated universe starts in the
    correct half-space.

    Surface source files store unsigned surface IDs, and the half-space is
    recovered when the file is read back. That recovery compares the root-frame
    particle direction against Surface::normal(), which is expressed in the
    local frame of the universe holding the surface. With a 135 degree
    rotation the two disagree in sign, so the particle used to be started in
    the cell on the wrong side of the plane.
    """
    model, _, _ = _rotated_universe_model(135.0)
    model.settings.surf_source_write = {
        'surface_ids': [7], 'max_particles': 50}
    model.run()

    # Read the surface source back in and see which cell the particles start in
    model, inner1, inner2 = _rotated_universe_model(135.0)
    model.settings.source = openmc.FileSource('surface_source.h5')
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter([inner1, inner2])]
    tally.scores = ['flux']
    model.tallies = [tally]
    model.run(apply_tally_results=True)

    # The particle direction in the plane's own frame points into the negative
    # half-space, so all of the track length belongs to the -x cell.
    assert tally.mean.flat[0] == pytest.approx(10.0)
    assert tally.mean.flat[1] == pytest.approx(0.0)


def test_surface_source_half_space_root_universe(run_in_tmpdir):
    """The half-space is still recovered for a surface in the root universe."""
    openmc.reset_auto_ids()

    xplane = openmc.XPlane(0.0, surface_id=7)
    sph = openmc.Sphere(r=10.0, boundary_type='vacuum')

    def build():
        cell1 = openmc.Cell(cell_id=1, region=-sph & -xplane)
        cell2 = openmc.Cell(cell_id=2, region=-sph & +xplane)
        model = openmc.Model()
        model.geometry = openmc.Geometry([cell1, cell2])
        src = openmc.IndependentSource()
        src.space = openmc.stats.Point((-5.0, 0.0, 0.0))
        src.angle = openmc.stats.Monodirectional((1.0, 0.0, 0.0))
        model.settings.run_mode = 'fixed source'
        model.settings.batches = 1
        model.settings.particles = 20
        model.settings.source = src
        return model, cell1, cell2

    model, _, _ = build()
    model.settings.surf_source_write = {
        'surface_ids': [7], 'max_particles': 50}
    model.run()

    openmc.reset_auto_ids()
    model, cell1, cell2 = build()
    model.settings.source = openmc.FileSource('surface_source.h5')
    tally = openmc.Tally()
    tally.filters = [openmc.CellFilter([cell1, cell2])]
    tally.scores = ['flux']
    model.tallies = [tally]
    model.run(apply_tally_results=True)

    # Travelling along +x through the plane, the particle enters the +x cell
    assert tally.mean.flat[0] == pytest.approx(0.0)
    assert tally.mean.flat[1] == pytest.approx(10.0)
