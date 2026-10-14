"""Test of the fission bank birth-mesh tags in continuous-energy transport.

Two slabs of a mixture of U-235 and H-1, 10 cm thick, are separated by 10 cm
of B-10 at 50 g/cm^3. The total cross section of B-10 is at least 1.4 b below
20 MeV, so the layer is more than 40 mean free paths thick at every energy. A
neutron that scatters in the layer slows down and is absorbed: at 1 keV, its
mean free path for absorption is about 0.02 cm. A neutron born in one slab
therefore causes fissions only in that slab, and the tags of the banked
fission sites can be checked without reference results, as in the multigroup
birth_mesh test:

- with one mesh bin per slab, every tag is the bin of its own site;
- with two bins per slab along x (and two along y), some tags differ from the
  bin of the site, but none leaves the slab of the site;
- with a mesh over the left half of the problem, the sites in the right slab
  are tagged -1.

A last test runs the problem through openmc.lib and compares each tag with the
bin of the source site from which the neutron that produced the site was
started, its birth position.

"""
import numpy as np
import openmc
import openmc.lib
import pytest

from tests.regression_tests import config

FUEL_WIDTH = 10.0  # cm, thickness of each fuel slab
HALF_GAP = 5.0  # cm, half-thickness of the B-10 layer
X_MAX = HALF_GAP + FUEL_WIDTH


def make_model(mesh):
    fuel = openmc.Material(name='fuel')
    fuel.add_nuclide('U235', 2.2e-4)
    fuel.add_nuclide('H1', 6.6e-2)
    fuel.set_density('sum')
    absorber = openmc.Material(name='absorber')
    absorber.add_nuclide('B10', 1.0)
    absorber.set_density('g/cm3', 50.0)

    x0 = openmc.XPlane(-X_MAX, boundary_type='reflective')
    x1 = openmc.XPlane(-HALF_GAP)
    x2 = openmc.XPlane(HALF_GAP)
    x3 = openmc.XPlane(X_MAX, boundary_type='reflective')
    y0 = openmc.YPlane(-5.0, boundary_type='reflective')
    y1 = openmc.YPlane(5.0, boundary_type='reflective')
    z0 = openmc.ZPlane(-5.0, boundary_type='reflective')
    z1 = openmc.ZPlane(5.0, boundary_type='reflective')
    yz = +y0 & -y1 & +z0 & -z1
    cells = [openmc.Cell(fill=fuel, region=+x0 & -x1 & yz),
             openmc.Cell(fill=absorber, region=+x1 & -x2 & yz),
             openmc.Cell(fill=fuel, region=+x2 & -x3 & yz)]

    model = openmc.Model()
    model.materials = openmc.Materials([fuel, absorber])
    model.geometry = openmc.Geometry(cells)
    model.settings.run_mode = 'eigenvalue'
    model.settings.particles = 1000
    model.settings.batches = 6
    model.settings.inactive = 2
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-X_MAX, -5.0, -5.0), (X_MAX, 5.0, 5.0)),
        constraints={'fissionable': True})
    model.settings.birth_mesh = mesh
    return model


def regular_mesh(x_max, dimension):
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-X_MAX, -5.0, -5.0)
    mesh.upper_right = (x_max, 5.0, 5.0)
    mesh.dimension = dimension
    return mesh


def regular_bin(mesh, xyz):
    """0-based bin (x fastest) of each position of xyz, an array of shape
    (n, 3), in a regular mesh, or -1 for a position outside the mesh."""
    ll = np.asarray(mesh.lower_left)
    dimension = np.asarray(mesh.dimension)
    width = (np.asarray(mesh.upper_right) - ll) / dimension
    ijk = np.floor((xyz - ll) / width).astype(int)
    inside = np.all((ijk >= 0) & (ijk < dimension), axis=1)
    bins = ijk[:, 0] + dimension[0] * (ijk[:, 1] + dimension[1] * ijk[:, 2])
    return np.where(inside, bins, -1)


def run_tagged(tmp_path, mesh):
    """Run the eigenvalue problem and return the final source bank, the tags,
    the bin of each site, and whether each site is in the left slab."""
    model = make_model(mesh)
    kwargs = {'cwd': tmp_path, 'openmc_exec': config['exe'],
              'event_based': config['event']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    sp_path = model.run(**kwargs)
    with openmc.StatePoint(sp_path) as sp:
        source = sp.source
        tags = sp.birth_mesh_bin
    assert tags.dtype == np.int32
    assert tags.shape == source.shape == (model.settings.particles,)

    xyz = np.column_stack([source['r'][c] for c in 'xyz'])
    left = xyz[:, 0] < -HALF_GAP
    right = xyz[:, 0] > HALF_GAP
    # All sites are in the fuel, and both slabs hold sites
    assert np.all(left | right)
    assert np.any(left) and np.any(right)
    return tags, regular_bin(mesh, xyz), left


def test_birth_mesh_ce_one_bin_per_slab(tmp_path):
    # Bin 0 holds the left slab and bin 1 the right one
    mesh = regular_mesh(X_MAX, (2, 1, 1))
    tags, own_bin, left = run_tagged(tmp_path, mesh)
    np.testing.assert_array_equal(own_bin, np.where(left, 0, 1))
    np.testing.assert_array_equal(tags, own_bin)


def test_birth_mesh_ce_two_bins_per_slab(tmp_path):
    # Bins of 5 cm along x and y: x indices 0 and 1 in the left slab, 4 and 5
    # in the right slab
    nx = 6
    mesh = regular_mesh(X_MAX, (nx, 2, 1))
    tags, own_bin, left = run_tagged(tmp_path, mesh)
    assert tags.min() >= 0
    for in_slab, slab in ((left, {0, 1}), (~left, {nx - 2, nx - 1})):
        assert set(np.unique(own_bin[in_slab] % nx)) == slab
        assert set(np.unique(tags[in_slab] % nx)) <= slab
    assert set(np.unique(tags // nx)) == {0, 1}
    assert np.count_nonzero(tags != own_bin) > 0


def test_birth_mesh_ce_left_half(tmp_path):
    # One bin over the left half of the problem; the right slab is outside
    mesh = regular_mesh(0.0, (1, 1, 1))
    tags, _, left = run_tagged(tmp_path, mesh)
    np.testing.assert_array_equal(tags, np.where(left, 0, -1))


def test_birth_mesh_ce_birth_position(tmp_path, monkeypatch):
    if config['mpi']:
        pytest.skip('openmc.lib runs in this process, without MPI')
    mesh = regular_mesh(X_MAX, (6, 2, 1))
    model = make_model(mesh)
    monkeypatch.chdir(tmp_path)
    model.init_lib(output=False)
    try:
        openmc.lib.simulation_init()
        for _ in range(model.settings.batches - 1):
            openmc.lib.next_batch()
        # Source sites of the last batch, then the sites that they produced
        parents = openmc.lib.source_bank().copy()
        openmc.lib.next_batch()
        bank = openmc.lib.source_bank().copy()
        openmc.lib.simulation_finalize()
    finally:
        model.finalize_lib()

    # parent_id is the index of the source site from which the history of the
    # neutron that produced a site was started. That neutron was born at the
    # source site, unless it was created by an (n,2n) or similar reaction in
    # that history, where it was born at the point of the reaction. Such
    # neutrons cause few of the fissions here, of the order of one in a
    # thousand, so at least 99% of the tags are the bin of the source site,
    # while the bin of the site itself differs in many of them.
    parent = bank['parent_id']
    assert parent.min() >= 0 and parent.max() < parents.size
    tags = bank['birth_mesh_bin']
    n_sites = tags.size
    birth_bin = regular_bin(mesh, parents['r'][parent])
    own_bin = regular_bin(mesh, bank['r'])
    assert np.count_nonzero(tags != birth_bin) <= 0.01 * n_sites
    assert np.count_nonzero(tags != own_bin) > 0.1 * n_sites
