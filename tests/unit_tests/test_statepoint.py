import h5py
import numpy as np
import pytest

import openmc


def test_get_tally_filter_type(run_in_tmpdir):
    """Test various ways of retrieving tallies from a StatePoint object."""

    mat = openmc.Material()
    mat.add_nuclide("H1", 1.0)
    mat.set_density("g/cm3", 10.0)

    sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=mat, region=-sphere)
    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.particles = 10
    settings.batches = 2
    settings.run_mode = "fixed source"

    reg_mesh = openmc.RegularMesh().from_domain(cell)
    tally1 = openmc.Tally(tally_id=1)
    mesh_filter = openmc.MeshFilter(reg_mesh)
    tally1.filters = [mesh_filter]
    tally1.scores = ["flux"]

    tally2 = openmc.Tally(tally_id=2, name="heating tally")
    cell_filter = openmc.CellFilter(cell)
    tally2.filters = [cell_filter]
    tally2.scores = ["heating"]

    tallies = openmc.Tallies([tally1, tally2])
    model = openmc.Model(
        geometry=geometry, materials=[mat], settings=settings, tallies=tallies
    )

    sp_filename = model.run()

    sp = openmc.StatePoint(sp_filename)
    assert 'path' not in sp._f.attrs
    assert sp.path == ''

    tally_found = sp.get_tally(filter_type=openmc.MeshFilter)
    assert tally_found.id == 1

    tally_found = sp.get_tally(filter_type=openmc.CellFilter)
    assert tally_found.id == 2

    tally_found = sp.get_tally(filters=[mesh_filter])
    assert tally_found.id == 1

    tally_found = sp.get_tally(filters=[cell_filter])
    assert tally_found.id == 2

    tally_found = sp.get_tally(scores=["heating"])
    assert tally_found.id == 2

    tally_found = sp.get_tally(name="heating tally")
    assert tally_found.id == 2

    tally_found = sp.get_tally(name=None)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=1)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=2)
    assert tally_found.id == 2


class _FakeStatePoint:
    """Stand-in exposing just what multiplication/source_efficiency use"""
    run_mode = 'eigenvalue'
    _k_generation_inactive = openmc.StatePoint._k_generation_inactive

    def __init__(self, k_generation, n_inactive, generations_per_batch, keff):
        self.k_generation = np.asarray(k_generation)
        self.n_inactive = n_inactive
        self.generations_per_batch = generations_per_batch
        self.keff = keff


def _multiplication(sp):
    return openmc.StatePoint.multiplication.fget(sp)


def _source_efficiency(sp):
    return openmc.StatePoint.source_efficiency.fget(sp)


def test_multiplication_no_inactive():
    sp = _FakeStatePoint([0.9, 0.8], 0, 1, 0.8)
    assert _multiplication(sp) == 1.0
    assert _source_efficiency(sp) == 1.0


def test_multiplication_generations_per_batch():
    # Two inactive batches of two generations each use four k estimates
    k_gen = [0.5, 0.5, 0.5, 0.5, 0.1, 0.1]
    keff = 0.5
    sp = _FakeStatePoint(k_gen, 2, 2, keff)
    expected = 1.0 + 0.5 + 0.25 + 0.125 + 0.0625/(1.0 - keff)
    assert _multiplication(sp) == pytest.approx(expected)
    assert _source_efficiency(sp) == pytest.approx(1.0)


def _write_minimal_statepoint(path, version, global_tallies, k_generation):
    with h5py.File(path, 'w') as f:
        f.attrs['filetype'] = np.bytes_('statepoint')
        f.attrs['version'] = np.array(version)
        f['run_mode'] = np.bytes_('eigenvalue')
        f['n_realizations'] = 4
        f['n_inactive'] = 2
        f['generations_per_batch'] = 1
        f['k_combined'] = [0.5, 0.01]
        f['global_tallies'] = global_tallies
        f['k_generation'] = k_generation


@pytest.mark.parametrize('version', [(18, 2), (18, 3)])
def test_global_tallies_layout(run_in_tmpdir, version):
    """Old files store (value, sum, sum_sq); 18.3+ store (sum, sum_sq)"""
    sums = np.array([4.0, 4.4, 3.6, 0.8])
    sums_sq = np.array([4.1, 4.9, 3.3, 0.2])
    if version < (18, 3):
        data = np.column_stack([np.full(4, 99.0), sums, sums_sq])
        k_gen = np.array([0.9, 0.95, 1.0, 1.0, 1.0, 1.0])
    else:
        data = np.column_stack([sums, sums_sq])
        k_gen = np.column_stack([[0.9, 0.95, 1.0, 1.0, 1.0, 1.0],
                                 np.full(6, 0.01)])
    _write_minimal_statepoint('statepoint.h5', version, data, k_gen)

    with openmc.StatePoint('statepoint.h5', autolink=False) as sp:
        gt = sp.global_tallies
        np.testing.assert_allclose(gt['sum'], sums)
        np.testing.assert_allclose(gt['sum_sq'], sums_sq)
        np.testing.assert_allclose(gt['mean'], sums/4)
        assert sp.keff.n == pytest.approx(0.5)
        assert len(sp.k_generation) == 6
        assert sp.source_efficiency.n == pytest.approx(0.9*0.95/0.5**2)
        expected = 1.0 + 0.9 + 0.9*0.95/(1.0 - 0.5)
        assert sp.multiplication.n == pytest.approx(expected)
