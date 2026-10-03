import h5py
import numpy as np
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


def _write_tagged_file(path, source_present, tags=None):
    """Write a minimal state point file with an optional birth tag dataset."""
    with h5py.File(path, 'w') as f:
        f.attrs['filetype'] = np.bytes_('statepoint')
        f.attrs['version'] = np.array(
            [openmc.statepoint._VERSION_STATEPOINT, 0], dtype=np.int32)
        f.attrs['source_present'] = np.int32(source_present)
        if tags is not None:
            f.attrs['birth_mesh_id'] = np.int32(3)
        if source_present:
            source = np.zeros(4, dtype=[('E', '<f8'), ('wgt', '<f8')])
            f.create_dataset('source_bank', data=source)
            if tags is not None:
                f.create_dataset('birth_mesh_bin', data=tags)


def test_birth_mesh_tags(run_in_tmpdir):
    tags = np.array([0, 3, -1, 1], dtype=np.int32)

    # Tags stored with the source bank are returned aligned with it
    _write_tagged_file('tagged.h5', True, tags)
    with openmc.StatePoint('tagged.h5', autolink=False) as sp:
        assert sp.birth_mesh_id == 3
        assert isinstance(sp.birth_mesh_id, int)
        np.testing.assert_array_equal(sp.birth_mesh_bin, tags)
        assert sp.birth_mesh_bin.dtype == np.int32
        assert sp.birth_mesh_bin.shape == sp.source.shape

    # A tagged run whose source bank is stored in a separate file
    _write_tagged_file('separate.h5', False, tags)
    with openmc.StatePoint('separate.h5', autolink=False) as sp:
        assert sp.birth_mesh_id == 3
        assert sp.birth_mesh_bin is None

    # A run without tagging
    _write_tagged_file('untagged.h5', True)
    with openmc.StatePoint('untagged.h5', autolink=False) as sp:
        assert sp.source_present
        assert sp.birth_mesh_id is None
        assert sp.birth_mesh_bin is None
