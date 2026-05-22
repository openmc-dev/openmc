import h5py
import openmc


def test_ww_roundtrip(request, run_in_tmpdir):
    # Load weight windows from a wwinp file
    wwinp_file = request.path.with_name('wwinp_n')
    wws = openmc.WeightWindowsList.from_wwinp(wwinp_file)

    # Roundtrip them, writing to HDF5 and reading back in
    wws.export_to_hdf5('ww.h5')
    wws_new = openmc.WeightWindowsList.from_hdf5('ww.h5')

    # Check that the new weight windows are the same as the original
    assert len(wws) == len(wws_new)
    for ww, ww_new in zip(wws, wws_new):
        assert ww.particle_type == ww_new.particle_type
        assert (ww.lower_ww_bounds == ww_new.lower_ww_bounds).all()
        assert (ww.upper_ww_bounds == ww_new.upper_ww_bounds).all()
        assert ww.survival_ratio == ww_new.survival_ratio
        assert ww.num_energy_bins == ww_new.num_energy_bins
        assert ww.max_split == ww_new.max_split
        assert ww.weight_cutoff == ww_new.weight_cutoff
        assert ww.mesh.id == ww_new.mesh.id


def test_export_hdf5_format(request, run_in_tmpdir):
    # C++ openmc_weight_windows_import expects this on-disk layout.
    wws = openmc.WeightWindowsList.from_wwinp(request.path.with_name('wwinp_n'))
    wws.export_to_hdf5('ww.h5')
    with h5py.File('ww.h5') as f:
        assert f.attrs['filetype'] == b'weight_windows'
        assert list(f.attrs['version']) == [1, 0]
        wws_grp = f['weight_windows']
        assert int(wws_grp.attrs['n_weight_windows']) == len(wws)
        for ww in wws:
            g = wws_grp[f'weight_windows_{ww.id}']
            # 2D shape is required by the C++ tensor::Tensor<double> reader.
            assert g['lower_ww_bounds'].ndim == 2
            assert g['lower_ww_bounds'].shape[0] == ww.num_energy_bins
            # max_lower_bound_ratio is read unconditionally by C++.
            assert 'max_lower_bound_ratio' in g
        m_grp = f['meshes']
        for name in m_grp:
            assert 'id' in m_grp[name].attrs
            assert 'type' in m_grp[name]
