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
