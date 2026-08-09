import h5py
import numpy as np
import pytest

import openmc


def test_photon_production_matrix_tallies():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    production = openmc.mgxs.PhotonProductionMatrixXS(
        material, 'material', groups)

    assert production.particle_type == openmc.ParticleType.PHOTON
    assert production.estimator == ['tracklength', 'analog', 'analog']
    assert production.tally_keys == [
        'flux', 'primary photon production', 'secondary photon production']

    tallies = production.tallies
    assert tallies['primary photon production'].scores == ['scatter']
    assert tallies['secondary photon production'].scores == ['events']
    assert tallies['secondary photon production'].contains_filter(
        openmc.ParticleProductionFilter)
    for tally in tallies.values():
        particle_filter = tally.find_filter(openmc.ParticleFilter)
        assert particle_filter.bins == ['photon']


def test_photon_production_matrix_constraints():
    with pytest.raises(ValueError, match='by nuclide'):
        openmc.mgxs.MGXS.get_mgxs(
            'nu-scatter matrix', by_nuclide=True, particle_type='photon')
    with pytest.raises(ValueError, match='isotropic'):
        openmc.mgxs.MGXS.get_mgxs(
            'nu-scatter matrix', num_polar=2, particle_type='photon')

    neutron_production = openmc.mgxs.MGXS.get_mgxs(
        'nu-scatter matrix', particle_type='neutron')
    assert isinstance(neutron_production, openmc.mgxs.ScatterMatrixXS)


def test_photon_production_matrix_combines_primary_and_secondary():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    production = openmc.mgxs.PhotonProductionMatrixXS(
        material, 'material', groups)

    for tally in production.tallies.values():
        values = np.arange(1, tally.num_filter_bins + 1, dtype=float)
        values.shape = (tally.num_filter_bins, 1, 1)
        tally._sum = 2.0 * values
        tally._sum_sq = 2.0 * values**2
        tally._num_realizations = 2
        tally._sp_filename = 'statepoint.h5'
        tally._results_read = True

    np.testing.assert_allclose(
        production.rxn_rate_tally.mean.ravel(), [2.0, 4.0, 6.0, 8.0])


def test_set_photon_nu_scatter_mgxs(monkeypatch):
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    production = openmc.mgxs.PhotonProductionMatrixXS(
        material, 'material', groups)
    values = np.array([[0.1, 0.2], [0.3, 0.4]])
    monkeypatch.setattr(production, 'get_xs', lambda **kwargs: values)

    xsdata = openmc.XSdata('photon', groups)
    xsdata.order = 0
    absorption = np.array([0.5, 0.6])
    xsdata.set_absorption(absorption)
    xsdata.set_scatter_matrix_mgxs(production)

    assert xsdata.scatter_format == 'legendre'
    assert xsdata.order == 0
    np.testing.assert_allclose(xsdata._scatter_matrix[0][:, :, 0], values)
    np.testing.assert_allclose(xsdata._absorption[0], absorption)
    assert xsdata._multiplicity_matrix[0] is None


def test_photon_mgxs_library(tmp_path):
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    library = openmc.mgxs.Library(
        geometry, mgxs_types=[
            'total', 'absorption', 'nu-scatter matrix'],
        particle_type='photon')
    library.domain_type = 'material'
    library.energy_groups = groups
    library.correction = None
    library.build_library()
    library.check_library_for_openmc_mgxs()

    for mgxs in library.all_mgxs[material.id].values():
        assert mgxs.particle_type == openmc.ParticleType.PHOTON

    mg_library = openmc.MGXSLibrary(groups, particle_type='photon')
    path = tmp_path / 'mgxs.h5'
    mg_library.export_to_hdf5(path)
    with h5py.File(path) as h5file:
        assert h5file.attrs['particle_type'] == b'photon'
    assert openmc.MGXSLibrary.from_hdf5(path).particle_type == \
        openmc.ParticleType.PHOTON


def test_all_mgxs_types_respects_particle_type():
    geometry = openmc.Geometry([openmc.Cell()])
    neutron_library = openmc.mgxs.Library(geometry, mgxs_types='all')
    photon_library = openmc.mgxs.Library(
        geometry, mgxs_types='all', particle_type='photon')

    assert photon_library.mgxs_types == (
        'total', 'absorption', 'nu-scatter matrix')
