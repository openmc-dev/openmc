import copy

import numpy as np
import openmc


def _add_tally_data(mgxs):
    for tally in mgxs.tallies.values():
        tally._derived = True
        tally._mean = np.ones(tally.shape)
        tally._std_dev = np.zeros(tally.shape)


def test_photon_production_matrix_tallies():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    production = openmc.mgxs.PhotonTransferMatrixXS(
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
    library.build_library()
    library.check_library_for_openmc_mgxs()
    assert library.correction is None

    clone = copy.deepcopy(library)
    clone.check_library_for_openmc_mgxs()
    assert clone.estimator is None

    for mgxs in library.all_mgxs[material.id].values():
        assert mgxs.particle_type == openmc.ParticleType.PHOTON
        _add_tally_data(mgxs)

    library._sp_filename = 'statepoint.h5'
    coarse_groups = openmc.mgxs.EnergyGroups([1.0, 100.0])
    condensed = library.get_condensed_library(coarse_groups)
    xsdata = condensed.create_mg_library().xsdatas[0]
    assert xsdata.scatter_matrix[0].shape == (1, 1, 1)

    mg_library = openmc.MGXSLibrary(groups, particle_type='photon')
    path = tmp_path / 'mgxs.h5'
    mg_library.export_to_hdf5(path)
    assert openmc.MGXSLibrary.from_hdf5(path).particle_type == \
        openmc.ParticleType.PHOTON


def test_all_mgxs_types_respects_particle_type():
    geometry = openmc.Geometry([openmc.Cell()])
    photon_library = openmc.mgxs.Library(
        geometry, mgxs_types='all', particle_type='photon')

    assert photon_library.mgxs_types == (
        'total', 'absorption', 'nu-scatter matrix')
    assert photon_library.correction is None


def test_photon_production_matrix_group_transformations():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    production = openmc.mgxs.PhotonTransferMatrixXS(
        material, 'material', groups)
    _add_tally_data(production)

    coarse_groups = openmc.mgxs.EnergyGroups([1.0, 100.0, 1000.0])
    condensed = production.get_condensed_xs(coarse_groups)
    assert condensed.get_xs().shape == (2, 2)
    assert condensed.tallies['secondary photon production'].contains_filter(
        openmc.EnergyoutFilter)

    sliced = production.get_slice(in_groups=[1, 2], out_groups=[1, 2])
    assert sliced.get_xs().shape == (2, 2)
    assert sliced.tallies['secondary photon production'].contains_filter(
        openmc.EnergyoutFilter)
    assert production.tallies[
        'secondary photon production'].contains_filter(
            openmc.ParticleProductionFilter)

    resliced = sliced.get_slice(in_groups=[1], out_groups=[1])
    assert resliced.get_xs().shape == (1, 1)
