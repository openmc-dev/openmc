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
    assert openmc.MGXSLibrary.from_hdf5(path).particle_type == \
        openmc.ParticleType.PHOTON


def test_all_mgxs_types_respects_particle_type():
    geometry = openmc.Geometry([openmc.Cell()])
    photon_library = openmc.mgxs.Library(
        geometry, mgxs_types='all', particle_type='photon')

    assert photon_library.mgxs_types == (
        'total', 'absorption', 'nu-scatter matrix')
