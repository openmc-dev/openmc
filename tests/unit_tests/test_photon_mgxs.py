import copy

import numpy as np
import openmc
import pytest


def _add_tally_data(mgxs):
    for tally in mgxs.tallies.values():
        tally._derived = True
        tally._mean = np.ones(tally.shape)
        tally._std_dev = np.zeros(tally.shape)


@pytest.fixture
def photon_library():
    geometry = openmc.Geometry([openmc.Cell(fill=openmc.Material())])
    library = openmc.mgxs.Library(
        geometry, mgxs_types='all', particle_type='photon')
    library.domain_type = 'material'
    library.energy_groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    library.build_library()
    library._sp_filename = 'statepoint.h5'
    for mgxs in library.all_mgxs[library.domains[0].id].values():
        _add_tally_data(mgxs)
    return library


def test_photon_multiplicity(photon_library):
    domain = photon_library.domains[0]
    total = photon_library.get_mgxs(domain, 'total')
    absorption = photon_library.get_mgxs(domain, 'absorption')
    production = photon_library.get_mgxs(domain, 'nu-scatter matrix')
    total.tallies['total']._mean[:, 0, 0] = [0.0, 3.0]
    absorption.tallies['absorption']._mean[:, 0, 0] = [0.0, 1.0]
    production.tallies['primary photon production']._mean[:, 0, 0] = [
        0.0, 0.0, 4.0, 2.0]
    production.tallies['secondary photon production']._mean.fill(0.0)

    xsdata = photon_library.create_mg_library().xsdatas[0]
    assert np.array_equal(
        xsdata.multiplicity_matrix[0], [[3.0, 3.0], [1.0, 1.0]])
    assert np.array_equal(
        xsdata.scatter_matrix[0][:, :, 0], [[2.0, 4.0], [0.0, 0.0]])


@pytest.mark.parametrize('absorption, production', [
    (1.0, 1.0),
    (0.0, 0.0),
])
def test_invalid_photon_multiplicity(photon_library, absorption, production):
    domain = photon_library.domains[0]
    photon_library.get_mgxs(domain, 'absorption').tallies[
        'absorption']._mean.fill(absorption)
    transfer = photon_library.get_mgxs(domain, 'nu-scatter matrix')
    transfer.tallies['primary photon production']._mean.fill(production)
    transfer.tallies['secondary photon production']._mean.fill(0.0)

    with pytest.raises(ValueError, match='Photon MGXS data'):
        photon_library.create_mg_library()


@pytest.mark.parametrize('mgxs_type', [
    'transport',
    'nu-transport',
    'diffusion-coefficient',
    'nu-diffusion-coefficient',
    'consistent scatter matrix',
    'consistent nu-scatter matrix',
])
def test_explicit_neutron_particle_filter(mgxs_type):
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    mgxs = openmc.mgxs.MGXS.get_mgxs(
        mgxs_type, material, 'material', groups,
        particle_type='neutron')
    _add_tally_data(mgxs)
    if 'total' in mgxs.tallies:
        mgxs.tallies['total']._mean.fill(2.0)

    for tally in mgxs.tallies.values():
        particle_filter = tally.find_filter(openmc.ParticleFilter)
        assert particle_filter.bins == ['neutron']

    assert np.all(np.isfinite(mgxs.get_xs()))


def test_photon_production_matrix_tallies():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    production = openmc.mgxs.PhotonTransferMatrixXS(
        material, 'material', groups)

    assert production.particle_type == openmc.ParticleType.PHOTON
    assert production.estimator == ['tracklength', 'analog', 'analog']
    with pytest.raises(ValueError, match='fixed tracklength and analog'):
        production.estimator = 'analog'
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


def test_photon_mgxs_library(photon_library, tmp_path):
    library = photon_library
    material = library.domains[0]
    groups = library.energy_groups
    library.check_library_for_openmc_mgxs()
    assert library.correction is None

    clone = copy.deepcopy(library)
    clone.check_library_for_openmc_mgxs()
    assert clone.estimator is None

    for mgxs_type, mgxs in library.all_mgxs[material.id].items():
        assert mgxs.particle_type == openmc.ParticleType.PHOTON
        if mgxs_type == 'total':
            mgxs.tallies['total']._mean.fill(3.0)

    coarse_groups = openmc.mgxs.EnergyGroups([1.0, 100.0])
    condensed = library.get_condensed_library(coarse_groups)
    xsdata = condensed.create_mg_library().xsdatas[0]
    assert xsdata.scatter_matrix[0].shape == (1, 1, 1)
    production = xsdata.scatter_matrix[0][:, :, 0]
    multiplicity = xsdata.multiplicity_matrix[0]
    physical_scatter = xsdata.total[0] - xsdata.absorption[0]
    assert multiplicity is not None
    assert np.allclose(np.sum(production / multiplicity, axis=1),
                       physical_scatter)

    mg_library = openmc.MGXSLibrary(groups, particle_type='photon')
    path = tmp_path / 'mgxs.h5'
    mg_library.export_to_hdf5(path)
    assert openmc.MGXSLibrary.from_hdf5(path).particle_type == \
        openmc.ParticleType.PHOTON


def test_all_mgxs_types_respects_particle_type():
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    photon_library = openmc.mgxs.Library(geometry, mgxs_types='all')
    photon_library.particle_type = 'photon'
    photon_library.domain_type = 'material'
    photon_library.energy_groups = openmc.mgxs.EnergyGroups([1.0, 10.0])
    photon_library.build_library()

    assert photon_library.mgxs_types == (
        'total', 'absorption', 'nu-scatter matrix')
    assert photon_library.correction is None

    clone = copy.deepcopy(photon_library)
    clone.particle_type = 'neutron'
    assert clone.mgxs_types == (
        openmc.mgxs.MGXS_TYPES + openmc.mgxs.MDGXS_TYPES +
        openmc.mgxs.ARBITRARY_VECTOR_TYPES +
        openmc.mgxs.ARBITRARY_MATRIX_TYPES)


@pytest.mark.parametrize('attribute, value, error', [
    ('correction', 'P0', 'correction'),
    ('scatter_format', 'histogram', 'scatter_format'),
    ('legendre_order', 1, 'legendre_order'),
    ('histogram_bins', 8, 'histogram_bins'),
    ('estimator', 'analog', 'estimator'),
    ('num_polar', 2, 'num_polar'),
])
def test_invalid_photon_library_setting(attribute, value, error):
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    library = openmc.mgxs.Library(
        geometry, mgxs_types='all', particle_type='photon')
    library.domain_type = 'material'
    library.energy_groups = openmc.mgxs.EnergyGroups([1.0, 10.0])
    setattr(library, attribute, value)

    with pytest.raises(ValueError, match=error):
        library.build_library()


@pytest.mark.parametrize('mgxs_types, by_nuclide, error', [
    (['total', 'absorption', 'scatter matrix'], False, 'unsupported'),
    (['total', 'absorption', 'nu-scatter matrix'], True, 'by_nuclide'),
])
def test_invalid_photon_library_scope(mgxs_types, by_nuclide, error):
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    library = openmc.mgxs.Library(
        geometry, by_nuclide=by_nuclide, mgxs_types=mgxs_types,
        particle_type='photon')
    library.domain_type = 'material'
    library.energy_groups = openmc.mgxs.EnergyGroups([1.0, 10.0])

    with pytest.raises(ValueError, match=error):
        library.build_library()


def test_photon_library_requires_complete_export_data():
    geometry = openmc.Geometry([openmc.Cell(fill=openmc.Material())])
    library = openmc.mgxs.Library(
        geometry, mgxs_types=['total'], particle_type='photon')
    library.domain_type = 'material'
    library.energy_groups = openmc.mgxs.EnergyGroups([1.0, 10.0])
    library.build_library()

    with pytest.raises(ValueError, match='missing required'):
        library.check_library_for_openmc_mgxs()


def test_explicit_neutron_mdgxs_particle_filter():
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0])
    library = openmc.mgxs.Library(
        geometry, mgxs_types=['delayed-nu-fission'],
        particle_type='neutron')
    library.domain_type = 'material'
    library.energy_groups = groups
    library.build_library()
    mdgxs = library.get_mgxs(material, 'delayed-nu-fission')

    assert mdgxs.particle_type == openmc.ParticleType.NEUTRON
    for tally in mdgxs.tallies.values():
        assert tally.find_filter(openmc.ParticleFilter).bins == ['neutron']


def test_explicit_neutron_library_matches_unfiltered_export():
    material = openmc.Material()
    geometry = openmc.Geometry([openmc.Cell(fill=material)])
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    xsdatas = []

    for particle_type in (None, 'neutron'):
        library = openmc.mgxs.Library(
            geometry,
            mgxs_types=['transport', 'absorption', 'scatter matrix'],
            particle_type=particle_type)
        library.domain_type = 'material'
        library.energy_groups = groups
        library.build_library()
        for mgxs in library.all_mgxs[material.id].values():
            _add_tally_data(mgxs)
            if 'total' in mgxs.tallies:
                mgxs.tallies['total']._mean.fill(2.0)
            for tally in mgxs.tallies.values():
                if particle_type is None:
                    assert not tally.contains_filter(openmc.ParticleFilter)
                else:
                    assert tally.find_filter(openmc.ParticleFilter).bins == [
                        'neutron']
        library._sp_filename = 'statepoint.h5'
        with pytest.warns((UserWarning, RuntimeWarning)):
            xsdatas.append(library.create_mg_library().xsdatas[0])

    for attribute in ('total', 'absorption', 'scatter_matrix'):
        assert np.array_equal(
            getattr(xsdatas[0], attribute)[0],
            getattr(xsdatas[1], attribute)[0])


def test_photon_production_matrix_group_transformations():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    production = openmc.mgxs.PhotonTransferMatrixXS(
        material, 'material', groups)
    for tally in production.tallies.values():
        tally._derived = True
        tally._mean = np.zeros(tally.shape)
        tally._std_dev = np.zeros(tally.shape)
    production.tallies['flux']._mean.fill(1.0)
    production.tallies[
        'primary photon production']._mean[:, 0, 0] = np.arange(1.0, 10.0)

    energy = production.tallies['flux'].find_filter(openmc.EnergyFilter)
    assert all(tally.find_filter(openmc.EnergyFilter) is energy
               for tally in production.tallies.values())

    coarse_groups = openmc.mgxs.EnergyGroups([1.0, 100.0, 1000.0])
    condensed = production.get_condensed_xs(coarse_groups)
    assert np.array_equal(condensed.get_xs(), [[9.0, 15.0], [4.5, 6.0]])
    assert condensed.tallies['secondary photon production'].contains_filter(
        openmc.EnergyoutFilter)
    assert np.array_equal(energy.values, groups.group_edges)
    assert np.array_equal(
        production.get_xs(), [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0],
                              [3.0, 2.0, 1.0]])

    sliced = production.get_slice(in_groups=[1, 2], out_groups=[1, 2])
    assert np.array_equal(sliced.get_xs(), [[9.0, 8.0], [6.0, 5.0]])
    assert sliced.tallies['secondary photon production'].contains_filter(
        openmc.EnergyoutFilter)
    assert production.tallies[
        'secondary photon production'].contains_filter(
            openmc.ParticleProductionFilter)

    resliced = sliced.get_slice(in_groups=[1], out_groups=[1])
    assert np.array_equal(resliced.get_xs(), [[9.0]])

    in_sliced = production.get_slice(in_groups=[1])
    assert np.array_equal(in_sliced.get_xs(), [[9.0, 8.0, 7.0]])

    out_sliced = production.get_slice(out_groups=[2, 3])
    assert np.array_equal(
        out_sliced.get_xs(), [[8.0, 7.0], [5.0, 4.0], [2.0, 1.0]])

    assert np.array_equal(
        production.get_xs(in_groups=[1], out_groups=[2, 3]),
        [[8.0, 7.0]])


def test_vector_mgxs_slice_order():
    material = openmc.Material()
    groups = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0, 1000.0])
    total = openmc.mgxs.TotalXS(material, 'material', groups)
    _add_tally_data(total)
    total.tallies['total']._mean[:, 0, 0] = [1.0, 2.0, 3.0]

    sliced = total.get_slice(groups=[1, 2])
    assert np.array_equal(sliced.get_xs(), [3.0, 2.0])
