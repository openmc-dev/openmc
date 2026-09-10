from pathlib import Path

import h5py
import pytest

import openmc
import openmc.lib
import openmc.stats


def test_export_to_xml(run_in_tmpdir):

    tmp_properties_file = 'properties_test.h5'

    s = openmc.Settings(run_mode='fixed source', batches=1000, seed=17)
    s.generations_per_batch = 10
    s.inactive = 100
    s.particles = 1000000
    s.max_lost_particles = 5
    s.rel_max_lost_particles = 1e-4
    s.keff_trigger = {'type': 'std_dev', 'threshold': 0.001}
    s.energy_mode = 'continuous-energy'
    s.max_order = 5
    s.max_tracks = 1234
    s.source = openmc.IndependentSource(space=openmc.stats.Point())
    s.output = {'summary': True, 'tallies': False, 'path': 'here'}
    s.verbosity = 7
    s.sourcepoint = {'batches': [50, 150, 500, 1000], 'separate': True,
                     'write': True, 'overwrite': True, 'mcpl': True}
    s.statepoint = {'batches': [50, 150, 500, 1000]}
    s.surf_source_read = {'path': 'surface_source_1.h5'}
    s.surf_source_write = {'surface_ids': [2], 'max_particles': 200}
    s.surface_grazing_ratio = 0.7
    s.surface_grazing_cutoff = 0.1
    s.confidence_intervals = True
    s.ptables = True
    s.plot_seed = 100
    s.survival_biasing = True
    s.cutoff = {'weight': 0.25, 'weight_avg': 0.5, 'energy_neutron': 1.0e-5,
                'survival_normalization': True,
                'energy_photon': 1000.0, 'energy_electron': 1.0e-5,
                'energy_positron': 1.0e-5, 'time_neutron': 1.0e-5,
                'time_photon': 1.0e-5, 'time_electron': 1.0e-5,
                'time_positron': 1.0e-5}
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10., -10., -10.)
    mesh.upper_right = (10., 10., 10.)
    mesh.dimension = (5, 5, 5)
    s.entropy_mesh = mesh
    s.trigger_active = True
    s.trigger_max_batches = 10000
    s.trigger_batch_interval = 50
    s.no_reduce = False
    s.tabular_legendre = {'enable': True, 'num_points': 50}
    s.temperature = {'default': 293.6, 'method': 'interpolation',
                     'multipole': True, 'range': (200., 1000.)}
    s.properties_file = tmp_properties_file
    s.trace = (10, 1, 20)
    s.track = [(1, 1, 1), (2, 1, 1)]
    s.ufs_mesh = mesh
    s.resonance_scattering = {'enable': True, 'method': 'rvs',
                              'energy_min': 1.0, 'energy_max': 1000.0,
                              'nuclides': ['U235', 'U238', 'Pu239']}
    s.volume_calculations = openmc.VolumeCalculation(
        domains=[openmc.Cell()], samples=1000, lower_left=(-10., -10., -10.),
        upper_right = (10., 10., 10.))
    s.create_fission_neutrons = True
    s.create_delayed_neutrons = False
    s.log_grid_bins = 2000
    s.photon_transport = False
    s.electron_treatment = 'led'
    s.atomic_relaxation = False
    s.write_initial_source = True
    s.weight_window_checkpoints = {'surface': True, 'collision': False}
    source_region_mesh = openmc.RegularMesh()
    source_region_mesh.dimension = [2, 2, 2]
    source_region_mesh.lower_left = [-2, -2, -2]
    source_region_mesh.upper_right = [2, 2, 2]
    root_universe = openmc.Universe()
    s.random_ray = {
        'distance_inactive': 10.0,
        'distance_active': 100.0,
        'ray_source': openmc.IndependentSource(
            space=openmc.stats.Box((-1., -1., -1.), (1., 1., 1.))
        ),
        'source_region_meshes': [(source_region_mesh, [root_universe])],
        'volume_estimator': 'hybrid',
        'source_shape': 'linear',
        'volume_normalized_flux_tallies': True,
        'adjoint': False,
        'sample_method': 'halton'
    }
    s.max_particle_events = 100
    s.max_secondaries = 1_000_000
    s.source_rejection_fraction = 0.01
    s.free_gas_threshold = 800.0
    s.delta_tracking = True

    # Make sure exporting XML works
    s.export_to_xml()


    # Generate settings from XML
    s = openmc.Settings.from_xml()
    assert s.run_mode == 'fixed source'
    assert s.batches == 1000
    assert s.generations_per_batch == 10
    assert s.inactive == 100
    assert s.particles == 1000000
    assert s.max_lost_particles == 5
    assert s.rel_max_lost_particles == 1e-4
    assert s.keff_trigger == {'type': 'std_dev', 'threshold': 0.001}
    assert s.energy_mode == 'continuous-energy'
    assert s.max_order == 5
    assert s.max_tracks == 1234
    assert isinstance(s.source[0], openmc.IndependentSource)
    assert isinstance(s.source[0].space, openmc.stats.Point)
    assert s.output == {'summary': True, 'tallies': False, 'path': 'here'}
    assert s.verbosity == 7
    assert s.sourcepoint == {'batches': [50, 150, 500, 1000], 'separate': True,
                             'write': True, 'overwrite': True, 'mcpl': True}
    assert s.statepoint == {'batches': [50, 150, 500, 1000]}
    assert s.surf_source_read['path'].name == 'surface_source_1.h5'
    assert s.surf_source_write == {'surface_ids': [2], 'max_particles': 200}
    assert s.surface_grazing_ratio == 0.7
    assert s.surface_grazing_cutoff == 0.1
    assert s.confidence_intervals
    assert s.ptables
    assert s.plot_seed == 100
    assert s.seed == 17
    assert s.survival_biasing
    assert s.cutoff == {'weight': 0.25, 'weight_avg': 0.5,
                        'survival_normalization': True,
                        'energy_neutron': 1.0e-5, 'energy_photon': 1000.0,
                        'energy_electron': 1.0e-5, 'energy_positron': 1.0e-5,
                        'time_neutron': 1.0e-5, 'time_photon': 1.0e-5,
                        'time_electron': 1.0e-5, 'time_positron': 1.0e-5}
    assert isinstance(s.entropy_mesh, openmc.RegularMesh)
    assert s.entropy_mesh.lower_left == [-10., -10., -10.]
    assert s.entropy_mesh.upper_right == [10., 10., 10.]
    assert s.entropy_mesh.dimension == (5, 5, 5)
    assert s.trigger_active
    assert s.trigger_max_batches == 10000
    assert s.trigger_batch_interval == 50
    assert not s.no_reduce
    assert s.tabular_legendre == {'enable': True, 'num_points': 50}
    assert s.temperature == {'default': 293.6, 'method': 'interpolation',
                             'multipole': True, 'range': [200., 1000.]}
    assert s.properties_file == Path(tmp_properties_file)
    assert s.trace == [10, 1, 20]
    assert s.track == [(1, 1, 1), (2, 1, 1)]
    assert isinstance(s.ufs_mesh, openmc.RegularMesh)
    assert s.ufs_mesh.lower_left == [-10., -10., -10.]
    assert s.ufs_mesh.upper_right == [10., 10., 10.]
    assert s.ufs_mesh.dimension == (5, 5, 5)
    assert s.resonance_scattering == {'enable': True, 'method': 'rvs',
                                      'energy_min': 1.0, 'energy_max': 1000.0,
                                      'nuclides': ['U235', 'U238', 'Pu239']}
    assert s.create_fission_neutrons
    assert not s.create_delayed_neutrons
    assert s.log_grid_bins == 2000
    assert not s.photon_transport
    assert s.electron_treatment == 'led'
    assert not s.atomic_relaxation
    assert s.write_initial_source
    assert len(s.volume_calculations) == 1
    vol = s.volume_calculations[0]
    assert vol.domain_type == 'cell'
    assert len(vol.ids) == 1
    assert vol.samples == 1000
    assert vol.lower_left == (-10., -10., -10.)
    assert vol.upper_right == (10., 10., 10.)
    assert s.weight_window_checkpoints == {'surface': True, 'collision': False}
    assert s.max_particle_events == 100
    assert s.random_ray['distance_inactive'] == 10.0
    assert s.random_ray['distance_active'] == 100.0
    assert s.random_ray['ray_source'].space.lower_left == [-1., -1., -1.]
    assert s.random_ray['ray_source'].space.upper_right == [1., 1., 1.]
    assert 'source_region_meshes' in s.random_ray
    assert len(s.random_ray['source_region_meshes']) == 1
    mesh_and_domains = s.random_ray['source_region_meshes'][0]
    recovered_mesh = mesh_and_domains[0]
    assert recovered_mesh.dimension == (2, 2, 2)
    assert recovered_mesh.lower_left == [-2., -2., -2.]
    assert recovered_mesh.upper_right == [2., 2., 2.]
    assert s.random_ray['volume_estimator'] == 'hybrid'
    assert s.random_ray['source_shape'] == 'linear'
    assert s.random_ray['volume_normalized_flux_tallies']
    assert not s.random_ray['adjoint']
    assert s.random_ray['sample_method'] == 'halton'
    assert s.max_secondaries == 1_000_000
    assert s.source_rejection_fraction == 0.01
    assert s.free_gas_threshold == 800.0
    assert s.delta_tracking == True


def test_properties_file_load(tmp_path, mpi_intracomm):
    model = openmc.examples.pwr_assembly()

    # Session 1: export a structurally valid properties file via the C++ API,
    # then collect the cell/material structure so we can patch it with h5py.
    cell_instances = {}   # {cell_id: n_instances} — material cells only
    mat_densities = {}    # {mat_id: original atom/b-cm density}

    props_path = tmp_path / 'properties.h5'
    with openmc.lib.TemporarySession(model, intracomm=mpi_intracomm):
        openmc.lib.export_properties(str(props_path))
        for cell_id, cell in openmc.lib.cells.items():
            try:
                cell.fill  # raises NotImplementedError for non-material cells
                cell_instances[cell_id] = cell.num_instances
            except NotImplementedError:
                pass
        for mat_id, mat in openmc.lib.materials.items():
            mat_densities[mat_id] = mat.get_density('atom/b-cm')

    assert any(n > 1 for n in cell_instances.values())

    # Patch the exported properties file overwriting temperatures
    # with per-instance values and scale material atom densities.
    density_factor = 0.75
    with h5py.File(props_path, 'r+') as f:
        cells_grp = f['geometry/cells']
        for cell_id, n in cell_instances.items():
            cell_grp = cells_grp[f'cell {cell_id}']
            del cell_grp['temperature']
            cell_grp.create_dataset(
                'temperature', data=[500.0 + 5.0 * i for i in range(n)]
            )

        for mat_id, orig_density in mat_densities.items():
            f['materials'][f'material {mat_id}'].attrs['atom_density'] = \
                orig_density * density_factor

    # now apply the newly patched properties file using the settings
    # and load the model again, checking that the new temperature and
    # density values match those in the new file
    model.settings.properties_file = props_path

    with openmc.lib.TemporarySession(model, intracomm=mpi_intracomm):
        for cell_id, n in cell_instances.items():
            cell = openmc.lib.cells[cell_id]
            for i in range(n):
                assert cell.get_temperature(i) == pytest.approx(500.0 + 5.0 * i)

        for mat_id, orig_density in mat_densities.items():
            mat = openmc.lib.materials[mat_id]
            assert mat.get_density('atom/b-cm') == pytest.approx(
                orig_density * density_factor, rel=1e-5
            )
