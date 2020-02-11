import openmc
import openmc.stats


def test_export_to_xml(run_in_tmpdir):
    s = openmc.Settings()
    s.run_mode = 'fixed source'
    s.batches = 1000
    s.generations_per_batch = 10
    s.inactive = 100
    s.particles = 1000000
    s.max_lost_particles = 5
    s.rel_max_lost_particles = 1e-4
    s.keff_trigger = {'type': 'std_dev', 'threshold': 0.001}
    s.energy_mode = 'continuous-energy'
    s.max_order = 5
    s.source = openmc.Source(space=openmc.stats.Point())
    s.output = {'summary': True, 'tallies': False, 'path': 'here'}
    s.verbosity = 7
    s.sourcepoint = {'batches': [50, 150, 500, 1000], 'separate': True,
                     'write': True, 'overwrite': True}
    s.statepoint = {'batches': [50, 150, 500, 1000]}
    s.confidence_intervals = True
    s.ptables = True
    s.seed = 17
    s.survival_biasing = True
    s.cutoff = {'weight': 0.25, 'weight_avg': 0.5, 'energy_neutron': 1.0e-5,
                'energy_photon': 1000.0, 'energy_electron': 1.0e-5,
                'energy_positron': 1.0e-5}
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
    s.trace = (10, 1, 20)
    s.track = [1, 1, 1, 2, 1, 1]
    s.ufs_mesh = mesh
    s.resonance_scattering = {'enable': True, 'method': 'rvs',
                              'energy_min': 1.0, 'energy_max': 1000.0,
                              'nuclides': ['U235', 'U238', 'Pu239']}
    s.volume_calculations = openmc.VolumeCalculation(
        domains=[openmc.Cell()], samples=1000, lower_left=(-10., -10., -10.),
        upper_right = (10., 10., 10.))
    s.create_fission_neutrons = True
    s.log_grid_bins = 2000
    s.photon_transport = False
    s.electron_treatment = 'led'
    s.dagmc = False

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
    assert isinstance(s.source[0], openmc.Source)
    assert isinstance(s.source[0].space, openmc.stats.Point)
    assert s.output == {'summary': True, 'tallies': False, 'path': 'here'}
    assert s.verbosity == 7
    assert s.sourcepoint == {'batches': [50, 150, 500, 1000], 'separate': True,
                             'write': True, 'overwrite': True}
    assert s.statepoint == {'batches': [50, 150, 500, 1000]}
    assert s.confidence_intervals
    assert s.ptables
    assert s.seed == 17
    assert s.survival_biasing
    assert s.cutoff == {'weight': 0.25, 'weight_avg': 0.5,
                        'energy_neutron': 1.0e-5, 'energy_photon': 1000.0,
                        'energy_electron': 1.0e-5, 'energy_positron': 1.0e-5}
    assert isinstance(s.entropy_mesh, openmc.RegularMesh)
    assert s.entropy_mesh.lower_left == [-10., -10., -10.]
    assert s.entropy_mesh.upper_right == [10., 10., 10.]
    assert s.entropy_mesh.dimension == [5, 5, 5]
    assert s.trigger_active
    assert s.trigger_max_batches == 10000
    assert s.trigger_batch_interval == 50
    assert not s.no_reduce
    assert s.tabular_legendre == {'enable': True, 'num_points': 50}
    assert s.temperature == {'default': 293.6, 'method': 'interpolation',
                             'multipole': True, 'range': [200., 1000.]}
    assert s.trace == [10, 1, 20]
    assert s.track == [1, 1, 1, 2, 1, 1]
    assert isinstance(s.ufs_mesh, openmc.RegularMesh)
    assert s.ufs_mesh.lower_left == [-10., -10., -10.]
    assert s.ufs_mesh.upper_right == [10., 10., 10.]
    assert s.ufs_mesh.dimension == [5, 5, 5]
    assert s.resonance_scattering == {'enable': True, 'method': 'rvs',
                                      'energy_min': 1.0, 'energy_max': 1000.0,
                                      'nuclides': ['U235', 'U238', 'Pu239']}
    assert s.create_fission_neutrons
    assert s.log_grid_bins == 2000
    assert not s.photon_transport
    assert s.electron_treatment == 'led'
    assert not s.dagmc
