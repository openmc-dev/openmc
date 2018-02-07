import openmc
import openmc.stats


def test_export_to_xml(run_in_tmpdir):
    s = openmc.Settings()
    s.run_mode = 'fixed source'
    s.batches = 1000
    s.generations_per_batch = 10
    s.inactive = 100
    s.particles = 1000000
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
    s.cross_sections = '/path/to/cross_sections.xml'
    s.multipole_library = '/path/to/wmp/'
    s.ptables = True
    s.run_cmfd = False
    s.seed = 17
    s.survival_biasing = True
    s.cutoff = {'weight': 0.25, 'weight_avg': 0.5, 'energy': 1.0e-5}
    mesh = openmc.Mesh()
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
    s.threads = 8
    s.trace = (10, 1, 20)
    s.track = [1, 1, 1, 2, 1, 1]
    s.ufs_mesh = mesh
    s.resonance_scattering = {'enable': True, 'method': 'ares',
                              'energy_min': 1.0, 'energy_max': 1000.0,
                              'nuclides': ['U235', 'U238', 'Pu239']}
    s.volume_calculations = openmc.VolumeCalculation(
        domains=[openmc.Cell()], samples=1000, lower_left=(-10., -10., -10.),
        upper_right = (10., 10., 10.))
    s.create_fission_neutrons = True
    s.log_grid_bins = 2000

    # Make sure exporting XML works
    s.export_to_xml()
