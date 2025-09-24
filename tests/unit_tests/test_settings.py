import openmc
import openmc.stats
import xml.etree.ElementTree as ET


def test_export_to_xml(run_in_tmpdir):
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
    s.write_initial_source = True
    s.weight_window_checkpoints = {'surface': True, 'collision': False}
    s.random_ray = {
        'distance_inactive': 10.0,
        'distance_active': 100.0,
        'ray_source': openmc.IndependentSource(
            space=openmc.stats.Box((-1., -1., -1.), (1., 1., 1.))
        )
    }
    s.max_particle_events = 100
    s.max_secondaries = 1_000_000
    s.source_rejection_fraction = 0.01

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
    assert s.max_secondaries == 1_000_000
    assert s.source_rejection_fraction == 0.01


def test_random_ray_source_region_meshes_export(run_in_tmpdir):
    """Test that random ray source region meshes are properly exported to XML
    
    This is a regression test for the bug where meshes were not included in
    settings.xml when using export_to_xml() but worked correctly with 
    export_to_model_xml().
    """
    # Create minimal settings
    s = openmc.Settings()
    s.particles = 1000
    s.batches = 10
    s.inactive = 5
    s.run_mode = 'eigenvalue'
    
    # Create a mesh for random ray source region meshes
    mesh = openmc.RegularMesh()
    mesh.dimension = [2, 2, 2]
    mesh.lower_left = [-1, -1, -1]
    mesh.upper_right = [1, 1, 1]
    
    # Create a simple universe 
    root_universe = openmc.Universe()
    
    # Set up random ray with source region meshes
    s.random_ray = {
        'source_region_meshes': [(mesh, [root_universe])]
    }
    
    # Test 1: Export to XML using export_to_xml (the problematic case)
    s.export_to_xml()
    
    # Parse the settings XML and check for mesh elements
    tree = ET.parse('settings.xml')
    root = tree.getroot()
    mesh_elements = root.findall('.//mesh')
    
    # There should be 2 mesh elements:
    # 1. The mesh reference inside source_region_meshes 
    # 2. The actual mesh definition at the top level
    assert len(mesh_elements) == 2, f"Expected 2 mesh elements, found {len(mesh_elements)}"
    
    # Check that we have the mesh definition with all required data
    mesh_defs = [elem for elem in mesh_elements if elem.find('dimension') is not None]
    assert len(mesh_defs) == 1, f"Expected 1 mesh definition, found {len(mesh_defs)}"
    
    mesh_def = mesh_defs[0]
    assert mesh_def.get('id') == str(mesh.id)
    assert mesh_def.find('dimension').text == '2 2 2'
    assert mesh_def.find('lower_left').text == '-1 -1 -1'
    assert mesh_def.find('upper_right').text == '1 1 1'
    
    # Test 2: Also verify that to_xml_element works with mesh_memo=None
    xml_element = s.to_xml_element(mesh_memo=None)
    mesh_elements_direct = xml_element.findall('.//mesh')
    assert len(mesh_elements_direct) == 2, f"Expected 2 mesh elements with mesh_memo=None, found {len(mesh_elements_direct)}"
