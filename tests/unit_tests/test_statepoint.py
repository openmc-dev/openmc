import openmc
import pytest

@pytest.mark.parametrize("solver", ["monte carlo", "random ray"])
def test_statepoint(run_in_tmpdir, solver):
    if solver == 'monte carlo':
        model = openmc.examples.pwr_pin_cell()
    else:
        model = openmc.examples.random_ray_lattice()
    sp_file = model.run()
    with openmc.StatePoint(sp_file) as sp:
        current_batch = sp.current_batch
        energy_mode = sp.energy_mode
        n_energy_groups = sp.n_energy_groups
        n_batches = sp.n_batches
        n_inactive = sp.n_inactive
        n_particles = sp.n_particles
        random_ray_dict = sp.random_ray
        run_mode = sp.run_mode
        solver_type = sp.solver_type

    assert current_batch == model.settings.batches
    assert n_batches == model.settings.batches
    assert n_inactive == model.settings.inactive
    assert n_particles == model.settings.particles
    assert run_mode == model.settings.run_mode
    assert solver_type == solver
    if solver_type == 'random ray':
        assert energy_mode == 'multi-group'
        assert n_energy_groups == 7
        assert random_ray_dict['adjoint_mode'] is False
        assert random_ray_dict['avg_miss_rate'] == 0.0
        assert random_ray_dict['distance_active'] == 100.0
        assert random_ray_dict['distance_inactive'] == 20.0
        assert random_ray_dict['sample_method'] == 'prng'
        assert random_ray_dict['source_shape'] == 'flat'
        assert random_ray_dict['n_external_source_regions'] == 0
        assert random_ray_dict['n_geometric_intersections'] == 882796
        assert random_ray_dict['n_integrations'] == 6179572
        assert random_ray_dict['n_source_regions'] == 244
        assert random_ray_dict['volume_estimator'] == 'hybrid'
        assert random_ray_dict['volume_normalized_flux_tallies'] is True
    else:
        assert energy_mode == 'continuous-energy'
        assert n_energy_groups is None
        assert random_ray_dict is None
