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


def test_get_tally_filter_type(run_in_tmpdir):
    """Test various ways of retrieving tallies from a StatePoint object."""

    mat = openmc.Material()
    mat.add_nuclide("H1", 1.0)
    mat.set_density("g/cm3", 10.0)

    sphere = openmc.Sphere(r=10.0, boundary_type="vacuum")
    cell = openmc.Cell(fill=mat, region=-sphere)
    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.particles = 10
    settings.batches = 2
    settings.run_mode = "fixed source"

    reg_mesh = openmc.RegularMesh().from_domain(cell)
    tally1 = openmc.Tally(tally_id=1)
    mesh_filter = openmc.MeshFilter(reg_mesh)
    tally1.filters = [mesh_filter]
    tally1.scores = ["flux"]

    tally2 = openmc.Tally(tally_id=2, name="heating tally")
    cell_filter = openmc.CellFilter(cell)
    tally2.filters = [cell_filter]
    tally2.scores = ["heating"]

    tallies = openmc.Tallies([tally1, tally2])
    model = openmc.Model(
        geometry=geometry, materials=[mat], settings=settings, tallies=tallies
    )

    sp_filename = model.run()

    sp = openmc.StatePoint(sp_filename)

    tally_found = sp.get_tally(filter_type=openmc.MeshFilter)
    assert tally_found.id == 1

    tally_found = sp.get_tally(filter_type=openmc.CellFilter)
    assert tally_found.id == 2

    tally_found = sp.get_tally(filters=[mesh_filter])
    assert tally_found.id == 1

    tally_found = sp.get_tally(filters=[cell_filter])
    assert tally_found.id == 2

    tally_found = sp.get_tally(scores=["heating"])
    assert tally_found.id == 2

    tally_found = sp.get_tally(name="heating tally")
    assert tally_found.id == 2

    tally_found = sp.get_tally(name=None)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=1)
    assert tally_found.id == 1

    tally_found = sp.get_tally(id=2)
    assert tally_found.id == 2
