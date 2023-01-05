import openmc


def test_slice_shape_and_normalization():
    """Simulates a simple mesh tally and checks slices of the result are
    shaped and normalization corrected"""
    mat_1 = openmc.Material()
    mat_1.add_element("Li", 1.0)
    mat_1.set_density("g/cm3", 0.01)
    my_materials = openmc.Materials([mat_1])

    central_sphere_surf = openmc.Sphere(r=200, boundary_type="vacuum")
    central_sphere_region = -central_sphere_surf
    central_sphere_cell = openmc.Cell(region=central_sphere_region)
    universe = openmc.Universe(cells=[central_sphere_cell])
    my_geometry = openmc.Geometry(universe)

    my_settings = openmc.Settings()
    my_settings.batches = 1
    my_settings.particles = 100
    my_settings.run_mode = "fixed source"
    my_source = openmc.Source()
    my_settings.source = my_source

    mesh = openmc.RegularMesh().from_domain(domain=my_geometry, dimension=(5, 7, 3))
    mesh_filter = openmc.MeshFilter(mesh)

    mesh_tally = openmc.Tally(name="test_tally")
    mesh_tally.filters = [mesh_filter]
    mesh_tally.scores = ["flux"]
    my_tallies = openmc.Tallies([mesh_tally])

    model = openmc.model.Model(my_geometry, my_materials, my_settings, my_tallies)
    sp_filename = model.run()

    statepoint = openmc.StatePoint(sp_filename)
    tally = statepoint.get_tally(name="test_tally")

    view_dirs = ["x", "-x", "y", "-y", "z", "-z"]
    shapes = [(3, 7), (3, 7), (3, 5), (3, 5), (7, 5), (7, 5)]
    for view_dir, shape in zip(view_dirs, shapes):

        norm_vals = mesh.slice_of_data(
            dataset=tally.mean,
            view_direction=view_dir,
            slice_index=1,
            volume_normalization=True,
        )
        assert norm_vals.shape == shape
        vals = mesh.slice_of_data(
            dataset=tally.mean,
            view_direction=view_dir,
            slice_index=1,
            volume_normalization=False,
        )
        assert vals.shape == shape
        for val, norm_val in zip(vals.flatten(), norm_vals.flatten()):
            # all mesh volumes in regular mesh are the same
            assert val / mesh.volumes[0][0][0] == norm_val
