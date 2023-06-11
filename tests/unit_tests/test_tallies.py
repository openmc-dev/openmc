import numpy as np
import matplotlib

import openmc


def test_xml_roundtrip(run_in_tmpdir):
    # Create a tally with all possible gizmos
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10., -10., -10.)
    mesh.upper_right = (10., 10., 10.,)
    mesh.dimension = (5, 5, 5)
    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.nuclides = ['U235', 'I135', 'Li6']
    tally.scores = ['total', 'fission', 'heating']
    tally.derivative = openmc.TallyDerivative(
        variable='nuclide_density', material=1, nuclide='Li6'
    )
    tally.triggers = [openmc.Trigger('rel_err', 0.025)]
    tally.triggers[0].scores = ['total', 'fission']
    tallies = openmc.Tallies([tally])

    # Roundtrip through XML and make sure we get what we started with
    tallies.export_to_xml()
    new_tallies = openmc.Tallies.from_xml()
    assert len(new_tallies) == 1
    new_tally = new_tallies[0]
    assert new_tally.id == tally.id
    assert len(new_tally.filters) == 1
    assert isinstance(new_tally.filters[0], openmc.MeshFilter)
    assert np.allclose(new_tally.filters[0].mesh.lower_left, mesh.lower_left)
    assert new_tally.nuclides == tally.nuclides
    assert new_tally.scores == tally.scores
    assert new_tally.derivative.variable == tally.derivative.variable
    assert new_tally.derivative.material == tally.derivative.material
    assert new_tally.derivative.nuclide == tally.derivative.nuclide
    assert len(new_tally.triggers) == 1
    assert new_tally.triggers[0].trigger_type == tally.triggers[0].trigger_type
    assert new_tally.triggers[0].threshold == tally.triggers[0].threshold
    assert new_tally.triggers[0].scores == tally.triggers[0].scores


def test_get_values_slice_from_mesh():

    material = openmc.Material()
    material.set_density('g/cm3', 7)
    material.add_nuclide('Fe56', 1)
    materials = openmc.Materials([material])

    surf1 = openmc.Sphere(r=100, boundary_type='vacuum')
    cell1 = openmc.Cell(region=-surf1)
    my_geometry = openmc.Geometry([cell1])

    settings = openmc.Settings()
    settings.batches = 2
    settings.particles = 1000
    settings.run_mode = 'fixed source'
    settings.source = openmc.Source()

    mesh = openmc.RegularMesh().from_domain(my_geometry, dimension=[15, 20, 25])
    mesh_filter = openmc.MeshFilter(mesh)

    mesh_tally_1 = openmc.Tally(name='flux_on_mesh')
    mesh_tally_1.filters = [mesh_filter]
    mesh_tally_1.scores = ['flux']
    tallies = openmc.Tallies([mesh_tally_1])

    model = openmc.Model(my_geometry, materials, settings, tallies)
    sp_filename = model.run()

    sp = openmc.StatePoint(sp_filename)

    tally = sp.get_tally(name='flux_on_mesh')

    for basis in ['xz', 'yz', 'xy']:
        slice_data = tally.get_values_slice_from_mesh(
            basis=basis,
            slice_index=10,
            value='std_dev'
        )
        returned_axis = mesh.plot_tally_values_slice(
            dataset=slice_data,
            basis=basis,
        )
        assert isinstance(returned_axis, matplotlib.axes.Axes)

        slice_data = tally.get_values_slice_from_mesh_where(basis=basis, slice_value=10)
        returned_axis = mesh.plot_tally_values_slice(
            dataset=slice_data,
            basis=basis,
            colorbar_label='Flux'
        )
        assert isinstance(returned_axis, matplotlib.axes.Axes)
        returned_axis.set_title('checking plot title can be set')
