import numpy as np
import openmc
from pytest import approx


def test_filter_mesh_material(run_in_tmpdir):
    # Create four identical materials
    openmc.reset_auto_ids()
    materials = []
    for i in range(4):
        mat = openmc.Material()
        mat.add_nuclide('Fe56', 1.0)
        materials.append(mat)

    # Create a slab model with four cells
    z_values = [-10., -5., 0., 5., 10.]
    planes = [openmc.ZPlane(z) for z in z_values]
    print(planes)
    planes[0].boundary_type = 'vacuum'
    planes[-1].boundary_type = 'vacuum'
    regions = [+left & -right for left, right in zip(planes[:-1], planes[1:])]
    cells = [openmc.Cell(fill=m, region=r) for r, m in zip(regions, materials)]
    model = openmc.Model()
    model.geometry = openmc.Geometry(cells)
    model.settings.particles = 1_000
    model.settings.batches = 5
    model.settings.run_mode = 'fixed source'

    # Create a mesh that does not align with all planar surfaces
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1., -1., -10.)
    mesh.upper_right = (1., 1., 10.)
    mesh.dimension = (1, 1, 5)

    # Determine material volumes in each mesh element and use result to create a
    # MeshMaterialFilter with corresponding bins
    vols = mesh.material_volumes(model)
    mmf = openmc.MeshMaterialFilter.from_volumes(mesh, vols)
    expected_bins = [(0, 1), (1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]
    np.testing.assert_equal(mmf.bins, expected_bins)

    # Create two tallies, one with a mesh filter and one with mesh-material
    mesh_tally = openmc.Tally()
    mesh_tally.filters = [openmc.MeshFilter(mesh)]
    mesh_tally.scores = ['flux']
    mesh_material_tally = openmc.Tally()
    mesh_material_tally.filters = [mmf]
    mesh_material_tally.scores = ['flux']
    model.tallies = [mesh_tally, mesh_material_tally]

    # Run model to get results on the two tallies
    model.run(apply_tally_results=True)

    # The sum of the flux in each mesh-material combination within a single mesh
    # element should be equal to the flux in that mesh element
    mesh_mean = mesh_tally.mean.ravel()
    meshmat_mean = mesh_material_tally.mean.ravel()
    assert mesh_mean[0] == approx(meshmat_mean[0])
    assert mesh_mean[1] == approx(meshmat_mean[1] + meshmat_mean[2])
    assert mesh_mean[2] == approx(meshmat_mean[3] + meshmat_mean[4])
    assert mesh_mean[3] == approx(meshmat_mean[5] + meshmat_mean[6])
    assert mesh_mean[4] == approx(meshmat_mean[7])
    assert mesh_tally.mean.sum() == approx(mesh_material_tally.mean.sum())

    # Make sure get_pandas_dataframe method works
    mesh_material_tally.get_pandas_dataframe()
