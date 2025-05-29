from itertools import product

import numpy as np
from pathlib import Path
import pytest

vtk = pytest.importorskip("vtk")
from vtk.util import numpy_support as nps

import openmc

@pytest.fixture
def model():
    openmc.reset_auto_ids()

    surf1 = openmc.Sphere(r=10, boundary_type='vacuum')
    surf2 = openmc.XPlane(x0=-0.001, boundary_type='vacuum')

    cell = openmc.Cell(region=-surf1 & -surf2)

    geometry = openmc.Geometry([cell])

    settings = openmc.Settings()
    settings.batches = 2
    settings.particles = 100
    settings.run_mode = 'fixed source'

    source = openmc.IndependentSource()
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([1.0e6], [1.0])
    source.space = openmc.stats.Point((-0.01, -0.01, -0.01))

    settings.source = source

    model = openmc.Model(geometry=geometry, settings=settings)

    return model


regular_mesh = openmc.RegularMesh()
regular_mesh.lower_left = (-10, -10, -10)
regular_mesh.upper_right = (10, 10, 10)
regular_mesh.dimension = [30, 20, 10]

rectilinear_mesh = openmc.RectilinearMesh()
rectilinear_mesh.x_grid = np.linspace(-10, 10, 6)
rectilinear_mesh.y_grid = np.logspace(0, 1, 7)
rectilinear_mesh.y_grid = \
    np.concatenate((-rectilinear_mesh.y_grid[::-1], rectilinear_mesh.y_grid))
rectilinear_mesh.z_grid = np.linspace(-10, 10, 11)

cylinder_mesh = openmc.CylindricalMesh(
    r_grid=np.linspace(0, 10, 23),
    z_grid=np.linspace(0, 1, 15)
)
cylinder_mesh.phi_grid = np.linspace(0, np.pi, 21)

spherical_mesh = openmc.SphericalMesh(
    r_grid=np.linspace(1, 10, 30),
    phi_grid=np.linspace(0, 0.8*np.pi, 25),
    theta_grid=np.linspace(0, np.pi / 2, 15),
)

MESHES = [cylinder_mesh, regular_mesh, rectilinear_mesh, spherical_mesh]

x_plane = openmc.XPlane(x0=-0.001, boundary_type='vacuum')
y_plane = openmc.YPlane(y0=-0.001, boundary_type='vacuum')
z_plane = openmc.ZPlane(z0=-0.001, boundary_type='vacuum')

SURFS = [x_plane, y_plane, z_plane]


def ids(mesh):
    if isinstance(mesh, openmc.CylindricalMesh):
        return 'cylindrical_mesh'
    elif isinstance(mesh, openmc.RegularMesh):
        return 'regular_mesh'
    elif isinstance(mesh, openmc.RectilinearMesh):
        return 'rectilinear_mesh'
    elif isinstance(mesh, openmc.SphericalMesh):
        return 'spherical_mesh'


@pytest.mark.parametrize("mesh", MESHES, ids=ids)
def test_write_data_to_vtk(mesh, tmpdir):
    # BUILD
    filename = Path(tmpdir) / "out.vtk"

    # use mesh element volumes as data to check volume-normalization ordering
    # kji (i changing fastest) orering is expected for input data
    # by using the volumes transposed as the data here, we can ensure the
    # normalization is happening correctly
    data = mesh.volumes

    # RUN
    mesh.write_data_to_vtk(filename=filename, datasets={"label1": data, "label2": data})

    # TEST
    assert filename.is_file()

    # read file
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(str(filename))
    reader.Update()

    # check name of datasets
    vtk_grid = reader.GetOutput()
    array1 = vtk_grid.GetCellData().GetArray(0)
    array2 = vtk_grid.GetCellData().GetArray(1)

    assert array1.GetName() == "label1"
    assert array2.GetName() == "label2"

    # check size of datasets
    data1 = nps.vtk_to_numpy(array1)
    data2 = nps.vtk_to_numpy(array2)
    assert data1.size == data.size
    assert data2.size == data.size

    assert all(data1 == data2)
    assert all(data1 == 1.0)


@pytest.mark.parametrize("mesh", MESHES, ids=ids)
def test_write_data_to_vtk_size_mismatch(mesh):
    """Checks that an error is raised when the size of the dataset
    doesn't match the mesh number of cells

    Parameters
    ----------
    mesh : openmc.StructuredMesh
        The mesh to test
    """
    right_size = mesh.num_mesh_cells
    data = np.random.random(right_size + 1)

    # Error message has \ in to escape characters that are otherwise recognized
    # by regex. These are needed to make the test string match the error message
    # string when using the match argument as that uses regular expression
    expected_error_msg = (
        fr"The size of the dataset 'label' \({len(data)}\) should be equal to "
        fr"the number of mesh cells \({mesh.num_mesh_cells}\)"
    )
    with pytest.raises(ValueError, match=expected_error_msg):
        mesh.write_data_to_vtk(filename="out.vtk", datasets={"label": data})

def test_write_data_to_vtk_round_trip(run_in_tmpdir):
    cmesh = openmc.CylindricalMesh(
        r_grid=(0.0, 1.0, 2.0),
        z_grid=(0.0, 2.0, 4.0, 5.0),
        phi_grid=(0.0, 3.0, 6.0),
    )

    smesh = openmc.SphericalMesh(
        r_grid=(0.0, 1.0, 2.0),
        theta_grid=(0.0, 0.5, 1.0, 2.0),
        phi_grid=(0.0, 3.0, 6.0),
    )
    rmesh = openmc.RegularMesh()
    rmesh.lower_left = (0.0, 0.0, 0.0)
    rmesh.upper_right = (1.0, 3.0, 5.0)
    rmesh.dimension = (2, 1, 6)

    for mesh in [smesh, cmesh, rmesh]:

        filename = "mesh.vtk"
        data = np.array([1.0] * 12)  # there are 12 voxels in each mesh
        mesh.write_data_to_vtk(
            filename=filename,
            datasets={"normalized": data},
            volume_normalization=True
        )

        reader = vtk.vtkStructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllFieldsOn()
        reader.Update()

        cell_data = reader.GetOutput().GetCellData()
        uniform_array = cell_data.GetArray("normalized")
        num_tuples = uniform_array.GetNumberOfTuples()
        vtk_values = [uniform_array.GetValue(i) for i in range(num_tuples)]

        # checks that the vtk cell values are equal to the data / mesh volumes
        assert np.allclose(vtk_values, data / mesh.volumes.T.flatten())

        mesh.write_data_to_vtk(
            filename=filename,
            datasets={"not_normalized": data},
            volume_normalization=False,
        )

        reader = vtk.vtkStructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllFieldsOn()
        reader.Update()

        cell_data = reader.GetOutput().GetCellData()
        uniform_array = cell_data.GetArray("not_normalized")
        num_tuples = uniform_array.GetNumberOfTuples()
        vtk_values = [uniform_array.GetValue(i) for i in range(num_tuples)]

        # checks that the vtk cell values are equal to the data
        assert np.array_equal(vtk_values, data)

def mesh_surf_id(param):
    if isinstance(param, openmc.MeshBase):
        return ids(param)
    elif isinstance(param, openmc.XPlane):
        return 'XPlane'
    elif isinstance(param, openmc.YPlane):
        return 'YPlane'
    elif isinstance(param, openmc.ZPlane):
        return 'ZPlane'


@pytest.mark.parametrize("mesh,surface", product(MESHES, SURFS), ids=mesh_surf_id)
def test_vtk_write_ordering(run_in_tmpdir, model, mesh, surface):

    tally = openmc.Tally()
    tally.scores = ['flux']
    # use the mesh on the specified tally
    mesh_filter = openmc.MeshFilter(mesh)
    tally.filters = [mesh_filter]

    model.tallies = openmc.Tallies([tally])

    # run the problem
    sp_filename = model.run()

    with openmc.StatePoint(sp_filename) as sp:
        mean = sp.tallies[tally.id].mean

    # write the data to a VTK file
    vtk_filename = 'test.vtk'
    mesh.write_data_to_vtk(vtk_filename, datasets={'mean': mean})

     # read file
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(str(vtk_filename))
    reader.Update()

    # check name of datasets
    vtk_grid = reader.GetOutput()
    array = vtk_grid.GetCellData().GetArray(0)
    vtk_data = nps.vtk_to_numpy(array)

    # convenience function for determining if a mesh
    # element has vertices in the geometry. This
    # particular geometry allows us to assume that tally results
    # in the element should be zero if none of its vertices lie in the geometry
    def in_geom(cell):
        point_ids = cell.GetPointIds()

        for i in range(point_ids.GetNumberOfIds()):
            p = vtk_grid.GetPoint(point_ids.GetId(i))
            if model.geometry.find(p):
                return True

        return False

    # reshape mean according to mesh dimensions
    mean = mean.reshape(mesh.dimension[::-1]).T
    centroid = [0.0, 0.0, 0.0]

    # check that tally and vtk array results are zero where expected
    for ijk in mesh.indices:
        ijk = tuple(n - 1 for n in ijk)
        # get the cell from the stuctured mesh object
        cell = vtk_grid.GetCell(*ijk)
        if not in_geom(cell):
            cell.GetCentroid(centroid)
            err_msg = f'IJK: {ijk} should be zero but is not. Centroid: {centroid}'
            assert mean[ijk] == 0.0, err_msg

            # need to get flat index with axes reversed due to ordering passed into the VTK file
            flat_idx = np.ravel_multi_index(tuple(ijk[::-1]), mesh.dimension[::-1])
            assert vtk_data[flat_idx] == 0.0, err_msg


def test_sphere_mesh_coordinates(run_in_tmpdir):
    mesh = openmc.SphericalMesh(
        r_grid=np.linspace(0.1, 10, 30),
        phi_grid=np.linspace(0, 1.5*np.pi, 25),
        theta_grid=np.linspace(0, np.pi / 2, 15),
    )
    # write the data to a VTK file (no data)
    vtk_filename = 'test.vtk'
    mesh.write_data_to_vtk(vtk_filename, {})

    # read file
    reader = vtk.vtkStructuredGridReader()
    reader.SetFileName(str(vtk_filename))
    reader.Update()

    vtk_grid = reader.GetOutput()

    # create a region that matches the spherical mesh description
    x = openmc.XPlane()
    z = openmc.ZPlane()
    y = openmc.YPlane()
    s = openmc.Sphere(r=10.0)

    region = +z & +y & -s | -x & -y & +z & -s

    # the VTK interface will update this list when GetCentroid is called
    centroid = np.zeros(3)

    # ensure all centroids of the sphere mesh are inside the cell region
    for i in range(vtk_grid.GetNumberOfCells()):
        # get the cell from the stuctured mesh object
        cell = vtk_grid.GetCell(i)
        cell.GetCentroid(centroid)

        # if the coordinate conversion is happening correctly,
        # every one of the cell centroids should be in the CSG region
        assert centroid in region, \
            f'Cell centroid {centroid} not in equivalent ' \
            f'CSG region for spherical mesh {mesh}'

