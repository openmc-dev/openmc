from math import pi

import numpy as np
import pytest
import openmc
from pathlib import Path
import tempfile
from vtkmodules import vtkIOLegacy, vtkIOXML


@pytest.mark.parametrize("val_left,val_right", [(0, 0), (-1., -1.), (2.0, 2)])
def test_raises_error_when_flat(val_left, val_right):
    """Checks that an error is raised when a mesh is flat"""
    mesh = openmc.RegularMesh()

    # Same X
    with pytest.raises(ValueError):
        mesh.lower_left = [val_left, -25, -25]
        mesh.upper_right = [val_right, 25, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [val_right, 25, 25]
        mesh.lower_left = [val_left, -25, -25]

    # Same Y
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, val_left, -25]
        mesh.upper_right = [25, val_right, 25]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, val_right, 25]
        mesh.lower_left = [-25, val_left, -25]

    # Same Z
    with pytest.raises(ValueError):
        mesh.lower_left = [-25, -25, val_left]
        mesh.upper_right = [25, 25, val_right]

    with pytest.raises(ValueError):
        mesh.upper_right = [25, 25, val_right]
        mesh.lower_left = [-25, -25, val_left]


def test_regular_mesh_bounding_box():
    mesh = openmc.RegularMesh()
    mesh.lower_left = [-2, -3, -5]
    mesh.upper_right = [2, 3, 5]
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-2, -3 ,-5))
    np.testing.assert_array_equal(bb.upper_right, (2, 3, 5))


def test_rectilinear_mesh_bounding_box():
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [0., 1., 5., 10.]
    mesh.y_grid = [-10., -5., 0.]
    mesh.z_grid = [-100., 0., 100.]
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (0., -10. ,-100.))
    np.testing.assert_array_equal(bb.upper_right, (10., 0., 100.))


def test_cylindrical_mesh_bounding_box():
    # test with mesh at origin (0, 0, 0)
    mesh = openmc.CylindricalMesh(
        r_grid=[0.1, 0.2, 0.5, 1.],
        z_grid=[0.1, 0.2, 0.4, 0.6, 1.],
        origin=(0, 0, 0)
    )
    np.testing.assert_array_equal(mesh.upper_right, (1, 1, 1))
    np.testing.assert_array_equal(mesh.lower_left, (-1, -1, 0.1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-1, -1, 0.1))
    np.testing.assert_array_equal(bb.upper_right, (1, 1, 1))

    # test with mesh at origin (3, 5, 7)
    mesh.origin = (3, 5, 7)
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 8))
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, 7.1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (2, 4, 7.1))
    np.testing.assert_array_equal(bb.upper_right, (4, 6, 8))

    # changing z grid to contain negative numbers
    mesh.z_grid = [-10, 0, 10]
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, -3))
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 17))


def test_spherical_mesh_bounding_box():
    # test with mesh at origin (0, 0, 0)
    mesh = openmc.SphericalMesh([0.1, 0.2, 0.5, 1.], origin=(0., 0., 0.))
    np.testing.assert_array_equal(mesh.upper_right, (1, 1, 1))
    np.testing.assert_array_equal(mesh.lower_left, (-1, -1, -1))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (-1, -1, -1))
    np.testing.assert_array_equal(bb.upper_right, (1, 1, 1))

    # test with mesh at origin (3, 5, 7)
    mesh.origin = (3, 5, 7)
    np.testing.assert_array_equal(mesh.upper_right, (4, 6, 8))
    np.testing.assert_array_equal(mesh.lower_left, (2, 4, 6))
    bb = mesh.bounding_box
    assert isinstance(bb, openmc.BoundingBox)
    np.testing.assert_array_equal(bb.lower_left, (2, 4, 6))
    np.testing.assert_array_equal(bb.upper_right, (4, 6, 8))


def test_SphericalMesh_initiation():
    # test defaults
    mesh = openmc.SphericalMesh(r_grid=(0, 10))
    assert (mesh.origin == np.array([0, 0, 0])).all()
    assert (mesh.r_grid == np.array([0, 10])).all()
    assert (mesh.theta_grid == np.array([0, pi])).all()
    assert (mesh.phi_grid == np.array([0, 2*pi])).all()

    # test setting on creation
    mesh = openmc.SphericalMesh(
        origin=(1, 2, 3),
        r_grid=(0, 2),
        theta_grid=(1, 3),
        phi_grid=(2, 4)
    )
    assert (mesh.origin == np.array([1, 2, 3])).all()
    assert (mesh.r_grid == np.array([0., 2.])).all()
    assert (mesh.theta_grid == np.array([1, 3])).all()
    assert (mesh.phi_grid == np.array([2, 4])).all()

    # test attribute changing
    mesh.r_grid = (0, 11)
    assert (mesh.r_grid == np.array([0., 11.])).all()

    # test invalid r_grid values
    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[1, 1])

    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[0])

    # test invalid theta_grid values
    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[1, 2], theta_grid=[1, 1])

    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[1, 2], theta_grid=[0])

    # test invalid phi_grid values
    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[1, 2], phi_grid=[1, 1])

    with pytest.raises(ValueError):
        openmc.SphericalMesh(r_grid=[1, 2], phi_grid=[0])

    # waffles and pancakes are unfortunately not valid radii
    with pytest.raises(TypeError):
        openmc.SphericalMesh(('🧇', '🥞'))


def test_CylindricalMesh_initiation():
    # test defaults
    mesh = openmc.CylindricalMesh(r_grid=(0, 10), z_grid=(0, 10))
    assert (mesh.origin == np.array([0, 0, 0])).all()
    assert (mesh.r_grid == np.array([0, 10])).all()
    assert (mesh.phi_grid == np.array([0, 2*pi])).all()
    assert (mesh.z_grid == np.array([0, 10])).all()

    # test setting on creation
    mesh = openmc.CylindricalMesh(
        origin=(1, 2, 3),
        r_grid=(0, 2),
        z_grid=(1, 3),
        phi_grid=(2, 4)
    )
    assert (mesh.origin == np.array([1, 2, 3])).all()
    assert (mesh.r_grid == np.array([0., 2.])).all()
    assert (mesh.z_grid == np.array([1, 3])).all()
    assert (mesh.phi_grid == np.array([2, 4])).all()

    # test attribute changing
    mesh.r_grid = (0., 10.)
    assert (mesh.r_grid == np.array([0, 10.])).all()
    mesh.z_grid = (0., 4.)
    assert (mesh.z_grid == np.array([0, 4.])).all()

    # waffles and pancakes are unfortunately not valid radii
    with pytest.raises(TypeError):
        openmc.SphericalMesh(('🧇', '🥞'))


def test_invalid_cylindrical_mesh_errors():
    # Test invalid r_grid values
    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[5, 1], phi_grid=[0, pi], z_grid=[0, 10])

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[1, 2, 4, 3], phi_grid=[0, pi], z_grid=[0, 10])

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[1], phi_grid=[0, pi], z_grid=[0, 10])

    # Test invalid phi_grid values
    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[0, 1, 2], phi_grid=[-1, 3], z_grid=[0, 10])

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(
            r_grid=[0, 1, 2],
            phi_grid=[0, 2*pi + 0.1],
            z_grid=[0, 10]
        )

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[0, 1, 2], phi_grid=[pi], z_grid=[0, 10])

    # Test invalid z_grid values
    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[0, 1, 2], phi_grid=[0, pi], z_grid=[5])

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[0, 1, 2], phi_grid=[0, pi], z_grid=[5, 1])

    with pytest.raises(ValueError):
        openmc.CylindricalMesh(r_grid=[1, 2, 4, 3], phi_grid=[0, pi], z_grid=[0, 10, 5])


def test_centroids():
    # regular mesh
    mesh = openmc.RegularMesh()
    mesh.lower_left = (1., 2., 3.)
    mesh.upper_right = (11., 12., 13.)
    mesh.dimension = (1, 1, 1)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [6., 7., 8.])

    # rectilinear mesh
    mesh = openmc.RectilinearMesh()
    mesh.x_grid = [1., 11.]
    mesh.y_grid = [2., 12.]
    mesh.z_grid = [3., 13.]
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [6., 7., 8.])

    # cylindrical mesh
    mesh = openmc.CylindricalMesh(r_grid=(0, 10), z_grid=(0, 10), phi_grid=(0, np.pi))
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [0.0, 5.0, 5.0])
    # ensure that setting an origin is handled correctly
    mesh.origin = (5.0, 0, -10)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [5.0, 5.0, -5.0])

    # spherical mesh, single element xyz-positive octant
    mesh = openmc.SphericalMesh(r_grid=[0, 10], theta_grid=[0, 0.5*np.pi], phi_grid=[0, np.pi])
    x = 5.*np.cos(0.5*np.pi)*np.sin(0.25*np.pi)
    y = 5.*np.sin(0.5*np.pi)*np.sin(0.25*np.pi)
    z = 5.*np.sin(0.25*np.pi)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [x, y, z])

    mesh.origin = (-5.0, -5.0, 5.0)
    np.testing.assert_array_almost_equal(mesh.centroids[0, 0, 0], [x-5.0, y-5.0, z+5.0])


@pytest.mark.parametrize('mesh_type', ('regular', 'rectilinear', 'cylindrical', 'spherical'))
def test_mesh_vertices(mesh_type):

    ijk = (2, 3, 2)

    # create a new mesh object
    if mesh_type == 'regular':
        mesh = openmc.RegularMesh()
        ll = np.asarray([0.]*3)
        width = np.asarray([0.5]*3)
        mesh.lower_left = ll
        mesh.width = width
        mesh.dimension = (5, 7, 9)

        # spot check that an element has the correct vertex coordinates asociated with it
        # (using zero-indexing here)
        exp_i_j_k = ll + np.asarray(ijk, dtype=float) * width
        np.testing.assert_equal(mesh.vertices[ijk], exp_i_j_k)

        # shift the mesh using the llc
        shift  = np.asarray((3.0, 6.0, 10.0))
        mesh.lower_left += shift
        np.testing.assert_equal(mesh.vertices[ijk], exp_i_j_k+shift)
    elif mesh_type == 'rectilinear':
        mesh = openmc.RectilinearMesh()
        w = np.asarray([0.5] * 3)
        ll = np.asarray([0.]*3)
        dims = (5, 7, 9)
        mesh.x_grid = np.linspace(ll[0], w[0]*dims[0], dims[0])
        mesh.y_grid = np.linspace(ll[1], w[1]*dims[1], dims[1])
        mesh.z_grid = np.linspace(ll[2], w[2]*dims[2], dims[2])
        exp_vert = np.asarray((mesh.x_grid[2], mesh.y_grid[3], mesh.z_grid[2]))
        np.testing.assert_equal(mesh.vertices[ijk], exp_vert)
    elif mesh_type == 'cylindrical':
        r_grid = np.linspace(0, 5, 10)
        z_grid = np.linspace(-10, 10, 20)
        phi_grid = np.linspace(0, 2*np.pi, 8)
        mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid, phi_grid=phi_grid)
        exp_vert = np.asarray((mesh.r_grid[2], mesh.phi_grid[3], mesh.z_grid[2]))
        np.testing.assert_equal(mesh.vertices_cylindrical[ijk], exp_vert)
    elif mesh_type == 'spherical':
        r_grid = np.linspace(0, 13, 14)
        theta_grid = np.linspace(0, np.pi, 11)
        phi_grid = np.linspace(0, 2*np.pi, 7)
        mesh = openmc.SphericalMesh(r_grid=r_grid, theta_grid=theta_grid, phi_grid=phi_grid)
        exp_vert = np.asarray((mesh.r_grid[2], mesh.theta_grid[3], mesh.phi_grid[2]))
        np.testing.assert_equal(mesh.vertices_spherical[ijk], exp_vert)


def test_CylindricalMesh_get_indices_at_coords():
    # default origin (0, 0, 0) and default phi grid (0, 2*pi)
    mesh = openmc.CylindricalMesh(r_grid=(0, 5, 10), z_grid=(0, 5, 10))
    assert mesh.get_indices_at_coords([1, 0, 1]) == (0, 0, 0)
    assert mesh.get_indices_at_coords([6, 0, 1]) == (1, 0, 0)
    assert mesh.get_indices_at_coords([9, 0, 1]) == (1, 0, 0)
    assert mesh.get_indices_at_coords([0, 6, 0]) == (1, 0, 0)
    assert mesh.get_indices_at_coords([0, 9, 6]) == (1, 0, 1)
    assert mesh.get_indices_at_coords([-2, -2, 9]) == (0, 0, 1)

    with pytest.raises(ValueError):
        assert mesh.get_indices_at_coords([8, 8, 1])  # resulting r value to large
    with pytest.raises(ValueError):
        assert mesh.get_indices_at_coords([-8, -8, 1])  # resulting r value to large
    with pytest.raises(ValueError):
        assert mesh.get_indices_at_coords([1, 0, -1])  # z value below range
    with pytest.raises(ValueError):
        assert mesh.get_indices_at_coords([1, 0, 11])  # z value above range

    assert mesh.get_indices_at_coords([1, 1, 1]) == (0, 0, 0)

    # negative range on z grid
    mesh = openmc.CylindricalMesh(
        r_grid=(0, 5, 10),
        phi_grid=(0, 0.5 * pi, pi, 1.5 * pi, 1.9 * pi),
        z_grid=(-5, 0, 5, 10),
    )
    assert mesh.get_indices_at_coords([1, 1, 1]) == (0, 0, 1)  # first angle quadrant
    assert mesh.get_indices_at_coords([2, 2, 6]) == (0, 0, 2)  # first angle quadrant
    assert mesh.get_indices_at_coords([-2, 0.1, -1]) == (0, 1, 0)  # second angle quadrant
    assert mesh.get_indices_at_coords([-2, -0.1, -1]) == (0, 2, 0)  # third angle quadrant
    assert mesh.get_indices_at_coords([2, -0.9, -1]) == (0, 3, 0)  # forth angle quadrant

    with pytest.raises(ValueError):
        assert mesh.get_indices_at_coords([2, -0.1, 1])  # outside of phi range

    # origin of mesh not default
    mesh = openmc.CylindricalMesh(
        r_grid=(0, 5, 10),
        phi_grid=(0, 0.5 * pi, pi, 1.5 * pi, 1.9 * pi),
        z_grid=(-5, 0, 5, 10),
        origin=(100, 200, 300),
    )
    assert mesh.get_indices_at_coords([101, 201, 301]) == (0, 0, 1)  # first angle quadrant
    assert mesh.get_indices_at_coords([102, 202, 306]) == (0, 0, 2)  # first angle quadrant
    assert mesh.get_indices_at_coords([98, 200.1, 299]) == (0, 1, 0)  # second angle quadrant
    assert mesh.get_indices_at_coords([98, 199.9, 299]) == (0, 2, 0)  # third angle quadrant
    assert mesh.get_indices_at_coords([102, 199.1, 299]) == (0, 3, 0)  # forth angle quadrant


def test_mesh_name_roundtrip(run_in_tmpdir):

    mesh = openmc.RegularMesh()
    mesh.name = 'regular-mesh'
    mesh.lower_left = (-1, -1, -1)
    mesh.width = (1, 1, 1)
    mesh.dimension = (1, 1, 1)

    mesh_filter = openmc.MeshFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['flux']

    openmc.Tallies([tally]).export_to_xml()

    xml_tallies = openmc.Tallies.from_xml()

    mesh = xml_tallies[0].find_filter(openmc.MeshFilter).mesh
    assert mesh.name == 'regular-mesh'


def test_umesh_roundtrip(run_in_tmpdir, request):
    umesh = openmc.UnstructuredMesh(request.path.parent / 'test_mesh_tets.e', 'moab')
    umesh.output = True

    # create a tally using this mesh
    mf = openmc.MeshFilter(umesh)
    tally = openmc.Tally()
    tally.filters = [mf]
    tally.scores = ['flux']

    tallies = openmc.Tallies([tally])
    tallies.export_to_xml()

    xml_tallies = openmc.Tallies.from_xml()
    xml_tally = xml_tallies[0]
    xml_mesh = xml_tally.filters[0].mesh

    assert umesh.id == xml_mesh.id


@pytest.mark.parametrize('export_type', ('.vtk','.vtu'))
def test_umesh(request, export_type):
    """Performs a minimal UnstructuredMesh simulation, reads in the resulting
    statepoint file and writes the mesh data to vtk and vtkhdf files. It is
    necessary to read in the unstructured mesh from a statepoint file to ensure
    it has all the required attributes
    """
    surf1 = openmc.Sphere(r=1.0, boundary_type="vacuum")
    material1 = openmc.Material()
    material1.add_components(components={"H" : 1.0})
    material1.set_density('g/cm3', 1.0)

    materials = openmc.Materials([material1])
    cell1 = openmc.Cell(region=-surf1, fill=material1)
    geometry = openmc.Geometry([cell1])

    umesh = openmc.UnstructuredMesh(
       filename= request.path.parent.parent
        / "regression_tests/unstructured_mesh/test_mesh_dagmc_tets.vtk",
       library= "moab",
       mesh_id = 1
    )
    # setting ID to make it easier to get the mesh from the statepoint later
    mesh_filter = openmc.MeshFilter(umesh)

    # Create flux mesh tally to score alpha production
    mesh_tally = openmc.Tally(name="test_tally")
    mesh_tally.filters = [mesh_filter]
    mesh_tally.scores = ["total"]

    tallies = openmc.Tallies([mesh_tally])

    source = openmc.IndependentSource()

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = 2
    settings.particles = 100
    settings.source = source

    model = openmc.Model(
        materials=materials, geometry=geometry, settings=settings, tallies=tallies
    )

    with tempfile.TemporaryDirectory() as tmpdir:

        statepoint_file = model.run(cwd=tmpdir)

        statepoint = openmc.StatePoint(statepoint_file)
        tally: openmc.Tally = statepoint.get_tally(name="test_tally")
        umesh_from_sp: openmc.UnstructuredMesh = statepoint.meshes[1]

        datasets={
            "mean": tally.mean.squeeze(),
            "std_dev": tally.std_dev.squeeze()
        }

        filename = Path(tmpdir) / ("test_mesh" + export_type)
        umesh_from_sp.write_data_to_vtk(datasets=datasets, filename=filename)

        assert Path(filename).exists()

        if export_type == ".vtk":
            reader = vtkIOLegacy.vtkGenericDataObjectReader()
            reader.SetFileName(str(filename))
            reader.Update()
            assert reader.GetOutput().GetCellData().GetArray("mean")

        elif export_type == ".vtk":
            reader = vtkIOXML.vtkGenericDataObjectReader()
            reader.SetFileName(str(filename))
            reader.Update()
            assert reader.GetOutput().GetCellData().GetArray("mean")


def test_mesh_get_homogenized_materials():
    """Test the get_homogenized_materials method"""
    # Simple model with 1 cm of Fe56 next to 1 cm of H1
    fe = openmc.Material()
    fe.add_nuclide('Fe56', 1.0)
    fe.set_density('g/cm3', 5.0)
    h = openmc.Material()
    h.add_nuclide('H1', 1.0)
    h.set_density('g/cm3', 1.0)

    x0 = openmc.XPlane(-1.0, boundary_type='vacuum')
    x1 = openmc.XPlane(0.0)
    x2 = openmc.XPlane(1.0)
    x3 = openmc.XPlane(2.0, boundary_type='vacuum')
    cell1 = openmc.Cell(fill=fe, region=+x0 & -x1)
    cell2 = openmc.Cell(fill=h, region=+x1 & -x2)
    cell_empty = openmc.Cell(region=+x2 & -x3)
    model = openmc.Model(geometry=openmc.Geometry([cell1, cell2, cell_empty]))
    model.settings.particles = 1000
    model.settings.batches = 10

    mesh = openmc.RegularMesh()
    mesh.lower_left = (-1., -1., -1.)
    mesh.upper_right = (1., 1., 1.)
    mesh.dimension = (3, 1, 1)
    m1, m2, m3 = mesh.get_homogenized_materials(model, n_samples=1_000_000)

    # Left mesh element should be only Fe56
    assert m1.get_mass_density('Fe56') == pytest.approx(5.0)

    # Middle mesh element should be 50% Fe56 and 50% H1
    assert m2.get_mass_density('Fe56') == pytest.approx(2.5, rel=1e-2)
    assert m2.get_mass_density('H1') == pytest.approx(0.5, rel=1e-2)

    # Right mesh element should be only H1
    assert m3.get_mass_density('H1') == pytest.approx(1.0)

    mesh_void = openmc.RegularMesh()
    mesh_void.lower_left = (0.5, 0.5, -1.)
    mesh_void.upper_right = (1.5, 1.5, 1.)
    mesh_void.dimension = (1, 1, 1)
    m4, = mesh_void.get_homogenized_materials(model, n_samples=1_000_000)

    # Mesh element that overlaps void should have half density
    assert m4.get_mass_density('H1') == pytest.approx(0.5, rel=1e-2)

    # If not including void, density of homogenized material should be same as
    # original material
    m5, = mesh_void.get_homogenized_materials(
        model, n_samples=1000, include_void=False)
    assert m5.get_mass_density('H1') == pytest.approx(1.0)
