import openmc
import pytest


def test_get_mesh():
    """tests getting the correct mesh from the statepoint file"""

    model = openmc.Model()

    h1 = openmc.Material()
    h1.add_nuclide("H1", 1.0)
    h1.set_density("g/cm3", 1.0)
    model.materials = openmc.Materials([h1])

    sphere = openmc.Sphere(r=10, boundary_type="vacuum")
    cell = openmc.Cell(fill=h1, region=-sphere)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = "fixed source"
    model.settings.particles = 10
    model.settings.batches = 2

    model.settings.source = openmc.IndependentSource()

    mesh = openmc.RegularMesh()
    mesh.id = 42
    mesh.name = "custom_name"
    mesh.dimension = (2, 2, 1)
    mesh.lower_left = (-10, -10, -10)
    mesh.upper_right = (10, 10, 10)

    tally = openmc.Tally()
    tally.scores = ["flux"]
    tally.filters = [openmc.MeshFilter(mesh)]
    model.tallies = openmc.Tallies([tally])

    statepoint_fn = model.run()

    statepoint = openmc.StatePoint(statepoint_fn)

    # checks that the mesh is not found in the statepoint file
    with pytest.raises(LookupError):
        statepoint.get_mesh(id=999)
    with pytest.raises(LookupError):
        statepoint.get_mesh(mesh_type=openmc.CylindricalMesh)
    with pytest.raises(LookupError):
        statepoint.get_mesh(id=42, mesh_type=openmc.CylindricalMesh)
    with pytest.raises(LookupError):
        statepoint.get_mesh(id=999, mesh_type=openmc.RegularMesh)

    # checks that the mesh returned is the one with the id 42
    assert statepoint.get_mesh(id=42).id == 42
    assert statepoint.get_mesh(mesh_type=openmc.RegularMesh).id == 42
    assert statepoint.get_mesh(id=42, mesh_type=openmc.RegularMesh).id == 42

    # TODO add this test when in a bug fix PR as fixed the problem that he mesh is
    # returned from the statepoint has no name set
    # assert statepoint.get_mesh(name='custom_name').id == 42
