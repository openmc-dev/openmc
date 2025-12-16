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
    with pytest.raises(LookupError):
         statepoint.get_mesh(name='non_existent_name')

    # checks that the mesh returned is the one with the id 42
    assert statepoint.get_mesh(id=42).id == 42
    assert statepoint.get_mesh(mesh_type=openmc.RegularMesh).id == 42
    assert statepoint.get_mesh(id=42, mesh_type=openmc.RegularMesh).id == 42
    assert statepoint.get_mesh(name='custom_name').id == 42


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
