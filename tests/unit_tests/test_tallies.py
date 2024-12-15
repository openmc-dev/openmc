import numpy as np

import openmc


def test_xml_roundtrip(run_in_tmpdir):
    # Create a tally with all possible gizmos
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10.0, -10.0, -10.0)
    mesh.upper_right = (
        10.0,
        10.0,
        10.0,
    )
    mesh.dimension = (5, 5, 5)
    mesh_filter = openmc.MeshFilter(mesh)
    meshborn_filter = openmc.MeshBornFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter, meshborn_filter]
    tally.nuclides = ["U235", "I135", "Li6"]
    tally.scores = ["total", "fission", "heating"]
    tally.derivative = openmc.TallyDerivative(
        variable="nuclide_density", material=1, nuclide="Li6"
    )
    tally.triggers = [openmc.Trigger("rel_err", 0.025)]
    tally.triggers[0].scores = ["total", "fission"]
    tallies = openmc.Tallies([tally])

    # Roundtrip through XML and make sure we get what we started with
    tallies.export_to_xml()
    new_tallies = openmc.Tallies.from_xml()
    assert len(new_tallies) == 1
    new_tally = new_tallies[0]
    assert new_tally.id == tally.id
    assert len(new_tally.filters) == 2
    assert isinstance(new_tally.filters[0], openmc.MeshFilter)
    assert np.allclose(new_tally.filters[0].mesh.lower_left, mesh.lower_left)
    assert isinstance(new_tally.filters[1], openmc.MeshBornFilter)
    assert np.allclose(new_tally.filters[1].mesh.lower_left, mesh.lower_left)
    assert new_tally.nuclides == tally.nuclides
    assert new_tally.scores == tally.scores
    assert new_tally.derivative.variable == tally.derivative.variable
    assert new_tally.derivative.material == tally.derivative.material
    assert new_tally.derivative.nuclide == tally.derivative.nuclide
    assert len(new_tally.triggers) == 1
    assert new_tally.triggers[0].trigger_type == tally.triggers[0].trigger_type
    assert new_tally.triggers[0].threshold == tally.triggers[0].threshold
    assert new_tally.triggers[0].scores == tally.triggers[0].scores
