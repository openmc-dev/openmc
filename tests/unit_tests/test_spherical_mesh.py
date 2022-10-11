import openmc
import numpy as np

def test_origin_read_write_to_xml():
    """Tests that the origin attribute can be written and read back to XML
    """
    # build
    mesh = openmc.SphericalMesh()
    mesh.phi_grid = [1, 2, 3]
    mesh.theta_grid = [1, 2, 3]
    mesh.r_grid = [1, 2, 3]

    mesh.origin = [0.1, 0.2, 0.3]

    tally = openmc.Tally()

    mesh_filter = openmc.MeshFilter(mesh)
    tally.filters.append(mesh_filter)

    tally.scores.append("heating")

    tallies = openmc.Tallies([tally])

    tallies.export_to_xml()

    # read back
    new_tallies = openmc.Tallies.from_xml()
    new_tally = new_tallies[0]
    new_mesh = new_tally.filters[0].mesh
    assert np.allclose(new_mesh.origin, mesh.origin)