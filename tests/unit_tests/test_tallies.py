import numpy as np

import openmc


def test_xml_roundtrip(run_in_tmpdir):
    # Create a tally with all possible gizmos
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10., -10., -10.)
    mesh.upper_right = (10., 10., 10.,)
    mesh.dimension = (5, 5, 5)
    mesh_filter = openmc.MeshFilter(mesh)
    meshborn_filter = openmc.MeshBornFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter, meshborn_filter]
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
    assert new_tally.multiply_density == tally.multiply_density

def test_tally_application(run_in_tmpdir):

    model = openmc.examples.pwr_pin_cell()

    # Create a tally with most possible gizmos
    tally = openmc.Tally()
    tally.name = 'test tally'
    ef = openmc.EnergyFilter.from_group_structure('CASMO-40')
    bbox = model.geometry.bounding_box
    mesh = openmc.RegularMesh()
    mesh.lower_left = bbox[0][:2]
    mesh.upper_right = bbox[1][:2]
    mesh.dimension = (10, 10)
    mf = openmc.MeshFilter(mesh)
    tally.filters = [ef, mf]
    tally.scores = ['flux', 'absorption', 'fission', 'scatter']
    model.tallies = [tally]

    # run the simulation and apply retsults
    sp_file = model.run(apply_tally_results=True)

    with openmc.StatePoint(sp_file) as sp:
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.mean == tally.mean).all()
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert sp_tally.nuclides == tally.nuclides

