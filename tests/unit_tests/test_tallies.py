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


def test_tally_equivalence():
    tally_a = openmc.Tally()
    tally_b = openmc.Tally(tally_id=tally_a.id)

    tally_a.name = 'new name'
    assert tally_a != tally_b
    tally_b.name = tally_a.name
    assert tally_a == tally_b

    assert tally_a == tally_b
    ef_a = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    ef_b = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])

    tally_a.filters = [ef_a]
    assert tally_a != tally_b
    tally_b.filters = [ef_b]
    assert tally_a == tally_b

    tally_a.scores = ['flux', 'absorption', 'fission', 'scatter']
    assert tally_a != tally_b
    tally_b.scores = ['flux', 'absorption', 'fission', 'scatter']
    assert tally_a == tally_b

    tally_a.nuclides = []
    tally_b.nuclides = []
    assert tally_a == tally_b

    tally_a.nuclides = ['total']
    assert tally_a == tally_b

    # a tally with an estimator set to None is equal to
    # a tally with an estimator specified
    tally_a.estimator = 'collision'
    assert tally_a == tally_b
    tally_b.estimator = 'collision'
    assert tally_a == tally_b

    tally_a.multiply_density = False
    assert tally_a != tally_b
    tally_b.multiply_density = False
    assert tally_a == tally_b

    trigger_a = openmc.Trigger('rel_err', 0.025)
    trigger_b = openmc.Trigger('rel_err', 0.025)

    tally_a.triggers = [trigger_a]
    assert tally_a != tally_b
    tally_b.triggers = [trigger_b]
    assert tally_a == tally_b


def test_tally_application(run_in_tmpdir):
    model = openmc.examples.pwr_pin_cell()

    # Create a tally with most possible gizmos
    tally = openmc.Tally()
    tally.name = 'test tally'
    ef = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    bbox = model.geometry.bounding_box
    mesh = openmc.RegularMesh()
    mesh.lower_left = bbox[0][:2]
    mesh.upper_right = bbox[1][:2]
    mesh.dimension = (2, 2)
    mf = openmc.MeshFilter(mesh)
    tally.filters = [ef, mf]
    tally.scores = ['flux', 'absorption', 'fission', 'scatter']
    model.tallies = [tally]

    # run the simulation and apply retsults
    sp_file = model.run(apply_tally_results=True)
    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.mean == tally.mean).all()
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert sp_tally.nuclides == tally.nuclides

