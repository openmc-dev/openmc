import openmc


def test_tally_trigger(run_in_tmpdir):
    pincell = openmc.examples.pwr_pin_cell()

    # create a tally filter on the materials
    mat_filter = openmc.MaterialFilter(pincell.materials)

    # create a tally with triggers applied
    tally = openmc.Tally()
    tally.filters = [mat_filter]
    tally.scores = ["scatter"]

    trigger = openmc.Trigger("rel_err", 0.05)
    trigger.scores = ["scatter"]

    tally.triggers = [trigger]

    pincell.tallies = [tally]

    pincell.settings.trigger_active = True
    pincell.settings.trigger_max_batches = 100
    pincell.settings.trigger_batch_interval = 5

    sp_file = pincell.run()
    with openmc.StatePoint(sp_file) as sp:
        expected_realizations = sp.n_realizations

    # adding other scores to the tally should not change the
    # number of batches required to satisfy the trigger
    tally.scores = ["total", "absorption", "scatter"]

    sp_file = pincell.run()

    with openmc.StatePoint(sp_file) as sp:
        realizations = sp.n_realizations

    assert realizations == expected_realizations


def test_tally_trigger_null_score(run_in_tmpdir):
    pincell = openmc.examples.pwr_pin_cell()

    # create a tally filter on the materials
    mat_filter = openmc.MaterialFilter(pincell.materials)

    # apply a tally with a score that be tallied in this model
    tally = openmc.Tally()
    tally.filters = [mat_filter]
    tally.scores = ["pair-production"]

    trigger = openmc.Trigger("rel_err", 0.05)
    trigger.scores = ["pair-production"]

    tally.triggers = [trigger]

    pincell.tallies = [tally]

    pincell.settings.trigger_active = True
    pincell.settings.trigger_max_batches = 50
    pincell.settings.trigger_batch_interval = 5

    sp_file = pincell.run()

    with openmc.StatePoint(sp_file) as sp:
        # verify that the tally mean is zero
        tally_out = sp.get_tally(id=tally.id)
        assert all(tally_out.mean == 0.0)

        # we expect that this simulation will run
        # up to the max allowed batches
        total_batches = sp.n_realizations + sp.n_inactive
        assert total_batches == pincell.settings.trigger_max_batches


def test_tally_trigger_zero_ignored(run_in_tmpdir):
    pincell = openmc.examples.pwr_pin_cell()

    # create an energy filter below and around the O-16(n,p) threshold (1.02e7 eV)
    e_filter = openmc.EnergyFilter([0.0, 1e7, 2e7])

    # create a tally with triggers applied
    tally = openmc.Tally()
    tally.filters = [e_filter]
    tally.scores = ["(n,p)"]
    tally.nuclides = ["O16"]

    # 100% relative error: should be immediately satisfied in nonzero bin
    trigger = openmc.Trigger("rel_err", 1.0)
    trigger.scores = ["(n,p)"]
    trigger.ignore_zeros = True

    tally.triggers = [trigger]

    pincell.tallies = [tally]

    pincell.settings.particles = 1000  # we need a few more particles for this
    pincell.settings.trigger_active = True
    pincell.settings.trigger_max_batches = 50
    pincell.settings.trigger_batch_interval = 20

    sp_file = pincell.run()

    with openmc.StatePoint(sp_file) as sp:
        # verify that the first bin is zero and the second is nonzero
        tally_out = sp.get_tally(id=tally.id)
        below, above = tally_out.mean.squeeze()
        assert below == 0.0, "Tally events observed below expected threshold"
        assert above > 0, "No tally events observed. Test with more particles."

        # we expect that the trigger fires before max batches are hit
        total_batches = sp.n_realizations + sp.n_inactive
        assert total_batches < pincell.settings.trigger_max_batches


def test_trigger_he3_production(run_in_tmpdir):
    li6 = openmc.Material()
    li6.set_density("g/cm3", 1.0)
    li6.add_nuclide("Li6", 1.0)

    sph = openmc.Sphere(r=20, boundary_type="vacuum")
    outer_cell = openmc.Cell(fill=li6, region=-sph)
    model = openmc.Model()
    model.geometry = openmc.Geometry([outer_cell])
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(14.1e6)
    )
    model.settings.batches = 10
    model.settings.particles = 100
    model.settings.run_mode = "fixed source"
    model.settings.trigger_active = True
    model.settings.trigger_batch_interval = 10
    model.settings.trigger_max_batches = 30

    # Define tally with trigger
    trigger = openmc.Trigger(trigger_type="rel_err", threshold=0.0001)
    trigger.scores = ["He3-production"]
    he3_production_tally = openmc.Tally()
    he3_production_tally.scores = ["He3-production"]
    he3_production_tally.triggers = [trigger]
    model.tallies = openmc.Tallies([he3_production_tally])

    # Run model to verify that trigger works
    model.run()
