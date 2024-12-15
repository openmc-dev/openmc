from random import uniform

import pytest
import openmc


def inf_medium_model(cutoff_energy, source_energy):
    """Infinite medium problem with a monoenergetic photon source"""
    model = openmc.Model()

    m = openmc.Material()
    m.add_nuclide("Zr90", 1.0)
    m.set_density("g/cm3", 1.0)

    sph = openmc.Sphere(r=100.0, boundary_type="reflective")
    cell = openmc.Cell(fill=m, region=-sph)
    model.geometry = openmc.Geometry([cell])

    model.settings.run_mode = "fixed source"
    model.settings.source = openmc.IndependentSource(
        particle="photon",
        energy=openmc.stats.Discrete([source_energy], [1.0]),
    )
    model.settings.particles = 100
    model.settings.batches = 10
    model.settings.cutoff = {"energy_photon": cutoff_energy}

    tally_flux = openmc.Tally(name="flux")
    tally_flux.filters = [
        openmc.EnergyFilter([0.0, cutoff_energy, source_energy]),
        openmc.ParticleFilter(["photon"]),
    ]
    tally_flux.scores = ["flux"]
    tally_heating = openmc.Tally(name="heating")
    tally_heating.scores = ["heating"]
    model.tallies = openmc.Tallies([tally_flux, tally_heating])

    return model


def test_energy_cutoff(run_in_tmpdir):
    # Pick a random cutoff energy between 1 and 5 keV
    cutoff_energy = uniform(1e3, 5e3)

    # Pick a random source energy some factor higher than cutoff energy
    source_energy = uniform(10, 20) * cutoff_energy

    # Create model and run simulation
    model = inf_medium_model(cutoff_energy, source_energy)
    statepoint_path = model.run()

    # Get resulting flux and heating values
    with openmc.StatePoint(statepoint_path) as sp:
        flux = sp.get_tally(name="flux").mean.ravel()
        heating = sp.get_tally(name="heating").mean.ravel()

    # There should be no flux below the cutoff energy (first bin in the tally)
    assert flux[0] == 0.0
    assert flux[1] > 0.0

    # Despite killing particles below the cutoff, the total heating should be
    # equal to the source energy
    assert heating[0] == pytest.approx(source_energy)
