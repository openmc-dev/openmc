import openmc
from pytest import approx


def test_nuclide_heating(run_in_tmpdir):
    mat = openmc.Material()
    mat.add_nuclide("Li6", 0.5)
    mat.add_nuclide("Li7", 0.5)
    mat.set_density("g/cm3", 1.0)

    sphere = openmc.Sphere(r=20, boundary_type="reflective")
    inside_sphere = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([inside_sphere])

    model.settings.particles = 1000
    model.settings.batches = 1
    model.settings.photon_transport = True
    model.settings.electron_treatment = "ttb"
    model.settings.cutoff = {"energy_photon": 1000}
    model.settings.run_mode = "fixed source"
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(10.0e6),
        particle="photon"
    )

    # Create two tallies, one with heating by nuclide and one with total heating
    tally1 = openmc.Tally()
    tally1.scores = ["heating"]
    tally1.nuclides = mat.get_nuclides()
    tally2 = openmc.Tally()
    tally2.scores = ["heating"]
    model.tallies = [tally1, tally2]

    # Run the model
    model.run(apply_tally_results=True)

    # Make sure the heating results are consistent
    assert tally1.mean.sum() == approx(tally2.mean.sum())
