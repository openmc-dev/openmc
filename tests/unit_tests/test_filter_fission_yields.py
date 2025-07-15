import openmc
import numpy as np


def test_fissionyieldsfilter(run_in_tmpdir):
    steel = openmc.Material(name='U')
    steel.set_density('g/cm3', 8.00)
    steel.add_nuclide('U235', 1.0)

    sphere = openmc.Sphere(r=50.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=steel)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 100
    model.settings.batches = 10

    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1e6),
    )
    model.settings.run_mode = "eigenvalue"

    fissionyields_filter = openmc.FissionYieldsFilter(['Xe135','I135'])

    tally = openmc.Tally()
    tally.filters = [fissionyields_filter]
    tally.scores = ['fission']
    model.tallies = openmc.Tallies([tally])

    model.export_to_model_xml()

    # Run OpenMC
    model.run(apply_tally_results=True)

    assert np.all(tally.mean.reshape(2)>0)
    


