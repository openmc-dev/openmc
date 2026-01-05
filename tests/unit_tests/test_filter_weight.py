import openmc
import numpy as np


def test_weightfilter(run_in_tmpdir):
    steel = openmc.Material(name='Stainless Steel')
    steel.set_density('g/cm3', 8.00)
    steel.add_nuclide('Fe56', 1.0)

    sphere = openmc.Sphere(r=50.0, boundary_type='vacuum')
    cell = openmc.Cell(region=-sphere, fill=steel)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.particles = 100
    model.settings.batches = 10

    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(14e6),
    )
    model.settings.run_mode = "fixed source"

    radius = list(range(1, 50))
    sphere_mesh = openmc.SphericalMesh(radius)
    mesh_filter = openmc.MeshFilter(sphere_mesh)
    weight_filter = openmc.WeightFilter(
        [0.999, 0.9999, 0.99999, 0.999999, 1.0, 1.000001 ,1.00001, 1.0001, 1.001]
    )

    tally = openmc.Tally()
    tally.filters = [mesh_filter, weight_filter]
    tally.estimator = 'analog'
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    model.run(apply_tally_results=True)

    # Get current binned by mu
    neutron_flux = tally.mean.reshape(48, 8)

    # All contributions should show up in the fourth bin
    assert np.all(neutron_flux[:, 3] != 0.0)
    neutron_flux[:, 3] = 0.0
    assert np.all(neutron_flux == 0.0)
