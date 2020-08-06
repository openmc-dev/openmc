from math import pi, log
from random import uniform

import openmc.deplete
import openmc
import pytest


@pytest.fixture
def model():
    """Sphere of single nuclide"""
    model = openmc.model.Model()

    w = openmc.Material(name='tungsten')
    w.add_nuclide('W186', 1.0)
    w.set_density('g/cm3', 19.3)
    w.depletable = True

    r = uniform(1.0, 10.0)
    w.volume = 4/3 * pi * r**3

    surf = openmc.Sphere(r=r, boundary_type='vacuum')
    cell = openmc.Cell(fill=w, region=-surf)
    model.geometry = openmc.Geometry([cell])

    model.settings.batches = 10
    model.settings.particles = 1000
    model.settings.source = openmc.Source(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0])
    )
    model.settings.run_mode = 'fixed source'

    rx_tally = openmc.Tally()
    rx_tally.scores = ['(n,gamma)']
    model.tallies.append(rx_tally)

    return model


def test_activation(run_in_tmpdir, model):
    # Determine (n.gamma) reaction rate using initial run
    sp = model.run()
    with openmc.StatePoint(sp) as sp:
        tally = sp.tallies[1]
        capture_rate = tally.mean.flat[0]

    # Create one-nuclide depletion chain
    chain = openmc.deplete.Chain()
    w186 = openmc.deplete.Nuclide('W186')
    w186.add_reaction('(n,gamma)', None, 0.0, 1.0)
    chain.add_nuclide(w186)
    chain.export_to_xml('test_chain.xml')

    # Create transport operator
    op = openmc.deplete.Operator(
        model.geometry, model.settings, 'test_chain.xml',
        normalization_mode="source-rate"
    )

    # To determine the source rate necessary to reduce W186 density in half, we
    # start with the single-nuclide transmutation equation:
    #
    #                 dn/dt = -f * sigma * phi * n
    #                  n(t) = n0 * exp(-f * sigma * phi * t)
    #
    # where f is the source rate. The capture rate, r, is sigma * phi * n0,
    # meaning that:
    #
    #                  n(t) = n0 * exp(-f * r * t / n0)
    #
    # To reduce the density by half, we would need:
    #
    #               n(t)/n0 = exp(-f * r * t / n0) = 1/2
    #                     f = n0 / (r * t) ln(2)
    #
    # So we need to know the initial number of atoms (n0), the capture rate (r),
    # and choose an irradiation time (t)

    w = model.geometry.get_materials_by_name('tungsten')[0]
    atom_densities = w.get_nuclide_atom_densities()
    atom_per_cc = 1e24 * atom_densities['W186'][1]  # Density in atom/cm^3
    n0 = atom_per_cc * w.volume  # Absolute number of atoms

    # Pick a random irradiation time and then determine necessary source rate to
    # reduce material by half
    t = uniform(1.0, 5.0) * 86400
    source_rates = [n0/(capture_rate*t) * log(2.0)]

    # Now activate the material
    integrator = openmc.deplete.PredictorIntegrator(
        op, [t], source_rates=source_rates
    )
    integrator.integrate()

    # Get resulting number of atoms
    results = openmc.deplete.ResultsList.from_hdf5('depletion_results.h5')
    _, atoms = results.get_atoms(str(w.id), "W186")

    assert atoms[0] == pytest.approx(n0)
    assert atoms[1] / atoms[0] == pytest.approx(0.5, rel=1e-3)
