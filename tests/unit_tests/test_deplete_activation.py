from math import pi, log, log10
from random import uniform, normalvariate

import numpy as np

import openmc.deplete
import openmc
import pytest


@pytest.fixture
def model():
    """Sphere of single nuclide"""
    model = openmc.Model()

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
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0])
    )
    model.settings.run_mode = 'fixed source'

    rx_tally = openmc.Tally(name='activation tally')
    rx_tally.scores = ['(n,gamma)']
    model.tallies.append(rx_tally)

    return model


ENERGIES = np.logspace(log10(1e-5), log10(150e6), 500)


@pytest.mark.parametrize("reaction_rate_mode,reaction_rate_opts,tolerance", [
    ("direct", {}, 1e-5),
    ("flux", {'energies': ENERGIES}, 0.01),
    ("flux", {'energies': ENERGIES, 'reactions': ['(n,gamma)']}, 1e-5),
    ("flux", {'energies': ENERGIES, 'reactions': ['(n,gamma)'], 'nuclides': ['W186', 'H3']}, 1e-2),
])
def test_activation(run_in_tmpdir, model, reaction_rate_mode, reaction_rate_opts, tolerance):
    # Determine (n.gamma) reaction rate using initial run
    sp = model.run()
    with openmc.StatePoint(sp) as sp:
        tally = sp.get_tally(name='activation tally')
        capture_rate = tally.mean.flat[0]

    # Create one-nuclide depletion chain
    chain = openmc.deplete.Chain()
    w186 = openmc.deplete.Nuclide('W186')
    w186.add_reaction('(n,gamma)', None, 0.0, 1.0)
    chain.add_nuclide(w186)
    chain.export_to_xml('test_chain.xml')

    # Create transport operator
    op = openmc.deplete.CoupledOperator(
        model, 'test_chain.xml',
        normalization_mode="source-rate",
        reaction_rate_mode=reaction_rate_mode,
        reaction_rate_opts=reaction_rate_opts,
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
    atom_per_cc = 1e24 * atom_densities['W186']  # Density in atom/cm^3
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
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(w, "W186")

    assert atoms[0] == pytest.approx(n0)
    assert atoms[1] / atoms[0] == pytest.approx(0.5, rel=tolerance)


def test_decay(run_in_tmpdir):
    """Test decay-only timesteps where no transport solve is performed"""

    # Create a model with a single nuclide, Sr89
    mat = openmc.Material()
    mat.add_nuclide('Sr89', 1.0)
    mat.set_density('g/cm3', 1.0)
    mat.depletable = True
    r = 5.0
    mat.volume = 4/3 * pi * r**3
    surf = openmc.Sphere(r=r, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-surf)
    geometry = openmc.Geometry([cell])
    settings = openmc.Settings()
    settings.batches = 10
    settings.particles = 1000
    settings.run_mode = 'fixed source'

    # Create depletion chain with only Sr89 and sample its half-life. Note that
    # currently at least one reaction has to exist in the depletion chain
    chain = openmc.deplete.Chain()
    sr89 = openmc.deplete.Nuclide('Sr89')
    sr89.half_life = normalvariate(4365792.0, 6048.0)
    sr89.add_decay_mode('beta-', None, 1.0)
    sr89.add_reaction('(n,gamma)', None, 0.0, 1.0)
    chain.add_nuclide(sr89)
    chain.export_to_xml('test_chain.xml')

    model = openmc.Model(geometry=geometry, settings=settings)
    # Create transport operator
    op = openmc.deplete.CoupledOperator(
        model, 'test_chain.xml', normalization_mode="source-rate"
    )

    # Deplete with two decay steps
    integrator = openmc.deplete.PredictorIntegrator(
        op, [sr89.half_life, 2*sr89.half_life], source_rates=[0.0, 0.0]
    )
    integrator.integrate()

    # Get resulting number of atoms
    results = openmc.deplete.Results('depletion_results.h5')
    _, atoms = results.get_atoms(mat, "Sr89")

    # Ensure density goes down by a factor of 2 after each half-life
    assert atoms[1] / atoms[0] == pytest.approx(0.5)
    assert atoms[2] / atoms[1] == pytest.approx(0.25)


def test_flux_rr_missing_nuclide(run_in_tmpdir, model):
    # Create two-nuclide depletion chain -- since W184 is not in the model, this
    # test ensures that FluxCollapseHelper loads missing nuclides appropriately
    chain = openmc.deplete.Chain()
    w184 = openmc.deplete.Nuclide('W184')
    w184.add_reaction('(n,gamma)', None, 0.0, 1.0)
    chain.add_nuclide(w184)
    w186 = openmc.deplete.Nuclide('W186')
    w186.add_reaction('(n,gamma)', None, 0.0, 1.0)
    chain.add_nuclide(w186)
    chain.export_to_xml('test_chain.xml')

    # Create transport operator
    op = openmc.deplete.CoupledOperator(
        model, 'test_chain.xml',
        normalization_mode="source-rate",
        reaction_rate_mode="flux",
        reaction_rate_opts={'energies': [0.0, 20.0e6]},
    )

    # Deplete with two decay steps
    integrator = openmc.deplete.PredictorIntegrator(
        op, [100.0], source_rates=[10.0]
    )
    integrator.integrate()
