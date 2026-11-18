from math import pi, log, log10
from random import uniform, normalvariate
from types import SimpleNamespace
from unittest import mock

import numpy as np

import openmc.deplete
import openmc.deplete.microxs as microxs_mod
import openmc
import openmc.lib
import pytest
from openmc.data import REACTION_MT
from openmc.deplete.helpers import FluxCollapseHelper
from openmc.mgxs import GROUP_STRUCTURES


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


ENERGIES = np.logspace(log10(1e-5), log10(2e7), 100)


@pytest.mark.parametrize("reaction_rate_mode,reaction_rate_opts,tolerance", [
    ("direct", {}, 1e-5),
    ("flux", {'energies': ENERGIES}, 0.1),
    ("flux", {'energies': ENERGIES, 'reactions': ['(n,gamma)']}, 1e-5),
    ("flux", {'energies': ENERGIES, 'reactions': ['(n,gamma)'], 'nuclides': ['W186', 'H3']}, 1e-2),
])
@pytest.mark.flaky(reruns=1)
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

    # Create transport operator
    op = openmc.deplete.CoupledOperator(
        model, chain,
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

    # Check that material name is preserved in depletion results
    step_result = results[0]
    mat_from_results = step_result.get_material(f"{w.id}")
    assert mat_from_results.name == 'tungsten'


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


CASMO40 = np.asarray(GROUP_STRUCTURES['CASMO-40'], dtype=float)
N_GROUPS = CASMO40.size - 1
TEMPERATURE = 293.6


def _wire_helper(nuclides, scores, flux, n_nucs=None,
                 reactions_direct=None, nuclides_direct=None):
    """Wire a FluxCollapseHelper for one material without building tallies."""
    helper = FluxCollapseHelper(
        n_nucs or len(nuclides), len(scores), CASMO40,
        reactions=reactions_direct, nuclides=nuclides_direct)
    helper._materials = [SimpleNamespace(temperature=TEMPERATURE)]
    helper._xs_tables = {}
    helper._mts = [REACTION_MT[s] for s in scores]
    helper._scores = list(scores)
    helper._flux_tally_means_cache = np.asarray(flux, dtype=float)
    helper.nuclides = list(nuclides)
    return helper


def test_flux_collapse_helper_cache_lifecycle(monkeypatch):
    """The FluxCollapseHelper per-temperature cache lifecycle: flux-path rates
    equal collapse_rate; the table is built once per temperature; re-setting the
    same nuclides (as happens each timestep) does not rebuild it; and growing the
    nuclide set clears the cache, so the table is rebuilt and the newly-added
    nuclide gets correct (non-stale) rates.
    """
    nuclides = ['U235', 'U238', 'O16']
    grown = nuclides + ['Pu239']  # larger set, added later
    scores = ['fission', '(n,gamma)', '(n,2n)']  # (n,2n) is a threshold reaction
    flux = np.random.default_rng(0).random(N_GROUPS)
    react_index = list(range(len(scores)))

    build_spy = mock.MagicMock(wraps=microxs_mod._build_xs_table_ce)
    monkeypatch.setattr(microxs_mod, '_build_xs_table_ce', build_spy)

    with openmc.lib.TemporarySession():
        # Size the result cache for the larger (grown) set from the start
        helper = _wire_helper(nuclides, scores, flux, n_nucs=len(grown))
        rates = helper.get_material_rates(
            0, list(range(len(nuclides))), react_index).copy()

        # Flux-path rate equals collapse_rate for every (nuclide, score) ...
        for i, name in enumerate(nuclides):
            nuc = openmc.lib.nuclides[name]
            for j, s in enumerate(scores):
                expected = nuc.collapse_rate(
                    REACTION_MT[s], TEMPERATURE, CASMO40, flux)
                assert rates[i, j] == pytest.approx(expected, rel=1e-10)

        # ... built once for the single temperature ...
        assert build_spy.call_count == 1

        # ... and re-setting the same nuclide list (each step) does not rebuild
        helper.nuclides = list(nuclides)
        helper.get_material_rates(0, list(range(len(nuclides))), react_index)
        assert build_spy.call_count == 1

        # Growing the set clears the cache -> exactly one rebuild at the same T
        helper.nuclides = list(grown)
        rates = helper.get_material_rates(
            0, list(range(len(grown))), react_index).copy()
        assert build_spy.call_count == 2

        # The newly-added nuclide gets correct rates, not stale zeros
        pu = openmc.lib.nuclides['Pu239']
        for j, s in enumerate(scores):
            expected = pu.collapse_rate(
                REACTION_MT[s], TEMPERATURE, CASMO40, flux)
            assert rates[len(nuclides), j] == pytest.approx(expected, rel=1e-10)
        assert rates[len(nuclides), 0] > 0.0  # Pu239 fission nonzero -> not stale


def test_flux_collapse_helper_direct_override():
    """A direct-tally (nuclide, reaction) pair overrides the flux-collapsed value."""
    nuclides = ['U235', 'U238']
    scores = ['fission', '(n,gamma)']
    flux = np.random.default_rng(1).random(N_GROUPS)

    with openmc.lib.TemporarySession():
        # Direct-tally U235 fission; everything else via flux collapse
        helper = _wire_helper(nuclides, scores, flux,
            reactions_direct=['fission'], nuclides_direct=['U235'])
        direct_value = 1.2345e-3
        helper._rate_tally = SimpleNamespace(nuclides=['U235'])
        helper._rate_tally_means_cache = np.array([[direct_value]])

        rates = helper.get_material_rates(
            0, list(range(len(nuclides))), list(range(len(scores)))).copy()

        # U235 fission comes from the direct tally, not the flux collapse
        assert rates[0, 0] == pytest.approx(direct_value)

        # U235 (n,gamma) still comes from the flux collapse
        nuc = openmc.lib.nuclides['U235']
        expected = nuc.collapse_rate(
            REACTION_MT['(n,gamma)'], TEMPERATURE, CASMO40, flux)
        assert rates[0, 1] == pytest.approx(expected, rel=1e-10)
