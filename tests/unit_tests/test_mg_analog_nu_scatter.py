"""Analog nu-scatter tallies in multigroup mode.

Multigroup Monte Carlo represents the extra particles produced by scattering
reactions such as (n,2n) through the weight multiplicity of each sampled
group-to-group transfer, so an analog nu-scatter tally has to score the
pre-collision weight times that multiplicity. In an infinite medium the
expected group-to-group scores follow exactly from a group balance.

"""
import numpy as np
import openmc
import pytest

GROUP_EDGES = [1.0e-5, 1.0, 1.0e3, 2.0e7]

# Macroscopic data in 1/cm, indexed [g_in, g_out] with index 0 the
# highest-energy group. The multiplicity depends on the incoming group only:
# multigroup Monte Carlo samples the outgoing group from the nu-scatter
# matrix before applying the multiplicity, which reproduces the library data
# only when the multiplicity is the same for every outgoing group.
NU_SCATTER = np.array([
    [0.10, 0.12, 0.06],
    [0.00, 0.20, 0.15],
    [0.00, 0.00, 0.30],
])
MULTIPLICITY = np.repeat([[1.4], [1.25], [1.0]], 3, axis=1)
SCATTER = NU_SCATTER / MULTIPLICITY
ABSORPTION = np.array([0.05, 0.10, 0.50])
TOTAL = ABSORPTION + SCATTER.sum(axis=1)


def _model(survival_biasing):
    groups = openmc.mgxs.EnergyGroups(GROUP_EDGES)
    xsdata = openmc.XSdata('mat', groups)
    xsdata.order = 0
    xsdata.set_total(TOTAL)
    xsdata.set_absorption(ABSORPTION)
    xsdata.set_scatter_matrix(NU_SCATTER[..., np.newaxis])
    xsdata.set_multiplicity_matrix(MULTIPLICITY)
    library = openmc.MGXSLibrary(groups)
    library.add_xsdata(xsdata)
    library.export_to_hdf5('mgxs.h5')

    material = openmc.Material()
    material.set_density('macro', 1.0)
    material.add_macroscopic('mat')

    # A reflective box with a uniform, isotropic source is an infinite medium
    box = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
    model = openmc.Model()
    model.geometry = openmc.Geometry([openmc.Cell(fill=material, region=-box)])
    model.materials = openmc.Materials([material])
    model.materials.cross_sections = 'mgxs.h5'

    # The outgoing energy filter makes nu-scatter an analog tally
    tally = openmc.Tally()
    tally.filters = [
        openmc.EnergyFilter(GROUP_EDGES),
        openmc.EnergyoutFilter(GROUP_EDGES),
    ]
    tally.scores = ['scatter', 'nu-scatter']
    model.tallies = [tally]

    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 20
    model.settings.particles = 1000
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))
    if survival_biasing:
        model.settings.survival_biasing = True
    return model


@pytest.mark.parametrize('survival_biasing', [False, True])
def test_analog_nu_scatter_multiplicity(run_in_tmpdir, survival_biasing):
    model = _model(survival_biasing)
    model.run(apply_tally_results=True)
    tally = model.tallies[0]

    # Group track lengths per source neutron born in the top group:
    # TOTAL[g] * L[g] = S[g] + sum_g' NU_SCATTER[g', g] * L[g']
    track_length = np.linalg.solve(
        np.diag(TOTAL) - NU_SCATTER.T, [1.0, 0.0, 0.0])
    expected = {
        'scatter': track_length[:, np.newaxis] * SCATTER,
        'nu-scatter': track_length[:, np.newaxis] * NU_SCATTER,
    }

    # Axes: energy in, energy out, nuclide, score. Energy filter bins run
    # from low to high energy, so reverse them to get library group order.
    means = tally.get_reshaped_data(value='mean')
    std_devs = tally.get_reshaped_data(value='std_dev')
    for score, reference in expected.items():
        i = tally.get_score_index(score)
        mean = means[::-1, ::-1, 0, i]
        std_dev = std_devs[::-1, ::-1, 0, i]
        possible = reference > 0.0
        assert np.all(mean[~possible] == 0.0), score
        deviation = np.abs(mean - reference)[possible] / std_dev[possible]
        assert np.all(deviation < 5.0), \
            f'{score}: {mean} vs. {reference} ({deviation} std. dev.)'