"""Multigroup reaction-rate tallies must honor cell density multipliers.

In MG mode, transport uses macroscopic cross sections scaled by the cell
density multiplier. Every macroscopic reaction-rate score in a tally must be
evaluated at that same density, so ratios between scores on the same tally
(e.g. scatter/total) are independent of the multiplier and equal to the ratios
of the input cross sections, history by history.
"""
import openmc
import pytest


TOTAL = 1.0
ABSORPTION = 0.5
SCATTER = 0.5
FISSION = 0.2
NU_FISSION = 0.25
KAPPA_FISSION = 200.0e6 * FISSION


def _model(density_mult, estimator):
    groups = openmc.mgxs.EnergyGroups([0.0, 20.0e6])
    xs = openmc.XSdata('mat', groups, representation='isotropic')
    xs.order = 0
    xs.set_total([TOTAL])
    xs.set_absorption([ABSORPTION])
    xs.set_scatter_matrix([[[SCATTER]]])
    xs.set_fission([FISSION])
    xs.set_nu_fission([NU_FISSION])
    xs.set_kappa_fission([KAPPA_FISSION])
    xs.set_chi([1.0])
    library = openmc.MGXSLibrary(groups)
    library.add_xsdata(xs)
    library.export_to_hdf5('mgxs.h5')

    mat = openmc.Material()
    mat.add_macroscopic('mat')
    mat.set_density('macro', 1.0)
    materials = openmc.Materials([mat])
    materials.cross_sections = 'mgxs.h5'

    sphere = openmc.Sphere(r=5.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    cell.density = density_mult

    settings = openmc.Settings()
    settings.energy_mode = 'multi-group'
    settings.run_mode = 'fixed source'
    settings.batches = 3
    settings.particles = 200
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([1.0e6], [1.0]),
    )

    tally = openmc.Tally()
    tally.estimator = estimator
    tally.scores = ['flux', 'total', 'absorption', 'scatter', 'nu-scatter',
                    'fission', 'nu-fission', 'kappa-fission']
    return openmc.Model(openmc.Geometry([cell]), materials, settings,
                        openmc.Tallies([tally])), tally.id


@pytest.mark.parametrize('estimator', ['tracklength', 'collision'])
@pytest.mark.parametrize('density_mult', [0.5, 1.0, 2.0])
def test_mg_scores_scale_with_cell_density(run_in_tmpdir, density_mult,
                                           estimator):
    model, tally_id = _model(density_mult, estimator)
    sp_path = model.run()
    with openmc.StatePoint(sp_path) as sp:
        tally = sp.tallies[tally_id]
        mean = {s: tally.get_values(scores=[s]).sum() for s in tally.scores}

    assert mean['total'] > 0.0

    # Closure on the same particle histories
    assert mean['absorption'] + mean['scatter'] == pytest.approx(
        mean['total'], rel=1e-10)

    # Every macroscopic reaction rate is consistent with the total
    for score, xs in [('absorption', ABSORPTION), ('scatter', SCATTER),
                      ('nu-scatter', SCATTER), ('fission', FISSION),
                      ('nu-fission', NU_FISSION),
                      ('kappa-fission', KAPPA_FISSION)]:
        assert mean[score] / mean['total'] == pytest.approx(
            xs / TOTAL, rel=1e-10), score

    # The total rate is the flux times the density-scaled total cross section
    if estimator == 'tracklength':
        assert mean['total'] / mean['flux'] == pytest.approx(
            density_mult * TOTAL, rel=1e-10)
