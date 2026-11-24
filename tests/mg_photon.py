"""Shared model and exact solutions for multigroup photon transport tests.

A reflective cube with a uniform, isotropic source represents an infinite
medium. The group balance gives the exact scalar flux. The P1 transfer tallies
also check the angular sampling, which has no effect on the scalar flux here.
"""
import numpy as np
import openmc

GROUP_EDGES = [1.0e3, 1.0e5, 5.0e5, 2.0e6]

# Macroscopic cross sections in 1/cm; index 0 is the highest-energy group
TOTAL = np.array([0.25, 0.40, 1.20])
ABSORPTION = np.array([0.05, 0.15, 1.00])
SCATTER_XS = TOTAL - ABSORPTION

# P0 nu-scatter matrix indexed [g_in, g_out]. Secondary photons raise the
# production cross section of the top group from 0.20/cm to 0.25/cm.
NU_SCATTER = np.array([
    [0.10, 0.08, 0.07],
    [0.00, 0.15, 0.10],
    [0.00, 0.00, 0.20],
])
MU_BAR = 0.2  # average scattering cosine, P1/P0
YIELD = NU_SCATTER.sum(axis=1) / SCATTER_XS
SPEED_OF_LIGHT = 2.99792458e10  # cm/s


def make_model(multiplicity='explicit', absorption=ABSORPTION):
    """Build the photon model and write its multigroup library."""
    groups = openmc.mgxs.EnergyGroups(group_edges=GROUP_EDGES)
    photon = openmc.XSdata('photon', groups)
    photon.order = 1
    photon.set_total(TOTAL)
    photon.set_absorption(absorption)
    photon.set_scatter_matrix(
        np.stack([NU_SCATTER, MU_BAR * NU_SCATTER], axis=-1))
    if multiplicity == 'explicit':
        photon.set_multiplicity_matrix(
            np.repeat(YIELD[:, None], len(TOTAL), axis=1))

    library = openmc.MGXSLibrary(groups, particle_type='photon')
    library.add_xsdata(photon)
    library.export_to_hdf5('mgxs.h5')

    material = openmc.Material()
    material.set_density('macro', 1.0)
    material.add_macroscopic('photon')

    box = openmc.model.RectangularParallelepiped(
        -5.0, 5.0, -5.0, 5.0, -5.0, 5.0, boundary_type='reflective')
    cell = openmc.Cell(fill=material, region=-box)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.materials = openmc.Materials([material])
    model.materials.cross_sections = 'mgxs.h5'

    energy_filter = openmc.EnergyFilter(GROUP_EDGES)
    particle_filter = openmc.ParticleFilter('photon')

    tally = openmc.Tally(name='photon')
    tally.filters = [energy_filter, particle_filter]
    tally.scores = ['flux', 'total', 'absorption', 'scatter', 'nu-scatter',
                    'inverse-velocity']

    # The outgoing energy and Legendre filters make this tally analog
    transfer = openmc.Tally(name='photon transfer')
    transfer.filters = [
        energy_filter,
        openmc.EnergyoutFilter(GROUP_EDGES),
        openmc.LegendreFilter(1),
        particle_filter,
    ]
    transfer.scores = ['scatter', 'nu-scatter']
    model.tallies = openmc.Tallies([tally, transfer])

    model.settings.energy_mode = 'multi-group'
    model.settings.run_mode = 'fixed source'
    model.settings.photon_transport = True
    model.settings.batches = 20
    model.settings.particles = 1000
    model.settings.source = openmc.IndependentSource(
        particle='photon',
        space=openmc.stats.Box((-5.0, -5.0, -5.0), (5.0, 5.0, 5.0)),
        energy=openmc.stats.Discrete([1.0e6], [1.0]))
    return model


def _assert_matches(name, mean, std_dev, reference):
    """Check scored bins within five sigma and empty bins exactly."""
    assert np.isfinite(mean).all() and np.isfinite(std_dev).all(), (
        f'{name} has nonfinite tally results')
    zero = reference == 0.0
    np.testing.assert_array_equal(mean[zero], 0.0, err_msg=name)
    assert np.all(std_dev[~zero] > 0.0), f'{name} needs positive uncertainty'
    deviation = np.abs(mean[~zero] - reference[~zero]) / std_dev[~zero]
    assert np.all(deviation <= 5.0), f'{name} deviations in sigma: {deviation}'


def check_exact_solution(statepoint_path):
    """Check every tally in the statepoint against the exact solution"""
    # Collisions balance the source in the top group and production from
    # scattering in all groups.
    exact = np.linalg.solve(np.diag(TOTAL) - NU_SCATTER.T, [1.0, 0.0, 0.0])

    with openmc.StatePoint(statepoint_path) as statepoint:
        tally = statepoint.get_tally(name='photon')
        expected = {
            'flux': exact,
            'total': TOTAL * exact,
            'absorption': ABSORPTION * exact,
            'scatter': SCATTER_XS * exact,
            'nu-scatter': NU_SCATTER.sum(axis=1) * exact,
        }

        # Energy filter bins run from low to high energy, so reverse them
        # to match the group ordering of the library
        for score, reference in expected.items():
            mean = tally.get_values(scores=[score]).ravel()[::-1]
            std_dev = tally.get_values(
                scores=[score], value='std_dev').ravel()[::-1]
            _assert_matches(score, mean, std_dev, reference)

        # Photons travel at the speed of light, so the inverse-velocity
        # score must equal flux / c in every group
        flux = tally.get_values(scores=['flux']).ravel()
        inverse_velocity = tally.get_values(
            scores=['inverse-velocity']).ravel()
        np.testing.assert_allclose(
            SPEED_OF_LIGHT * inverse_velocity, flux, rtol=1e-10, atol=1e-8,
            equal_nan=False, err_msg='Photon inverse velocity must be 1/c')

        # Every transfer has mean cosine MU_BAR, so P1 is MU_BAR times P0.
        nu_scatter = exact[:, None] * NU_SCATTER
        transfer_expected = {
            'scatter': nu_scatter / YIELD[:, None],
            'nu-scatter': nu_scatter,
        }
        transfer = statepoint.get_tally(name='photon transfer')
        # Axes: energy in, energy out, Legendre order, particle, nuclide, score
        means = transfer.get_reshaped_data(value='mean')
        std_devs = transfer.get_reshaped_data(value='std_dev')
        for score, p0 in transfer_expected.items():
            i = transfer.get_score_index(score)
            for order, reference in enumerate([p0, MU_BAR * p0]):
                _assert_matches(
                    f'{score} P{order} matrix',
                    means[::-1, ::-1, order, 0, 0, i],
                    std_devs[::-1, ::-1, order, 0, 0, i],
                    reference)
