"""Multigroup photon transport with the Monte Carlo solver.

A reflective cube with a spatially uniform, isotropic source behaves as an
infinite medium, so the group track lengths (volume-integrated fluxes) per
source photon follow from a group balance and can be checked exactly.
The library uses P1 scattering and a top group that produces more photons
than it scatters, represented by a weight multiplicity that is either given
explicitly or derived by OpenMC when left out of the library. Besides the
tracklength reaction rates, an analog tally checks the group-to-group scatter
and nu-scatter matrices and their P1 moments: the scalar flux in an infinite
medium does not depend on the scattering anisotropy, so it cannot detect
errors in the angular sampling on its own.

"""
import numpy as np
import openmc

GROUP_EDGES = [1.0e3, 1.0e5, 5.0e5, 2.0e6]

# Macroscopic cross sections in 1/cm; index 0 is the highest-energy group
TOTAL = np.array([0.25, 0.40, 1.20])
ABSORPTION = np.array([0.05, 0.15, 1.00])

# P0 production (nu-scatter) matrix indexed [g_in, g_out]. The top group
# produces 0.25/cm but only scatters TOTAL - ABSORPTION = 0.20/cm; the extra
# photons stand in for fluorescence and annihilation photons.
PRODUCTION = np.array([
    [0.10, 0.08, 0.07],
    [0.00, 0.15, 0.10],
    [0.00, 0.00, 0.20],
])
MU_BAR = 0.2  # average scattering cosine, P1/P0

SPEED_OF_LIGHT = 2.99792458e10  # cm/s


def _multiplicity():
    """Row-constant multiplicity so each collision yields PRODUCTION/TOTAL"""
    ratio = PRODUCTION.sum(axis=1) / (TOTAL - ABSORPTION)
    return np.repeat(ratio[:, np.newaxis], len(TOTAL), axis=1)


def _exact_track_lengths():
    """Group track length per source photon in an infinite medium.

    Each group balances collisions against source and in-production:
    TOTAL[g] * L[g] = S[g] + sum_g' PRODUCTION[g', g] * L[g']

    """
    source = np.array([1.0, 0.0, 0.0])  # 1 MeV photons are born in group 1
    return np.linalg.solve(np.diag(TOTAL) - PRODUCTION.T, source)


def make_model(multiplicity='explicit', survival_biasing=False,
               absorption=ABSORPTION):
    """Build the photon model and write its multigroup library."""
    groups = openmc.mgxs.EnergyGroups(group_edges=GROUP_EDGES)
    photon = openmc.XSdata('photon', groups)
    photon.order = 1
    photon.set_total(TOTAL)
    photon.set_absorption(absorption)
    photon.set_scatter_matrix(
        np.stack([PRODUCTION, MU_BAR * PRODUCTION], axis=-1))
    if multiplicity == 'explicit':
        photon.set_multiplicity_matrix(_multiplicity())
    # Otherwise the multiplicity is derived when OpenMC reads the library.
    # No inverse velocity is given, so photons must default to 1/c

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
    if survival_biasing:
        model.settings.survival_biasing = True

    return model


def _assert_matches(name, mean, std_dev, reference):
    """Compare tally means to exact values within 5 standard deviations.

    Bins where the exact value is zero (upscatter) can never receive a
    score, so they must be exactly zero. Every other bin must be scored with
    a nonzero spread, so the deviation below is always finite.

    """
    if not (np.all(np.isfinite(mean)) and np.all(np.isfinite(std_dev))):
        raise AssertionError(
            f"Multigroup photon '{name}' has non-finite results: {mean} "
            f"+/- {std_dev}")
    zero = reference == 0.0
    if np.any(mean[zero] != 0.0):
        raise AssertionError(
            f"Multigroup photon '{name}' scored bins that must be empty: "
            f"{mean}")
    if np.any(std_dev[~zero] == 0.0):
        raise AssertionError(
            f"Multigroup photon '{name}' has no spread in bins that must be "
            f"scored: {mean} +/- {std_dev}")
    deviation = np.abs(mean[~zero] - reference[~zero]) / std_dev[~zero]
    if np.any(deviation > 5.0):
        raise AssertionError(
            f"Multigroup photon '{name}' disagrees with the infinite-medium "
            f"solution: {mean} vs. {reference} ({deviation} standard "
            "deviations)")


def check_exact_solution(statepoint_path):
    """Check every tally in the statepoint against the exact solution"""
    exact = _exact_track_lengths()

    with openmc.StatePoint(statepoint_path) as statepoint:
        tally = statepoint.get_tally(name='photon')
        expected = {
            'flux': exact,
            'total': TOTAL * exact,
            'absorption': ABSORPTION * exact,
            'scatter': (TOTAL - ABSORPTION) * exact,
            'nu-scatter': PRODUCTION.sum(axis=1) * exact,
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
        if not np.allclose(SPEED_OF_LIGHT * inverse_velocity, flux,
                           rtol=1e-10):
            raise AssertionError(
                'Multigroup photons are not moving at the speed of light: '
                f'flux / inverse-velocity = {flux / inverse_velocity}')

        # Group-to-group transfer. Each analog nu-scatter event scores the
        # multiplicity of the sampled transfer, and every transfer has the
        # mean cosine MU_BAR, so the P1 moments are MU_BAR times the P0 ones.
        nu_scatter = exact[:, np.newaxis] * PRODUCTION
        transfer_expected = {
            'scatter': nu_scatter / _multiplicity(),
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
