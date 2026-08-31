import numpy as np
import pytest
import openmc


@pytest.fixture
def model():
    openmc.reset_auto_ids()
    model = openmc.Model()

    # Define materials
    NaI = openmc.Material()
    NaI.set_density('g/cm3', 3.7)
    NaI.add_element('Na', 1.0)
    NaI.add_element('I', 1.0)

    # Define geometry: two spheres in each other
    s1 = openmc.Sphere(r=1)
    s2 = openmc.Sphere(r=2, boundary_type='vacuum')
    inner_sphere = openmc.Cell(name='inner sphere', fill=NaI, region=-s1)
    outer_sphere = openmc.Cell(name='outer sphere', fill=NaI, region=+s1 & -s2)
    model.geometry = openmc.Geometry([inner_sphere, outer_sphere])

    # Define settings
    model.settings.run_mode = 'fixed source'
    model.settings.batches = 1
    model.settings.particles = 10000
    model.settings.photon_transport = True
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1e6),
        particle='photon'
    )

    # Define tallies
    energy_filter = openmc.EnergyFilter([1e3, 1e7])    

    tally1 = openmc.Tally()
    tally1.scores = ['pulse-height']
    cell_filter1 = openmc.CellFilter([inner_sphere, outer_sphere])
    tally1.filters = [cell_filter1, energy_filter]

    tally2 = openmc.Tally()
    tally2.scores = ['pulse-height']
    cell_filter2 = openmc.CellFilter([outer_sphere, inner_sphere])
    tally2.filters = [cell_filter2, energy_filter]

    model.tallies = [tally1, tally2]
    return model


def test_pulse_height(model, run_in_tmpdir):
    sp_path = model.run()
    sp = openmc.StatePoint(sp_path)
    t1 = sp.tallies[1].mean.squeeze()
    t2 = sp.tallies[2].mean.squeeze()
    
    np.testing.assert_array_equal(t1, t2[::-1])


# ---------------------------------------------------------------------------
# Shared secondary bank
#
# A pulse-height tally scores exactly one count per source history, in the bin
# containing the total energy that history's entire particle tree deposited in a
# cell. In the default transport modes the whole tree lives in one Particle
# object, so that total is available directly at particle death. Under the
# shared secondary bank each secondary generation is transported as a fresh set
# of Particle objects, redistributed across MPI ranks between generations, so a
# history's deposition is spread over many Particle objects and has to be
# reassembled before scoring.
#
# The comparison between modes cannot be exact. compute_particle_id() and
# compute_transport_seed() both take a different branch when the shared bank is
# active, so the two modes sample different random number streams and produce
# different realizations of the same distribution. Only count conservation is
# exact; everything else is compared statistically.
#
# These run single-rank. The cross-rank aggregation in
# finalize_pulse_height_tallies() is only exercised under MPI.
#
# One trap when adding tests here: never place an energy filter edge at a
# deposition value a history can hit exactly, such as the source energy of a
# monoenergetic source in a fully absorbing detector. Floating point in the
# accumulated sum then decides which side of the edge each history falls on,
# roughly 6% land above it and are dropped from the tally, and what looks like
# a physics discrepancy is only a difference in rounding.
# ---------------------------------------------------------------------------


def _detector_model(particle, radius, energy_bounds, particles=1000,
                    batches=10, shared_secondary=False):
    """NaI sphere in a void, with a pulse-height tally on the detector cell."""
    openmc.reset_auto_ids()
    model = openmc.Model()

    NaI = openmc.Material()
    NaI.set_density('g/cm3', 3.7)
    NaI.add_element('Na', 1.0)
    NaI.add_element('I', 1.0)

    detector_surf = openmc.Sphere(r=radius)
    outer_surf = openmc.Sphere(r=radius + 1.0, boundary_type='vacuum')
    detector = openmc.Cell(name='detector', fill=NaI, region=-detector_surf)
    outside = openmc.Cell(name='outside', region=+detector_surf & -outer_surf)
    model.geometry = openmc.Geometry([detector, outside])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = batches
    model.settings.particles = particles
    model.settings.photon_transport = True
    model.settings.shared_secondary_bank = shared_secondary
    model.settings.source = openmc.IndependentSource(
        energy=openmc.stats.delta_function(1e6),
        particle=particle
    )

    tally = openmc.Tally(name='pht')
    tally.scores = ['pulse-height']
    tally.filters = [
        openmc.CellFilter(detector),
        openmc.EnergyFilter(energy_bounds),
    ]
    model.tallies = openmc.Tallies([tally])

    return model


def _pulse_height(statepoint_path):
    """Return the pulse-height spectrum and its per-bin standard deviation."""
    with openmc.StatePoint(statepoint_path) as sp:
        tally = sp.get_tally(name='pht')
        return tally.mean.ravel().copy(), tally.std_dev.ravel().copy()


def _assert_spectra_consistent(a, a_err, b, b_err, n_sigma=5.0):
    """Compare two spectra bin by bin against their combined standard errors.

    Only bins holding at least 1% of histories are compared. Bins in the tail
    carry few counts per batch, so their batch-to-batch spread is a poor
    estimate of their true uncertainty and would drive spurious failures.
    """
    sigma = np.hypot(a_err, b_err)
    significant = (0.5 * (a + b) > 0.01) & (sigma > 0.0)
    assert significant.any(), "no bins carry enough counts to compare"
    z = np.abs(a[significant] - b[significant]) / sigma[significant]
    assert z.max() < n_sigma, f"largest per-bin discrepancy is {z.max():.1f} sigma"


def _mean_deposition(spectrum, std_dev, bounds):
    """Mean deposited energy per history, and its standard error."""
    centers = 0.5 * (bounds[:-1] + bounds[1:])
    mean = float(centers @ spectrum)
    err = float(np.sqrt(np.sum((centers * std_dev) ** 2)))
    return mean, err


@pytest.mark.parametrize('shared_secondary', [False, True])
@pytest.mark.parametrize('particle', ['photon', 'neutron'])
def test_pulse_height_count_conservation(particle, shared_secondary,
                                         run_in_tmpdir):
    """Every history scores exactly once, in both transport modes.

    Fixed-source results are normalized by the source strength divided by the
    number of source particles (Tally::accumulate), so summing a pulse-height
    tally over all its energy bins gives the number of scores per source
    particle, which must be exactly one. The energy filter is wide enough that
    no history can deposit outside it.

    This is the direct test for scoring per track rather than per history: if
    each secondary were scored separately the sum would exceed one by the mean
    number of tracks per history. A history whose deposition is dropped instead
    makes the sum fall short.
    """
    # Neutron capture in iodine releases several MeV of prompt gammas, so
    # deposition is not bounded by the source energy in the neutron case.
    upper = 1.1e6 if particle == 'photon' else 20.0e6
    model = _detector_model(
        particle, radius=1.0, energy_bounds=np.linspace(0.0, upper, 101),
        shared_secondary=shared_secondary,
    )

    spectrum, _ = _pulse_height(model.run())

    assert spectrum.sum() == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize('particle', ['photon', 'neutron'])
def test_shared_secondary_matches_local(particle, run_in_tmpdir):
    """Both modes give the same pulse-height distribution in a thin detector.

    A 1 cm NaI sphere at 1 MeV produces short cascades, so this mostly exercises
    the aggregation bookkeeping rather than deep secondary trees. Compared are
    the mean deposited energy per history and the spectrum shape bin by bin.
    """
    upper = 1.1e6 if particle == 'photon' else 20.0e6
    bounds = np.linspace(0.0, upper, 51)

    local = _detector_model(particle, radius=1.0, energy_bounds=bounds,
                            particles=2000, shared_secondary=False)
    local_spec, local_err = _pulse_height(local.run())

    shared = _detector_model(particle, radius=1.0, energy_bounds=bounds,
                             particles=2000, shared_secondary=True)
    shared_spec, shared_err = _pulse_height(shared.run())

    # Count conservation must hold in both before the shapes are compared
    assert local_spec.sum() == pytest.approx(1.0, abs=1e-9)
    assert shared_spec.sum() == pytest.approx(1.0, abs=1e-9)

    local_mean, local_mean_err = _mean_deposition(local_spec, local_err, bounds)
    shared_mean, shared_mean_err = _mean_deposition(
        shared_spec, shared_err, bounds)
    mean_sigma = np.hypot(local_mean_err, shared_mean_err)
    assert mean_sigma > 0.0
    assert abs(local_mean - shared_mean) < 4.0 * mean_sigma

    _assert_spectra_consistent(local_spec, local_err, shared_spec, shared_err)


def test_shared_secondary_matches_local_thick_detector(run_in_tmpdir):
    """Both modes agree when the secondary cascade is deep.

    A 50 cm NaI sphere is tens of mean free paths thick at 1 MeV, so nothing
    escapes and each history spawns roughly 1.7 secondaries through Compton
    scattering, fluorescence and bremsstrahlung. That makes this far more
    sensitive than the thin-detector case to the subtraction performed on the
    parent in create_secondary() being matched by the descendants' own
    contributions once they are reassembled.

    The top bin deliberately extends past the source energy. Pulse height
    accumulates as a telescoping sum of E_last() - E() over a history, which
    equals the source energy in exact arithmetic but not in floating point:
    measured on this geometry, about 94% of histories sum to exactly 1 MeV or
    just below and about 6% land a few ULPs above it. EnergyFilter matches on
    E >= bins.front() && E <= bins.back(), so a top edge sitting exactly at
    1 MeV drops that 6% entirely. They are then absent from the tally, count
    conservation silently breaks, and the comparison below degenerates into a
    comparison of rounding behaviour between two different random number
    streams rather than of spectra.
    """
    # Top bin straddles the full-energy peak so histories whose floating point
    # sum lands marginally above 1 MeV are still binned
    bounds = np.concatenate([np.linspace(0.0, 0.99e6, 20), [1.1e6]])

    local = _detector_model('photon', radius=50.0, energy_bounds=bounds,
                            shared_secondary=False)
    local_spec, local_err = _pulse_height(local.run())

    shared = _detector_model('photon', radius=50.0, energy_bounds=bounds,
                             shared_secondary=True)
    shared_spec, shared_err = _pulse_height(shared.run())

    assert local_spec.sum() == pytest.approx(1.0, abs=1e-9)
    assert shared_spec.sum() == pytest.approx(1.0, abs=1e-9)

    # The full-energy peak holds the great majority of histories, so its
    # fraction is the sharpest single statistic available here.
    peak_sigma = np.hypot(local_err[-1], shared_err[-1])
    assert peak_sigma > 0.0
    assert abs(local_spec[-1] - shared_spec[-1]) < 4.0 * peak_sigma

    _assert_spectra_consistent(local_spec, local_err, shared_spec, shared_err)
