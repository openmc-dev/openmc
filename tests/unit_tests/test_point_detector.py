"""Tests for the next-event (point detector) estimator.

These cover three properties that are easy to break and hard to notice:

* scoring a point detector must not disturb the transport random sequence,
  nor may one detector disturb another (the estimator draws random numbers
  to sample the outgoing energy, so it needs its own RNG substreams);
* the emitting particle's weight has to reach the score;
* the ray to the detector must accumulate optical depth over the distance to
  the detector, not out to the next surface behind it.
"""

import numpy as np
import pytest

import openmc
import openmc.data


DETECTOR = ((10.0, 0.0, 0.0), 1.0)


def _point_tally(detectors, scores=('flux',), energy_bins=None):
    """Build a tally filtered on the given point detectors."""
    filters = [openmc.PointFilter(list(detectors))]
    if energy_bins is not None:
        filters.append(openmc.EnergyFilter(energy_bins))
    tally = openmc.Tally(name='detector')
    tally.filters = filters
    tally.scores = list(scores)
    return tally


def _hydrogen_model(detectors, energy_bins=None, survival=False, absorber=False,
                    radius=100.0, particles=200, batches=10, seed=1):
    """Point source at the origin inside a homogeneous hydrogen sphere.

    Hydrogen is deliberate: awr < 1 puts elastic scattering on the
    double-valued CM->lab branch, which is the path that consumes random
    numbers inside the estimator. The optional B10 gives histories somewhere to
    terminate -- without it neutrons thermalize and wander for a very long time,
    since H1 capture alone barely removes anything.
    """
    h = openmc.Material()
    h.add_nuclide('H1', 1.0)
    if absorber:
        h.add_nuclide('B10', 0.02)
    h.set_density('atom/b-cm', 0.05)

    sphere = openmc.Sphere(r=radius, boundary_type='vacuum')
    cell = openmc.Cell(fill=h, region=-sphere)

    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.run_mode = 'fixed source'
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.seed = seed
    model.settings.survival_biasing = survival
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.delta_function(2.0e6),
    )
    if detectors is not None:
        model.tallies = openmc.Tallies(
            [_point_tally(detectors, energy_bins=energy_bins)])
    return model, h


def _fissile_model(detectors, particles=200, batches=10, inactive=2, seed=1):
    """Eigenvalue model exercising elastic, inelastic, fission and S(a,b)."""
    fuel = openmc.Material()
    fuel.add_nuclide('U235', 1.0)
    fuel.add_nuclide('O16', 2.0)
    fuel.set_density('g/cm3', 10.0)

    water = openmc.Material()
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    water.add_s_alpha_beta('c_H_in_H2O')

    inner = openmc.Sphere(r=8.0)
    outer = openmc.Sphere(r=30.0, boundary_type='vacuum')
    core = openmc.Cell(fill=fuel, region=-inner)
    reflector = openmc.Cell(fill=water, region=+inner & -outer)

    model = openmc.Model()
    model.geometry = openmc.Geometry([core, reflector])
    model.settings.particles = particles
    model.settings.batches = batches
    model.settings.inactive = inactive
    model.settings.seed = seed

    # A tracklength tally is the witness for "the transport did not change"
    witness = openmc.Tally(name='witness')
    witness.filters = [openmc.CellFilter([core, reflector])]
    witness.scores = ['flux', 'fission']
    witness.estimator = 'tracklength'

    tallies = [witness]
    if detectors is not None:
        tallies.append(_point_tally(detectors))
    model.tallies = openmc.Tallies(tallies)
    return model


def _run(model):
    """Run and return (keff, {tally name: (mean, std_dev)}).

    Pinned to one thread. Tally contributions are accumulated with atomics, so
    with several threads the summation order varies between runs and results
    agree only to the last couple of bits -- which would defeat the exact
    comparisons below. The properties under test here (whose random numbers a
    detector draws) are exactly reproducible; the floating-point summation
    order is not, and is not what these tests are about.

    Results are read eagerly. openmc.Tally loads them lazily from whichever
    statepoint it was linked to, and successive runs in one working directory
    overwrite that file, so anything not copied out here is gone by the time
    the next model has run.
    """
    sp_path = model.run(threads=1)
    results = {}
    with openmc.StatePoint(sp_path) as sp:
        keff = sp.keff if sp.run_mode == 'eigenvalue' else None
        for tally in sp.tallies.values():
            results[tally.name] = (tally.mean.copy(), tally.std_dev.copy())
    return keff, results


def test_point_detector_does_not_perturb_transport(run_in_tmpdir):
    """Adding a point tally must leave the transport bit-for-bit unchanged.

    The estimator samples an outgoing energy for every detector at every
    collision. If those draws come off STREAM_TRACKING, the particle histories
    depend on whether a point tally happens to be defined, which quietly
    invalidates any comparison against a reference run.
    """
    keff_ref, ref = _run(_fissile_model(None))
    keff_det, det = _run(_fissile_model([DETECTOR]))

    assert keff_det.nominal_value == keff_ref.nominal_value
    assert keff_det.std_dev == keff_ref.std_dev
    np.testing.assert_array_equal(det['witness'][0], ref['witness'][0])


def test_point_detector_independent_of_other_detectors(run_in_tmpdir):
    """A detector's score must not depend on how many others are present.

    Each detector is served from its own substream, offset deterministically
    from the per-event base seed, so the number of random numbers consumed by
    the preceding detectors is irrelevant.
    """
    second = ((0.0, 25.0, 0.0), 1.0)

    _, one = _run(_fissile_model([DETECTOR]))
    _, two = _run(_fissile_model([DETECTOR, second]))

    # Transport is untouched either way -- this part is exact
    np.testing.assert_array_equal(two['witness'][0], one['witness'][0])

    # ...and so is the first detector's bin. This is exact: every input to the
    # score -- the substream, the ray's own RNG state, the geometry -- is a
    # function of the detector's position alone.
    np.testing.assert_array_equal(
        two['detector'][0].reshape(2, -1)[0],
        one['detector'][0].reshape(1, -1)[0])


def test_point_detector_uncollided_attenuation(run_in_tmpdir):
    """Uncollided flux at the detector must match exp(-Sigma_t R)/(4 pi R^2).

    The source contribution is deterministic -- every source particle scores
    exactly this value -- so the comparison is essentially exact. It is also
    the quantity that breaks if the ray accumulates optical depth out to the
    next surface (here r = 100 cm) instead of stopping at the detector.
    """
    # Only the source term lands in a narrow window at the source energy:
    # elastic scattering on hydrogen always removes energy. The cutoff keeps
    # marginally-degraded neutrons from contributing a second time.
    energy_bins = [1.99e6, 2.01e6]
    model, material = _hydrogen_model(
        [DETECTOR], energy_bins=energy_bins, particles=100, batches=5)
    model.settings.cutoff = {'energy_neutron': 1.99e6}

    _, results = _run(model)
    mean, _ = results['detector']

    # Total macroscopic cross section at the source energy
    library = openmc.data.DataLibrary.from_xml()
    h1 = openmc.data.IncidentNeutron.from_hdf5(
        library.get_by_material('H1')['path'])
    sigma_t = h1[1].xs[h1.temperatures[0]](2.0e6) * 0.05

    distance = np.linalg.norm(DETECTOR[0])
    expected = np.exp(-sigma_t * distance) / (4.0 * np.pi * distance**2)

    assert mean.sum() == pytest.approx(expected, rel=0.02)


def test_point_detector_applies_particle_weight(run_in_tmpdir):
    """Survival biasing must not inflate the detector response.

    Under implicit capture the colliding particle carries a reduced weight;
    if the estimator scores at unit weight instead, the detector reads high by
    roughly the inverse of the surviving weight fraction. Analog and
    survival-biased runs are compared statistically, so the tolerance is loose
    -- the failure mode this guards against is order-unity, not marginal.
    """
    kwargs = dict(absorber=True, radius=30.0, particles=1000, batches=20)
    analog, _ = _hydrogen_model([DETECTOR], survival=False, **kwargs)
    biased, _ = _hydrogen_model([DETECTOR], survival=True, **kwargs)

    _, analog_results = _run(analog)
    _, biased_results = _run(biased)

    a_mean, a_std = analog_results['detector']
    b_mean, b_std = biased_results['detector']

    sigma = np.hypot(a_std.sum(), b_std.sum())
    difference = abs(b_mean.sum() - a_mean.sum())

    assert difference < max(4.0 * sigma, 0.05 * a_mean.sum())


# ---------------------------------------------------------------------------
# PointFilter API. These need no nuclear data.
# ---------------------------------------------------------------------------

def test_point_filter_xml_round_trip():
    """A PointFilter must survive a write/read cycle.

    The generic Filter.from_xml_element reads bins as integers, which cannot
    represent detector coordinates, so PointFilter has to override it. Without
    that, openmc.Tallies.from_xml -- and so Model.from_xml -- raises on any
    model containing a point detector.
    """
    detectors = [((10.0, 0.0, 0.0), 1.0), ((0.0, 25.0, -3.5), 2.0)]
    original = openmc.PointFilter(detectors)

    restored = openmc.Filter.from_xml_element(original.to_xml_element())

    assert isinstance(restored, openmc.PointFilter)
    assert restored.bins == original.bins
    assert restored.num_bins == len(detectors)


def test_point_filter_hdf5_round_trip(run_in_tmpdir):
    """The same for the statepoint representation.

    The group is built the way PointFilter::to_statepoint writes it: a flat
    'bins' dataset of four doubles per detector, plus 'n_bins'.
    """
    h5py = pytest.importorskip('h5py')

    detectors = [((1.5, -2.0, 3.0), 0.5), ((0.0, 0.0, 7.25), 0.0)]
    flat = [v for pos, r0 in detectors for v in (*pos, r0)]

    with h5py.File('filter.h5', 'w') as f:
        group = f.create_group('filter 1')
        group.create_dataset('type', data=np.bytes_('point'))
        group.create_dataset('n_bins', data=len(detectors))
        group.create_dataset('bins', data=np.array(flat))
    with h5py.File('filter.h5', 'r') as f:
        restored = openmc.PointFilter.from_hdf5(f['filter 1'])

    assert restored.bins == detectors
    assert restored.num_bins == len(detectors)


@pytest.mark.parametrize('bad_bins', [
    [((0.0, 0.0, 0.0), -1.0)],   # negative exclusion radius
    [((0.0, 0.0), 1.0)],         # position is not three-dimensional
    [((0.0, 0.0, 0.0),)],        # missing radius
])
def test_point_filter_rejects_invalid_bins(bad_bins):
    with pytest.raises(ValueError):
        openmc.PointFilter(bad_bins)


def test_point_filter_rejects_ragged_xml():
    """Bin values must come in groups of four (x, y, z, radius)."""
    import lxml.etree as ET

    elem = ET.fromstring(
        '<filter id="1" type="point"><bins>1.0 2.0 3.0</bins></filter>')
    with pytest.raises(ValueError, match='multiple of four'):
        openmc.PointFilter.from_xml_element(elem)


def test_point_filter_registered_in_lib():
    """openmc.lib must be able to construct a point filter by type name.

    openmc.lib._get_filter() looks the C++ filter type up in a table; a missing
    entry raises KeyError for any model that defines a point detector.
    """
    import openmc.lib
    assert openmc.lib.filter._FILTER_TYPE_MAP['point'] is openmc.lib.PointFilter


# ---------------------------------------------------------------------------
# Setup-time rejection of configurations the estimator cannot represent.
# These run the transport solver and so require nuclear data.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('boundary', ['reflective', 'white'])
def test_point_detector_rejects_nonvacuum_boundary(run_in_tmpdir, boundary):
    """Non-vacuum boundaries must be refused rather than silently mis-tallied.

    A next-event estimator draws a straight line from the emitting event to
    the detector, so it cannot account for particles that arrive after
    reflecting or being re-emitted by a white boundary. The match is on the
    specific message: a bare `raises(RuntimeError)` would also be satisfied by
    an unrelated failure such as missing nuclear data, and would keep passing
    if the guard were removed.
    """
    model, _ = _hydrogen_model([DETECTOR], particles=10, batches=2)
    for surface in model.geometry.get_all_surfaces().values():
        surface.boundary_type = boundary

    with pytest.raises(RuntimeError, match='non-vacuum boundary'):
        model.run()


def test_point_detector_rejects_photon_transport(run_in_tmpdir):
    """Photon transport must be refused for a tally that can see photons.

    The estimator has no scoring hooks in photon physics -- neither photon
    collisions nor secondary photon production contribute -- so a photon
    response would be missing every collided term while still looking
    plausible.
    """
    model, _ = _hydrogen_model([DETECTOR], particles=10, batches=2)
    model.settings.photon_transport = True

    with pytest.raises(RuntimeError, match='photon transport'):
        model.run()


def test_point_detector_allows_neutron_only_with_photon_transport(run_in_tmpdir):
    """...but a neutron-restricted tally stays legal.

    Enabling photon transport does not alter the neutron random walk, and
    every neutron emission path is hooked, so a detector filtered to neutrons
    is complete. The guard above must not be so broad that it blocks this.
    """
    model, _ = _hydrogen_model([DETECTOR], particles=10, batches=2)
    model.settings.photon_transport = True
    tally = model.tallies[0]
    tally.filters = tally.filters + [openmc.ParticleFilter(['neutron'])]

    model.run()  # must not raise


def test_point_detector_rejects_monodirectional_source(run_in_tmpdir):
    """A delta-function angular distribution has no density to evaluate.

    Sampling a position and then asking for the angular density toward the
    detector gives zero for almost every history, while the true uncollided
    flux is not zero: substituting r = D - s*u0 turns the contribution into
    the line integral of the spatial source density back along the beam, the
    1/distance^2 having cancelled against the volume element. Until that is
    implemented the run has to stop rather than drop the uncollided term.
    """
    model, _ = _hydrogen_model([DETECTOR], particles=10, batches=2)
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
        energy=openmc.stats.delta_function(1.0e6))

    with pytest.raises(RuntimeError, match='monodirectional'):
        model.run()


def test_point_detector_rejects_non_independent_source(run_in_tmpdir):
    """Only an independent source has an angular density to evaluate."""
    model, _ = _hydrogen_model([DETECTOR], particles=10, batches=2)

    # A real source file, so that reading settings.xml succeeds and the tally
    # check is actually reached
    openmc.write_source_file(
        [openmc.SourceParticle(r=(0., 0., 0.), E=1.0e6)], 'source.h5')
    model.settings.source = openmc.FileSource('source.h5')

    with pytest.raises(RuntimeError, match='independent source'):
        model.run()
