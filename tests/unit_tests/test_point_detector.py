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


def _hydrogen_model(detectors, energy_bins=None, survival=False,
                    particles=200, batches=10, seed=1):
    """Point source at the origin inside a homogeneous hydrogen sphere.

    Hydrogen is deliberate: awr < 1 puts elastic scattering on the
    double-valued CM->lab branch, which is the path that consumes random
    numbers inside the estimator.
    """
    h = openmc.Material()
    h.add_nuclide('H1', 1.0)
    h.set_density('atom/b-cm', 0.05)

    sphere = openmc.Sphere(r=100.0, boundary_type='vacuum')
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
    """Run and return (keff, {tally name: mean array})."""
    sp_path = model.run()
    results = {}
    with openmc.StatePoint(sp_path) as sp:
        keff = sp.keff if sp.run_mode == 'eigenvalue' else None
        for tally in sp.tallies.values():
            results[tally.name] = tally.mean.copy()
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
    np.testing.assert_array_equal(det['witness'], ref['witness'])


def test_point_detector_independent_of_other_detectors(run_in_tmpdir):
    """A detector's score must not depend on how many others are present.

    Each detector is served from its own substream, offset deterministically
    from the per-event base seed, so the number of random numbers consumed by
    the preceding detectors is irrelevant.
    """
    second = ((0.0, 25.0, 0.0), 1.0)

    _, one = _run(_fissile_model([DETECTOR]))
    _, two = _run(_fissile_model([DETECTOR, second]))

    # Transport is untouched either way
    np.testing.assert_array_equal(two['witness'], one['witness'])

    # ...and so is the first detector's bin
    np.testing.assert_array_equal(
        two['detector'].reshape(2, -1)[0], one['detector'].reshape(1, -1)[0])


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

    # Total macroscopic cross section at the source energy
    library = openmc.data.DataLibrary.from_xml()
    h1 = openmc.data.IncidentNeutron.from_hdf5(
        library.get_by_material('H1')['path'])
    sigma_t = h1[1].xs[h1.temperatures[0]](2.0e6) * 0.05

    distance = np.linalg.norm(DETECTOR[0])
    expected = np.exp(-sigma_t * distance) / (4.0 * np.pi * distance**2)

    assert results['detector'].sum() == pytest.approx(expected, rel=0.02)


def test_point_detector_applies_particle_weight(run_in_tmpdir):
    """Survival biasing must not inflate the detector response.

    Under implicit capture the colliding particle carries a reduced weight;
    if the estimator scores at unit weight instead, the detector reads high by
    roughly the inverse of the surviving weight fraction. Analog and
    survival-biased runs are compared statistically, so the tolerance is loose
    -- the failure mode this guards against is order-unity, not marginal.
    """
    kwargs = dict(particles=4000, batches=20)
    analog, _ = _hydrogen_model([DETECTOR], survival=False, **kwargs)
    biased, _ = _hydrogen_model([DETECTOR], survival=True, **kwargs)

    analog.run(apply_tally_results=True)
    biased.run(apply_tally_results=True)

    a = analog.tallies[0]
    b = biased.tallies[0]
    sigma = np.hypot(a.std_dev.sum(), b.std_dev.sum())
    difference = abs(b.mean.sum() - a.mean.sum())

    assert difference < max(4.0 * sigma, 0.02 * a.mean.sum())
