import numpy as np
import pytest
import openmc
import openmc.tally_stats as ts
import os, h5py, pprint
import shutil
from fractions import Fraction


def test_xml_roundtrip(run_in_tmpdir):
    # Create a tally with all possible gizmos
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10.0, -10.0, -10.0)
    mesh.upper_right = (
        10.0,
        10.0,
        10.0,
    )
    mesh.dimension = (5, 5, 5)
    mesh_filter = openmc.MeshFilter(mesh)
    meshborn_filter = openmc.MeshBornFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter, meshborn_filter]
    tally.nuclides = ["U235", "I135", "Li6"]
    tally.scores = ["total", "fission", "heating"]
    tally.derivative = openmc.TallyDerivative(
        variable="nuclide_density", material=1, nuclide="Li6"
    )
    tally.triggers = [openmc.Trigger("rel_err", 0.025)]
    tally.triggers[0].scores = ["total", "fission"]
    tallies = openmc.Tallies([tally])

    # Roundtrip through XML and make sure we get what we started with
    tallies.export_to_xml()
    new_tallies = openmc.Tallies.from_xml()
    assert len(new_tallies) == 1
    new_tally = new_tallies[0]
    assert new_tally.id == tally.id
    assert len(new_tally.filters) == 2
    assert isinstance(new_tally.filters[0], openmc.MeshFilter)
    assert np.allclose(new_tally.filters[0].mesh.lower_left, mesh.lower_left)
    assert isinstance(new_tally.filters[1], openmc.MeshBornFilter)
    assert np.allclose(new_tally.filters[1].mesh.lower_left, mesh.lower_left)
    assert new_tally.nuclides == tally.nuclides
    assert new_tally.scores == tally.scores
    assert new_tally.derivative.variable == tally.derivative.variable
    assert new_tally.derivative.material == tally.derivative.material
    assert new_tally.derivative.nuclide == tally.derivative.nuclide
    assert len(new_tally.triggers) == 1
    assert new_tally.triggers[0].trigger_type == tally.triggers[0].trigger_type
    assert new_tally.triggers[0].threshold == tally.triggers[0].threshold
    assert new_tally.triggers[0].scores == tally.triggers[0].scores
    assert new_tally.multiply_density == tally.multiply_density


def test_tally_equivalence():
    tally_a = openmc.Tally()
    tally_b = openmc.Tally(tally_id=tally_a.id)

    tally_a.name = "new name"
    assert tally_a != tally_b
    tally_b.name = tally_a.name
    assert tally_a == tally_b

    assert tally_a == tally_b
    ef_a = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    ef_b = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])

    tally_a.filters = [ef_a]
    assert tally_a != tally_b
    tally_b.filters = [ef_b]
    assert tally_a == tally_b

    tally_a.scores = ["flux", "absorption", "fission", "scatter"]
    assert tally_a != tally_b
    tally_b.scores = ["flux", "absorption", "fission", "scatter"]
    assert tally_a == tally_b

    tally_a.nuclides = []
    tally_b.nuclides = []
    assert tally_a == tally_b

    tally_a.nuclides = ["total"]
    assert tally_a == tally_b

    # a tally with an estimator set to None is equal to
    # a tally with an estimator specified
    tally_a.estimator = "collision"
    assert tally_a == tally_b
    tally_b.estimator = "collision"
    assert tally_a == tally_b

    tally_a.multiply_density = False
    assert tally_a != tally_b
    tally_b.multiply_density = False
    assert tally_a == tally_b

    trigger_a = openmc.Trigger("rel_err", 0.025)
    trigger_b = openmc.Trigger("rel_err", 0.025)

    tally_a.triggers = [trigger_a]
    assert tally_a != tally_b
    tally_b.triggers = [trigger_b]
    assert tally_a == tally_b


def test_figure_of_merit(sphere_model, run_in_tmpdir):
    # Run model with a few simple tally scores
    tally = openmc.Tally()
    tally.scores = ["total", "absorption", "scatter"]
    sphere_model.tallies = [tally]
    sp_path = sphere_model.run(apply_tally_results=True)

    # Get execution time and relative error
    with openmc.StatePoint(sp_path) as sp:
        time = sp.runtime["simulation"]
    rel_err = tally.std_dev / tally.mean

    # Check that figure of merit is calculated correctly
    assert tally.figure_of_merit == pytest.approx(1 / (rel_err**2 * time))


def test_tally_normality_functions():
    values = np.arange(1, 15, dtype=float)
    n = values.size
    mean = np.array([values.mean()])
    sum_sq = np.array([np.sum(values**2)])
    sum_third = np.array([np.sum(values**3)])
    sum_fourth = np.array([np.sum(values**4)])

    sqrt_b1, b2 = ts._calc_b1_b2(n, mean, sum_sq, sum_third, sum_fourth)
    assert sqrt_b1.shape == mean.shape
    assert b2.shape == mean.shape

    Zb1, p_skew, _ = ts.skewness_test(n, mean, sum_sq, sum_third, sum_fourth)
    assert Zb1.shape == mean.shape
    assert p_skew.shape == mean.shape
    assert np.all((0.0 <= p_skew) & (p_skew <= 1.0))

    Zb2, p_kurt, _ = ts.kurtosis_test(n, mean, sum_sq, sum_third, sum_fourth)
    assert Zb2.shape == mean.shape
    assert p_kurt.shape == mean.shape
    assert np.all((0.0 <= p_kurt) & (p_kurt <= 1.0))

    K2, p_omni = ts.k2_test(Zb1, Zb2)
    assert K2.shape == mean.shape
    assert p_omni.shape == mean.shape
    assert np.all((0.0 <= p_omni) & (p_omni <= 1.0))


@pytest.mark.parametrize(
    "x, expected",
    [
        ([-1, 1], Fraction(0, 1)),
        ([-1, 0, 1], Fraction(1, 6)),
        ([0, 1, 2], Fraction(1, 6)),
        ([0, 1, 2, 3], Fraction(4, 25)),
        ([0, 1, 2, 3, 4], Fraction(7, 50)),
        ([0, 0, 1, 1], Fraction(0, 1)),
    ],
)
def test_vov_deterministic(x, expected, run_in_tmpdir):
    x = np.asarray(x, dtype=float)
    N = x.size
    S1 = x.sum()
    S2 = (x**2).sum()
    S3 = (x**3).sum()
    S4 = (x**4).sum()

    t = openmc.Tally()
    t._sp_filename = run_in_tmpdir / "dummy.sp"  # allow property access
    t.num_realizations = N
    t.sum = S1
    t.sum_sq = S2
    t.sum_third = S3
    t.sum_fourth = S4

    vov = t.vov[0]
    assert np.isclose(vov, float(expected), rtol=0, atol=1e-15)


def test_tally_normality_stats(sphere_model, run_in_tmpdir):

    tally = openmc.Tally()
    tally.scores = ["flux"]
    tally.vov_enabled = True
    tally.normality_tests = True
    sphere_model.tallies = [tally]

    sphere_model.settings.particles = 500
    sphere_model.settings.batches = 30
    sphere_model.settings.run_mode = "fixed source"

    sp_file = sphere_model.run(apply_tally_results=True)

    n = tally.num_realizations
    mu = tally.mean
    s2 = tally.sum_sq
    s3 = tally.sum_third
    s4 = tally.sum_fourth

    assert n >= 20
    assert mu is not None

    Zg1, p_skew, _ = ts.skewness_test(n, mu, s2, s3, s4)
    Zg2, p_kurt, _ = ts.kurtosis_test(n, mu, s2, s3, s4)
    K2, p_omni = ts.k2_test(Zg1, Zg2)

    for arr in (Zg1, Zg2, K2, p_skew, p_kurt, p_omni):
        assert arr.shape == mu.shape

    for p in (p_skew, p_kurt, p_omni):
        assert ((0.0 <= p) & (p <= 1.0)).all()


def test_vov_stochastic(sphere_model, run_in_tmpdir):
    # Create a tally with all the gizmos
    tally = openmc.Tally(name="test tally")
    ef = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    mesh = openmc.RegularMesh.from_domain(sphere_model.geometry, (2, 2, 2))
    mf = openmc.MeshFilter(mesh)
    tally.filters = [ef, mf]
    tally.scores = ["flux", "absorption", "fission", "scatter"]
    tally.vov_enabled = True
    sphere_model.tallies = [tally]

    sp_file = sphere_model.run(apply_tally_results=True)

    assert tally._mean is None
    assert tally._std_dev is None
    assert tally._sum is None
    assert tally._sum_sq is None
    assert tally._sum_third is None
    assert tally._sum_fourth is None
    assert tally._num_realizations == 0
    assert tally._sp_filename == sp_file

    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    assert np.all(sp_tally.std_dev == tally.std_dev)
    assert np.all(sp_tally.mean == tally.mean)
    assert np.all(sp_tally.vov == tally.vov)
    assert sp_tally.nuclides == tally.nuclides

    n = sp_tally.num_realizations
    mean = sp_tally.mean
    sum_ = sp_tally._sum
    sum_sq = sp_tally._sum_sq
    sum_third = sp_tally._sum_third
    sum_fourth = sp_tally._sum_fourth

    expected_vov = np.zeros_like(mean)
    nonzero = np.abs(mean) > 0

    num = (
        sum_fourth
        - (4.0 * sum_third * sum_) / n
        + (6.0 * sum_sq * sum_**2) / (n**2)
        - (3.0 * sum_**4) / (n**3)
    )
    den = (sum_sq - (1.0 / n) * sum_**2) ** 2

    expected_vov[nonzero] = num[nonzero] / den[nonzero] - 1.0 / n

    assert np.allclose(expected_vov, sp_tally.vov, rtol=1e-7, atol=0.0)


def test_tally_application(sphere_model, run_in_tmpdir):
    # Create a tally with most possible gizmos
    tally = openmc.Tally(name="test tally")
    ef = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    mesh = openmc.RegularMesh.from_domain(sphere_model.geometry, (2, 2, 2))
    mf = openmc.MeshFilter(mesh)
    tally.filters = [ef, mf]
    tally.scores = ["flux", "absorption", "fission", "scatter"]
    tally.vov_enabled = True
    sphere_model.tallies = [tally]

    # FIRST RUN
    # run the simulation and apply results
    sp_file = sphere_model.run(apply_tally_results=True)
    # before calling for any property requiring results (including the equivalence check below),
    # the following internal attributes of the original should be unset
    assert tally._mean is None
    assert tally._std_dev is None
    assert tally._sum is None
    assert tally._sum_sq is None
    assert tally._sum_third is None
    assert tally._sum_fourth is None
    assert tally._num_realizations == 0
    # the statepoint file property should be set, however
    assert tally._sp_filename == sp_file

    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert (sp_tally.mean == tally.mean).all()
    assert (sp_tally.vov == tally.vov).all()
    assert sp_tally.nuclides == tally.nuclides

    # SECOND RUN
    # change the number of particles and ensure that the results are different
    sphere_model.settings.particles += 1
    sp_file = sphere_model.run(apply_tally_results=True)

    assert (sp_tally.std_dev != tally.std_dev).any()
    assert (sp_tally.mean != tally.mean).any()
    assert (sp_tally.vov != tally.vov).any()

    # now re-read data from the new stateopint file and
    # ensure that the new results match those in
    # the latest statepoint
    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert (sp_tally.mean == tally.mean).all()
    assert (sp_tally.vov == tally.vov).all()
    assert sp_tally.nuclides == tally.nuclides
