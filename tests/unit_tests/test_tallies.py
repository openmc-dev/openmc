from math import sqrt
import numpy as np
import pytest
import openmc
import scipy.stats as sps


def test_xml_roundtrip(run_in_tmpdir):
    # Create a tally with all possible gizmos
    mesh = openmc.RegularMesh()
    mesh.lower_left = (-10., -10., -10.)
    mesh.upper_right = (10., 10., 10.,)
    mesh.dimension = (5, 5, 5)
    mesh_filter = openmc.MeshFilter(mesh)
    meshborn_filter = openmc.MeshBornFilter(mesh)
    tally = openmc.Tally()
    tally.filters = [mesh_filter, meshborn_filter]
    tally.nuclides = ['U235', 'I135', 'Li6']
    tally.scores = ['total', 'fission', 'heating']
    tally.derivative = openmc.TallyDerivative(
        variable='nuclide_density', material=1, nuclide='Li6'
    )
    tally.triggers = [openmc.Trigger('rel_err', 0.025)]
    tally.triggers[0].scores = ['total', 'fission']
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

    tally_a.name = 'new name'
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

    tally_a.scores = ['flux', 'absorption', 'fission', 'scatter']
    assert tally_a != tally_b
    tally_b.scores = ['flux', 'absorption', 'fission', 'scatter']
    assert tally_a == tally_b

    tally_a.nuclides = []
    tally_b.nuclides = []
    assert tally_a == tally_b

    tally_a.nuclides = ['total']
    assert tally_a == tally_b

    # a tally with an estimator set to None is equal to
    # a tally with an estimator specified
    tally_a.estimator = 'collision'
    assert tally_a == tally_b
    tally_b.estimator = 'collision'
    assert tally_a == tally_b

    tally_a.multiply_density = False
    assert tally_a != tally_b
    tally_b.multiply_density = False
    assert tally_a == tally_b

    trigger_a = openmc.Trigger('rel_err', 0.025)
    trigger_b = openmc.Trigger('rel_err', 0.025)

    tally_a.triggers = [trigger_a]
    assert tally_a != tally_b
    tally_b.triggers = [trigger_b]
    assert tally_a == tally_b


def test_figure_of_merit(sphere_model, run_in_tmpdir):
    # Run model with a few simple tally scores
    tally = openmc.Tally()
    tally.scores = ['total', 'absorption', 'scatter']
    sphere_model.tallies = [tally]
    sp_path = sphere_model.run(apply_tally_results=True)

    # Get execution time and relative error
    with openmc.StatePoint(sp_path) as sp:
        time = sp.runtime['simulation']
    rel_err = tally.std_dev / tally.mean

    # Check that figure of merit is calculated correctly
    assert tally.figure_of_merit == pytest.approx(1 / (rel_err**2 * time))


def test_tally_application(sphere_model, run_in_tmpdir):
    # Create a tally with most possible gizmos
    tally = openmc.Tally(name='test tally')
    ef = openmc.EnergyFilter([0.0, 0.1, 1.0, 10.0e6])
    mesh = openmc.RegularMesh.from_domain(sphere_model.geometry, (2, 2, 2))
    mf = openmc.MeshFilter(mesh)
    tally.filters = [ef, mf]
    tally.scores = ['flux', 'absorption', 'fission', 'scatter']
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
    assert tally._num_realizations == 0
    # the statepoint file property should be set, however
    assert tally._sp_filename == sp_file

    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert (sp_tally.mean == tally.mean).all()
    assert sp_tally.nuclides == tally.nuclides

    # SECOND RUN
    # change the number of particles and ensure that the results are different
    sphere_model.settings.particles += 1
    sp_file = sphere_model.run(apply_tally_results=True)

    assert (sp_tally.std_dev != tally.std_dev).any()
    assert (sp_tally.mean != tally.mean).any()

    # now re-read data from the new stateopint file and
    # ensure that the new results match those in
    # the latest statepoint
    with openmc.StatePoint(sp_file) as sp:
        assert tally in sp.tallies.values()
        sp_tally = sp.tallies[tally.id]

    # at this point the tally information regarding results should be the same
    assert (sp_tally.std_dev == tally.std_dev).all()
    assert (sp_tally.mean == tally.mean).all()
    assert sp_tally.nuclides == tally.nuclides

def _tally_from_data(x, *, vov_enabled=True, normality=True):
    t = openmc.Tally()
    t.scores = ["flux"]  # 1 score
    t.nuclides = [openmc.Nuclide("H1")]  # 1 nuclide
    t._sp_filename = "dummy.h5"  # mark "results available"
    t._results_read = True  # don't try to read from disk
    t._num_realizations = int(len(x))  # n
    t.vov_enabled = bool(vov_enabled)
    t.normality_tests = bool(normality)

    x = np.asarray(x, dtype=float)
    # (num_filter_bins=1, num_nuclides=1, num_scores=1) -> (1,1,1) arrays
    t._sum = np.array([[[np.sum(x)]]], dtype=float)
    t._sum_sq = np.array([[[np.sum(x**2)]]], dtype=float)
    if vov_enabled:
        t._sum_third = np.array([[[np.sum(x**3)]]], dtype=float)
        t._sum_fourth = np.array([[[np.sum(x**4)]]], dtype=float)
    return t

@pytest.mark.parametrize(
    "x, skew_true, kurt_true",
    [   # Rademacher distribution
        (np.array([1.0, -1.0] * 200), 0.0, 1.0),
        # Two-point {0,3} with p(0)=3/4, p(3)=1/4
        (np.concatenate([np.zeros(600), np.full(200, 3.0)]), 2.0 / sqrt(3.0), 7.0 / 3.0),
        # Bernoulli distribution
        (np.concatenate([np.ones(300), np.zeros(700)]), (1 - 2 * 0.3) / sqrt(0.3 * 0.7), (1 - 3 * 0.3 + 3 * 0.3**2) / (0.3 * 0.7)),
    ],
)
def test_b1_b2_analytical_against_tally(x, skew_true, kurt_true):
    t = _tally_from_data(x, vov_enabled=True, normality=False)

    g1 = t.skew(bias=True)[0, 0, 0]
    b2 = t.kurtosis(bias=True, fisher=False)[0, 0, 0]

    assert np.isclose(g1, skew_true, rtol=0, atol=1e-12)
    assert np.isclose(b2, kurt_true, rtol=0, atol=1e-12)

@pytest.mark.parametrize(
    "draw, skew_true, kurt_true",
    [(lambda rng, n: rng.normal(0, 1, n), 0.0, 3.0),  # Normal
     (lambda rng, n: rng.random(n), 0.0, 1.8),  # Uniform(0,1)
     (lambda rng, n: rng.exponential(1.0, n), 2.0, 9.0),  # Exp(1)
     (lambda rng, n: (rng.random(n) < 0.3).astype(float),
            (1 - 2 * 0.3) / sqrt(0.3 * 0.7),
            (1 - 3 * 0.3 + 3 * 0.3**2) / (0.3 * 0.7),),],)

def test_b1_b2_scipy_and_theory(draw, skew_true, kurt_true):
    rng = np.random.default_rng(12345)
    N = 200_000
    x = draw(rng, N)

    # Tally outputs
    t = _tally_from_data(x, vov_enabled=True, normality=False)
    g1_t = t.skew(bias=True)[0, 0, 0]
    b2_t = t.kurtosis(bias=True, fisher=False)[0, 0, 0]

    # SciPy (population, bias=True to match population-moment style)
    skew_sp = sps.skew(x, bias=True)
    kurt_sp = sps.kurtosis(x, fisher=False, bias=True)

    # Compare to SciPy numerically
    assert np.isclose(g1_t, skew_sp, rtol=0, atol=5e-3)
    assert np.isclose(b2_t, kurt_sp, rtol=0, atol=5e-3)

    # Compare to analytical targets with size-dependent tolerances
    tol_skew = 0.02 if abs(skew_true) < 0.5 else 0.05
    tol_kurt = 0.03 if kurt_true < 4 else 0.1
    assert abs(g1_t - skew_true) < tol_skew
    assert abs(b2_t - kurt_true) < tol_kurt


def test_kurtosis_bias_fisher_combinations():
    """Test that all combinations of bias and fisher match scipy.stats.kurtosis"""
    rng = np.random.default_rng(42)
    x = rng.normal(0, 1, 10000)

    t = _tally_from_data(x, vov_enabled=True, normality=False)

    # Test all four combinations
    # 1. bias=True, fisher=False (Pearson's kurtosis, b2)
    b2_tally = t.kurtosis(bias=True, fisher=False)[0, 0, 0]
    b2_scipy = sps.kurtosis(x, fisher=False, bias=True)
    assert np.isclose(b2_tally, b2_scipy, rtol=0, atol=1e-10)
    assert np.isclose(b2_tally, 3.0, rtol=0.05, atol=0.1)  # Should be ~3 for normal

    # 2. bias=True, fisher=True (excess kurtosis, g2)
    g2_tally = t.kurtosis(bias=True, fisher=True)[0, 0, 0]
    g2_scipy = sps.kurtosis(x, fisher=True, bias=True)
    assert np.isclose(g2_tally, g2_scipy, rtol=0, atol=1e-10)
    assert np.isclose(g2_tally, 0.0, rtol=0, atol=0.1)  # Should be ~0 for normal
    assert np.isclose(g2_tally, b2_tally - 3.0, rtol=0, atol=1e-10)  # g2 = b2 - 3

    # 3. bias=False, fisher=True (adjusted excess kurtosis, G2)
    G2_tally = t.kurtosis(bias=False, fisher=True)[0, 0, 0]
    G2_tally_default = t.kurtosis()[0, 0, 0]  # Should be same as default
    G2_scipy = sps.kurtosis(x, fisher=True, bias=False)
    assert np.isclose(G2_tally, G2_tally_default, rtol=0, atol=1e-10)
    assert np.isclose(G2_tally, G2_scipy, rtol=0, atol=1e-10)
    assert np.isclose(G2_tally, 0.0, rtol=0, atol=0.1)  # Should be ~0 for normal

    # 4. bias=False, fisher=False (adjusted Pearson's kurtosis)
    adj_b2_tally = t.kurtosis(bias=False, fisher=False)[0, 0, 0]
    adj_b2_scipy = sps.kurtosis(x, fisher=False, bias=False)
    assert np.isclose(adj_b2_tally, adj_b2_scipy, rtol=0, atol=1e-10)
    assert np.isclose(adj_b2_tally, 3.0, rtol=0.05, atol=0.1)  # Should be ~3 for normal
    assert np.isclose(adj_b2_tally, G2_tally + 3.0, rtol=0, atol=1e-10)  # adj_b2 = G2 + 3


def test_ztests_scipy_comparison():
    rng = np.random.default_rng(987)
    x_norm = rng.normal(size=50_000)
    x_exp = rng.exponential(size=50_000)

    # -------- Normal dataset (should not reject) --------
    t0 = _tally_from_data(x_norm, vov_enabled=True, normality=True)
    stats0 = t0.normality_test(alternative="two-sided")

    Zb1_0 = stats0["Zb1"][0, 0, 0]
    p_skew_0 = stats0["p_skew"][0, 0, 0]
    Zb2_0 = stats0["Zb2"][0, 0, 0]
    p_kurt_0 = stats0["p_kurt"][0, 0, 0]
    K2_0 = stats0["K2"][0, 0, 0]
    p_omni_0 = stats0["p_K2"][0, 0, 0]

    z_skew_sp0, p_skew_sp0 = sps.skewtest(x_norm)
    z_kurt_sp0, p_kurt_sp0 = sps.kurtosistest(x_norm)
    k2_sp0, p_omni_sp0 = sps.normaltest(x_norm)

    assert np.isclose(Zb1_0, z_skew_sp0, atol=0.15)
    assert np.isclose(Zb2_0, z_kurt_sp0, atol=0.15)
    assert np.isclose(K2_0, k2_sp0, atol=0.30)
    assert np.isclose(p_skew_0, p_skew_sp0, atol=5e-3)
    assert np.isclose(p_kurt_0, p_kurt_sp0, atol=5e-3)
    assert np.isclose(p_omni_0, p_omni_sp0, atol=5e-3)

    # -------- Exponential dataset (should strongly reject) --------
    t1 = _tally_from_data(x_exp, vov_enabled=True, normality=True)
    stats1 = t1.normality_test(alternative="two-sided")

    Zb1_1 = stats1["Zb1"][0, 0, 0]
    p_skew_1 = stats1["p_skew"][0, 0, 0]
    Zb2_1 = stats1["Zb2"][0, 0, 0]
    p_kurt_1 = stats1["p_kurt"][0, 0, 0]
    K2_1 = stats1["K2"][0, 0, 0]
    p_omni_1 = stats1["p_K2"][0, 0, 0]

    z_skew_sp1, p_skew_sp1 = sps.skewtest(x_exp)
    z_kurt_sp1, p_kurt_sp1 = sps.kurtosistest(x_exp)
    k2_sp1, p_omni_sp1 = sps.normaltest(x_exp)

    # Both pipelines should reject very strongly
    assert p_skew_1 < 1e-6 and p_skew_sp1 < 1e-6
    assert p_kurt_1 < 1e-6 and p_kurt_sp1 < 1e-6
    assert p_omni_1 < 1e-6 and p_omni_sp1 < 1e-6

    # Right-skewed and heavy-tailed → large positive Z-statistics
    assert Zb1_1 > 30 and z_skew_sp1 > 30
    assert Zb2_1 > 30 and z_kurt_sp1 > 30
    assert K2_1 > 2000 and k2_sp1 > 2000

def test_vov_stochastic(sphere_model, run_in_tmpdir):
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

    num = (sum_fourth - (4.0*sum_third*sum_)/n + (6.0*sum_sq*sum_**2)/(n**2)
           - (3.0*sum_**4)/(n**3))
    den = (sum_sq - (1.0/n)*sum_**2)**2

    expected_vov[nonzero] = num[nonzero]/den[nonzero] - 1.0/n

    assert np.allclose(expected_vov, sp_tally.vov, rtol=1e-7, atol=0.0)
