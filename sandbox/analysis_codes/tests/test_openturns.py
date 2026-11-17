"""Analytic checks independent of OpenMC and saved transport observations."""
import builtins
import numpy as np
import pytest
from sandbox.analysis_codes.pce import PCE, matrix, indices
from sandbox.analysis_codes.sobol import analyse
from sandbox.analysis_codes import openturns_pce as backend


def test_missing_dependency_is_actionable(monkeypatch):
    original = builtins.__import__

    def missing(name, *args, **kwargs):
        if name == 'openturns':
            raise ImportError('simulated optional package absence')
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', missing)
    with pytest.raises(ImportError, match='pip install openturns'):
        backend.require_openturns()


def test_known_coefficients_predictions_and_native_sobol(tmp_path):
    pytest.importorskip('openturns')
    rng = np.random.default_rng(19)
    x = rng.random((180, 2))
    powers = indices(2, 2)
    c = np.array([1., .2, .4, .1, .15, .05])
    y = matrix(x, powers) @ c
    fit = backend.fit(x, y, np.full(len(x), .001), order=2)
    np.testing.assert_allclose(fit.coefficients[:, 0], c, atol=1e-10)
    np.testing.assert_allclose(fit.predict(x)[:, 0], y, atol=1e-10)
    stats = analyse(fit)
    native = fit.diagnostics['native_sobol'][0]
    np.testing.assert_allclose(stats['first'][:, 0], native['first'], atol=1e-12)
    np.testing.assert_allclose(stats['total'][:, 0], native['total'], atol=1e-12)
    fit.save(tmp_path/'model.npz')
    restored = PCE.load(tmp_path/'model.npz')
    np.testing.assert_array_equal(restored.predict(x), fit.predict(x))
    # Covariance must propagate supplied errors through the selected WLS fit.
    a = matrix(x, powers)
    expected = np.linalg.pinv(a) @ np.linalg.pinv(a).T * 1e-6
    np.testing.assert_allclose(fit.covariance[0], expected, atol=1e-14)


def test_heteroscedastic_errors_and_invalid_input():
    pytest.importorskip('openturns')
    x = np.linspace(.01, .99, 100)[:, None]
    y = 1 + .2*x[:, 0]
    y[-1] += 10
    errors = np.full(100, .001); errors[-1] = 100
    fit = backend.fit(x, y, errors, order=1)
    np.testing.assert_allclose(fit.predict(x)[:-1, 0], 1+.2*x[:-1, 0], atol=1e-5)
    for bad in (-np.ones(100), np.full(100, np.nan)):
        with pytest.raises(ValueError):
            backend.fit(x, y, bad)
    with pytest.raises(ValueError):
        backend.fit(x+2, y, errors)
    with pytest.raises(ValueError):
        backend.fit(x, y*np.inf, errors)
    unavailable = backend.fit(x, y, np.full(100, np.inf))
    assert unavailable.diagnostics['unavailable_bins'] == [0]


def test_training_only_order_selection():
    pytest.importorskip('openturns')
    x = np.random.default_rng(29).random((100, 1))
    y = 1 + .2*x[:, 0]**2
    order, scores = backend.cross_validate(x, y, np.full(100, .001), orders=(1, 2))
    assert order == 2 and scores[2] < scores[1]


def test_truncated_basis_and_retained_main_effects():
    pytest.importorskip('openturns')
    x=np.random.default_rng(23).random((150,3))
    powers=indices(3,3,q=.75,interaction=2)
    c=np.zeros(len(powers))
    for p,v in [((0,0,0),2),((1,0,0),.3),((0,1,0),1e-4),((1,1,0),.01)]:
        c[np.flatnonzero((powers==p).all(axis=1))[0]]=v
    model=backend.fit(x,matrix(x,powers)@c,np.full(150,1e-6),order=3,
        q=.75,interaction=2,keep_main_effects=True,ridge=1e-8)
    np.testing.assert_allclose(model.predict(x)[:,0],matrix(x,powers)@c,atol=1e-8)
    result=analyse(model)
    np.testing.assert_allclose(result['first'][:,0],model.diagnostics['native_sobol'][0]['first'],atol=1e-12)
