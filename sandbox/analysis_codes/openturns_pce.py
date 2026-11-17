"""Optional OpenTURNS weighted sparse PCE, adapted to existing diagnostics.

No transport is performed here. Coefficients use the same orthonormal Legendre
basis as the existing PCE/Sobol machinery. Noise covariance is conditional on
the selected basis; it excludes uncertainty in LARS selection and truncation.
"""
import numpy as np
from sandbox.analysis_codes.pce import PCE, indices, matrix


def require_openturns():
    """Import the optional dependency, with an actionable error message."""
    try:
        import openturns as ot
    except ImportError as error:
        raise ImportError(
            'OpenTURNS was requested. Install it with '
            '`python -m pip install openturns`, or use --use_openturns off.'
        ) from error
    return ot


def fit(unit, response, errors=None, *, order=2, q=1., interaction=None, ridge=0., sparse=True, keep_main_effects=False):
    """Fit each response with inverse-variance weighted OpenTURNS LARS/CLOO.

    Parameters
    ----------
    unit : numpy.ndarray
        Independent uniform input coordinates in [0, 1], shape (n, d).
    response : numpy.ndarray
        Training observations, shape (n,) or (n, responses).
    errors : numpy.ndarray, optional
        Transport standard errors. Infinite errors exclude observations;
        zero errors receive the existing adapter's positive weight floor.
    order : int
        Maximum total polynomial degree.

    Returns
    -------
    analysis_codes.pce.PCE
        Portable coefficients, conditional noise covariance and OT diagnostics.
    """
    ot = require_openturns()
    x = np.asarray(unit, float)
    if x.ndim != 2 or x.shape[1] < 1 or len(x) < 2:
        raise ValueError('Need a two-dimensional input array and >=2 samples')
    powers = indices(x.shape[1], order, q=q, interaction=interaction)
    a = matrix(x, powers)  # Also validates finite inputs and the unit domain.
    y = np.asarray(response, float)
    if y.ndim == 1:
        y = y[:, None]
    if y.ndim != 2 or len(y) != len(x) or not np.all(np.isfinite(y)):
        raise ValueError('Training responses must be finite and match inputs')
    s = np.ones_like(y) if errors is None else np.asarray(errors, float)
    if s.size != y.size:
        raise ValueError('Standard errors must match responses')
    s = s.reshape(y.shape)
    if np.any(np.isnan(s)) or np.any(s < 0):
        raise ValueError('Standard errors must be nonnegative and not NaN')
    distribution = ot.JointDistribution([ot.Uniform(0., 1.)] * x.shape[1])
    enumeration = ot.LinearEnumerateFunction(x.shape[1])
    basis = ot.OrthogonalProductPolynomialFactory(
        [ot.LegendreFactory()] * x.shape[1], enumeration)
    native_basis_ids = [int(enumeration.inverse(p.tolist())) for p in powers]
    size = len(powers)
    lookup = {tuple(p): i for i, p in enumerate(powers)}
    coef = np.zeros((len(powers), y.shape[1]))
    covariance = np.zeros((y.shape[1], len(powers), len(powers)))
    selected_counts, conditions, unavailable, direct_sobol = [], [], [], []
    prediction_discrepancies = []
    for column in range(y.shape[1]):
        finite = np.isfinite(s[:, column])
        if finite.sum() < 2:
            unavailable.append(column)
            selected_counts.append(0)
            conditions.append(float('inf'))
            direct_sobol.append(None)
            prediction_discrepancies.append(None)
            continue
        sigma = s[finite, column]
        positive = sigma[sigma > 0]
        floor = max(.1 * np.median(positive) if len(positive) else 1., 1e-300)
        # Rescale first to avoid overflow; common scaling leaves WLS unchanged.
        root_weights = np.minimum(np.maximum(sigma, floor).min()
                                  / np.maximum(sigma, floor), 1.)
        # Prewhiten BOTH design and observations. Explicitly exposing the
        # whitened design avoids the weighted sparse-selection degeneracy
        # observed with the installed OpenTURNS 1.27.post1 high-level path.
        # The extra coordinate carries known row weights, not an uncertain
        # nuclear-data parameter. It is absent from the final PCE and Sobol basis.
        root_weights /= np.sqrt(np.mean(root_weights**2))
        symbols = [f'z{i}' for i in range(x.shape[1])] + ['wgt']
        coordinate_projection = ot.SymbolicFunction(symbols, symbols[:-1])
        weight_function = ot.SymbolicFunction(symbols, ['wgt'])
        product = ot.SymbolicFunction(['p', 'w'], ['p*w'])
        weighted_basis = [ot.ComposedFunction(product, ot.AggregatedFunction([
            ot.ComposedFunction(basis.build(native_basis_ids[i]), coordinate_projection),
            weight_function])) for i in range(size)]
        augmented = np.column_stack([2*x[finite]-1, root_weights])
        algorithm = ot.LeastSquaresMetaModelSelectionFactory(
            ot.LARS(), ot.CorrectedLeaveOneOut()).build(
                ot.Sample(augmented),
                ot.Sample(y[finite, column:column+1]*root_weights[:, None]),
                ot.Point(int(finite.sum()), 1.), weighted_basis, list(range(size)))
        algorithm.run()
        native_coefficients = np.asarray(algorithm.getCoefficients())
        native_ids = np.flatnonzero(native_coefficients).tolist()
        if 0 not in native_ids:
            native_ids.insert(0, 0)
        selected = native_ids if sparse else list(range(len(powers)))
        if keep_main_effects:
            selected = sorted(set(selected) | set(np.flatnonzero((powers>0).sum(axis=1)<=1)))
        design = a[finite][:,selected] * root_weights[:,None]
        target = y[finite,column] * root_weights
        if ridge:
            penalty=np.sqrt(ridge)*np.eye(len(selected));penalty[0,0]=0
            design_aug=np.vstack([design,penalty]);target_aug=np.r_[target,np.zeros(len(selected))]
        else:design_aug,target_aug=design,target
        estimator_aug=np.linalg.pinv(design_aug)[:,:len(target)]
        fitted=estimator_aug @ target
        coef[selected,column]=fitted
        native_ids=[native_basis_ids[i] for i in selected]
        native_coefficients=fitted
        # Reconstruct the ordinary orthonormal PCE, excluding the auxiliary
        # regression-weight coordinate. OT independently evaluates it and its
        # Sobol decomposition below, checking the portable coefficient mapping.
        result = ot.FunctionalChaosResult(
            ot.Sample(x[finite]), ot.Sample(y[finite, column:column+1]), distribution,
            ot.DistributionTransformation(distribution, basis.getMeasure()),
            ot.DistributionTransformation(basis.getMeasure(), distribution),
            basis, native_ids, ot.Sample(native_coefficients[:, None]),
            [basis.build(i) for i in native_ids])
        result.setIsLeastSquares(True)
        result.setInvolvesModelSelection(True)
        # Propagate supplied independent transport noise through the selected
        # weighted LS estimator. Do not replace the OpenTURNS coefficients.
        design = a[finite][:, selected]
        weighted = design * root_weights[:, None]
        estimator = estimator_aug * root_weights[None, :]
        influence = estimator * sigma[None, :]
        covariance[column][np.ix_(selected, selected)] = influence @ influence.T
        selected_counts.append(len(selected))
        conditions.append(float(np.linalg.cond(weighted)))
        native_prediction = np.asarray(result.getMetaModel()(ot.Sample(x)))[:, 0]
        discrepancy = float(np.max(abs(a @ coef[:, column] - native_prediction)))
        if not np.allclose(a @ coef[:, column], native_prediction, rtol=1e-10, atol=1e-12):
            raise RuntimeError('OpenTURNS coefficient/basis mapping mismatch')
        prediction_discrepancies.append(discrepancy)
        if np.sum(coef[1:, column]**2) > 0:
            native_indices = ot.FunctionalChaosSobolIndices(result)
            direct_sobol.append(dict(
                first=[native_indices.getSobolIndex(i) for i in range(x.shape[1])],
                total=[native_indices.getSobolTotalIndex(i) for i in range(x.shape[1])]))
        else:
            direct_sobol.append(None)
    return PCE(powers, coef, covariance, dict(
        backend='openturns', version=ot.__version__, terms=len(powers),
        samples=len(x), selected_terms=selected_counts, condition_numbers=conditions,
        unavailable_bins=unavailable, sparse=sparse,
        method='OpenTURNS LARS + corrected leave-one-out on explicitly prewhitened design',
        weighting='Inverse transport-error variance; auxiliary row-weight coordinate excluded from PCE/Sobol',
        covariance='Transport noise conditional on selected support; excludes model selection/truncation',
        native_sobol=direct_sobol,
        maximum_native_prediction_discrepancy=prediction_discrepancies))


def cross_validate(unit, response, errors, *, orders=(2, 3), folds=4, seed=20261020, **options):
    """Select degree using training-only folds; never consume validation data."""
    x, y, s = map(np.asarray, (unit, response, errors))
    if y.ndim == 1:
        y=y[:,None]
    s=s.reshape(y.shape)
    if not 2 <= folds <= len(x):
        raise ValueError('Need 2 <= folds <= training samples')
    splits = np.array_split(np.random.default_rng(seed).permutation(len(x)), folds)
    scores = {}
    for order in orders:
        residuals = []
        for held in splits:
            train = np.setdiff1d(np.arange(len(x)), held)
            model = fit(x[train], y[train], s[train], order=order, **options)
            if model.diagnostics['unavailable_bins']:
                raise ValueError('Insufficient finite training errors in a response')
            residuals.append(np.mean((model.predict(x[held]) - y[held])**2))
        scores[order] = float(np.mean(residuals))
    return min(scores, key=scores.get), scores
