# 14 MeV pencil-beam Li6 UQ scan

**Default updated:** the scan now uses the symmetric fractional ±10% profile.
See [li6_symmetric_uq.md](li6_symmetric_uq.md) for its exact shape-budget,
antithetic design, paired transport estimator and noise-aware validation.
The positive log-tilt configuration described below is retained via
`--profile log-tilt`; it is not the new default.

Run from the repository root:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan
```

The default scan uses 64 training points and 32 independent validation points,
800,000 histories per point (40,000 particles × 20 batches), a quadratic
21-term Legendre PCE, and 100,000 cheap surrogate draws for the output PDF/CDF.
It also runs an unperturbed baseline and a same-seed restored baseline.
One library session is used throughout, with exact native data restoration,
explicit RNG/tally reset and one transport thread. Runtime progress is printed
at every eight scan points. No nuclear-data HDF5 is rewritten or reprocessed.

The physical model is the existing pure Li6 sphere (20 cm radius, 0.46 g/cm3),
with a 14 MeV zero-width beam entering near z=-20 cm and travelling along +z.
The code reuses `make_model`, `TrialData`, and `run_trial` from
`fusion_li6_runtime_uq.py`, which apply the new native runtime APIs.

## What uncertainty is represented?

Default coefficients are independent uniform random variables:

| Parameter | Range | Effect |
|---|---|---|
| a | 0 to 0.5 | MT2 elastic XS multiplied by exp(a) |
| b1 | 0 to 4 | MT24 outgoing-energy linear shape tilt |
| b2 | 0 to 2 | MT24 outgoing-energy quadratic shape tilt |
| c1 | 0 to 2 | MT2 angular Legendre P1 tilt |
| c2 | 0 to 1 | MT2 angular Legendre P2 tilt |

The incident control-point window remains 10–14 MeV, with native interpolation
transitions at neighbouring unchanged points. MT2 outgoing energy is coupled
to angle through native elastic kinematics. The independently mutable outgoing
energy law is the discovered neutron ContinuousTabular law, MT24 in the
configured library. Product identity, yields, frames and native grids remain.

These ranges are a **declared sensitivity prior**, not an evaluated covariance
library. Change them with five-value `--lower-bounds` and `--upper-bounds`:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan \
  --lower-bounds 0 0 0 0 0 --upper-bounds 0.25 2 1 1 0.5 \
  --output sandbox/examples/fusion_li6_uq_scan_narrow
```

The scan draws reproducible Latin-hypercube training and validation designs
with independent seeds. Transport uses independent trial seeds. Arbitrary
correlated/evaluated priors need a corresponding basis/design and are not
silently represented by this independent-uniform PCE.

## Spatial outputs

- `beam_mean.png`: central x–z slices of baseline, prior-average production,
  and average minus baseline.
- `beam_variation.png`: estimated parameter standard deviation, coefficient
  of variation, and RMS transport standard error. Parameter variance is
  estimated as variance of sampled tally means minus average transport
  estimator variance, clipped at zero when the signal is unresolved.
- `beam_sample_extremes.png`: actual native mesh results at the sampled minimum,
  middle and maximum total production rates. These are observed sample extremes,
  not population confidence bounds or pointwise quantile fields.
- `beam_depth_variation.png`: slab-integrated production profile with estimated
  parameter standard deviation across the assumed prior.
- `scan_meshes.npz`: all sample meshes, native standard errors, parameter vectors,
  total rates/errors and baseline mesh.
- `spatial_statistics.npz`: full 3D mean, corrected variation, observed variance,
  transport-error diagnostic, and x/y/z edges.

Mesh values are `(n,t)` reactions per voxel per source neutron, not rates per
cm3. The plots use one central finite-thickness y=0 voxel slice; they do not sum
y-bin uncertainties as though those bins were independent. Zero-score bins
may have nonfinite native errors. Corrected statistics for those voxels remain
NaN; low-production cells are masked in relative-variation images. Spatial
maps are estimated from actual transport samples, not invented by a surrogate.

## Total-production distribution and PCE validation

PCE fits **log(total MT105 reaction rate per source neutron)** to preserve
positive predictions. A total-degree product Legendre basis is orthonormal
under the standardized independent uniform inputs. Training uses weighted
least squares with relative tally errors. Coefficient covariance estimates
account for training tally noise only, not all surrogate approximation error.

A separate 32-point transport design validates the model before a probability
plot is written. Acceptance requires held-out R² >= 0.90 and relative RMSE <=2%.
The metadata reports held-out error alongside RMS transport standard error.
If validation fails, the program saves diagnostics and mesh data, stops with
an error, and does not publish a newly validated density. Increase histories,
design size or polynomial order as indicated by the diagnostics.

`total_nt_distribution.png` shows the propagated-prior PCE PDF and empirical
CDF. The PDF is a KDE of 100,000 surrogate predictions; its histogram overlay
contains actual scan means, which still have transport noise. `metadata.json`
reports the predicted mean, standard deviation and 5th/50th/95th percentiles.
This is the distribution of the modeled expected production rate under the
assumed parameter prior, not a distribution of individual reaction counts.
Monte Carlo fitting error and PCE truncation error remain finite.

Other files:

- `scan.csv`: all input vectors, split, native total/scalar rates and errors,
  selected reactions/products, seeds, histories and status.
- `pce_model.npz`: polynomial indices, coefficients, coefficient covariance,
  physical input bounds and log-response flag.
- `pce_validation.csv`, `pce_validation.png`: independent prediction checks.
- `total_nt_distribution.csv`: PDF/CDF grid.
- `pce_distribution_samples.npz`: raw surrogate rate draws.
- `metadata.json`: prior, physics, native selections/interpolation windows,
  design/RNG settings, library index hash, OpenMC/Git identity and validation.

Evaluate the saved PCE without rerunning transport:

```python
import numpy as np
from sandbox.examples.fusion_li6_uq_scan import predict

model = np.load('sandbox/examples/fusion_li6_uq_scan_output/pce_model.npz')
x = np.array([[0.2, 2., 1., 1., 0.5]])  # a, b1, b2, c1, c2
u = (x - model['lower_bounds']) / (model['upper_bounds'] - model['lower_bounds'])
if np.any((u < 0) | (u > 1)):
    raise ValueError('Do not extrapolate outside the validated prior bounds')
rate = predict(u, model['multi_indices'], model['coefficients'])
print(rate)  # Expected (n,t) reactions per source neutron
```

For a larger cubic model (56 terms), use at least 112 training samples:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan \
  --degree 3 --train-samples 112 --validation-samples 40 \
  --output sandbox/examples/fusion_li6_uq_scan_cubic
```

Deterministic math tests:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg pytest -q tests/unit_tests/test_li6_uq_pce.py
```

They verify orthonormality by exact Gaussian quadrature, recovery of a known
log-polynomial response on unseen points, rank-deficient design rejection and
transport-noise subtraction. Full native scan validation and restoration
results are recorded in output metadata and run logs.

## More samples and cheaper PCE propagation

To double the native design and use one million surrogate draws:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan \
  --train-samples 128 --validation-samples 64 --surrogate-draws 1000000 \
  --output sandbox/examples/fusion_li6_uq_scan_refined_output
```

More native design points improve PCE fitting and spatial-prior estimates;
they do not reduce the statistical error of each transport point unless its
histories are also increased. More PCE draws reduce integration/sampling error
in the predicted mean, quantiles and density, conditioned on the fitted model.
They cannot fix model bias, polynomial truncation error or noisy training data.

After fitting/validation, reuse the saved PCE without any more transport:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.propagate_li6_pce \
  --model-directory sandbox/examples/fusion_li6_uq_scan_refined_output \
  --draws 1000000
```

This writes a separate `pce_propagation/` directory containing the PDF/CDF
figure, samples, CSV and `summary.json`. The summary distinguishes:

- Parameter-induced output standard deviation under the chosen prior.
- Monte Carlo standard error of the propagated mean, SD/sqrt(draw count).
- Approximate quantile sampling errors from the estimated density.
- Predicted-mean fit uncertainty from native training tally noise, propagated
  through the PCE coefficient covariance. This excludes truncation/model bias.

The utility refuses models that have not passed independent native validation.
It evaluates in bounded blocks to control memory use. Increasing the number of
cheap draws is useful until their sampling error becomes small relative to
model fitting/approximation uncertainty. Report the held-out validation metrics
alongside the PDF rather than interpreting a smooth density as proof of accuracy.
