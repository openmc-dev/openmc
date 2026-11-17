# Symmetric ±10% Li6 scan

This is the default profile of `fusion_li6_uq_scan.py`:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan \
  --profile symmetric10 --surrogate-draws 1000000
```

Outputs go to `sandbox/examples/fusion_li6_uq_symmetric10_output/`, preserving earlier
one-sided scans. Use `--profile log-tilt` to explicitly request the old shape
parameterization. The new default prior is five independent uniform fractional
parameters in [-0.1,+0.1]. Training and validation designs are separate
antithetic Latin-hypercube designs: every input vector has its exact negative.
No inputs are selected or rejected based on the direction of the tally response.
The intended prior is independent; the finite experimental design is deliberately
paired to balance its sampled first moments exactly.

## Meaning of ±10%

- Li6 MT2 XS: `sigma_new = sigma_original*(1+delta_xs)` at native incident grid
  points from 10 through 14 MeV. The multiplier is uniform from 0.9 to 1.1 and
  has mean one. This avoids the upward arithmetic-mean shift of exp(a) with a
  symmetric prior for a.
- Native MT24 neutron outgoing-energy PDF and native MT2 angular PDF: bounded
  shape perturbations in the same incident-table window. For each table,
  native-interpolation-centered modes g1 and g2 have zero PDF-weighted mean
  and max absolute value at most one. The new PDF is
  `p_new=p_baseline*[1+(delta1*g1+delta2*g2)/2]`.
- Each delta is uniform in [-0.1,+0.1]. The two modes share a common 10% budget,
  so their *combined* pointwise change is bounded by ±10%, rather than allowing
  two 10% changes to add to 20%. Every positive coefficient vector has an
  equal negative counterpart. PDF expectation is the baseline up to native
  evaluated-table normalization roundoff. Zero-density points remain zero;
  relative changes are undefined at those points.
- Energy modes are centered/scaled versions of z and z², where
  z=E_out/max(table E_out). Angular modes are centered/scaled P1 and P2.
- Continuous mass is preserved with native histogram/linear-linear integration;
  discrete outgoing-energy masses are unchanged. Native loaded PDFs/CDFs may
  differ at roundoff precision, so in-window PDFs are normalized to their
  correct continuous mass. Native readback tests permit only 1e-6 rounding
  tolerance on the 0.9/1.1 ratio limits.

An outgoing energy value itself is not independently multiplied by a random
number; that could violate kinematics. MT2 has no separate energy PDF: its
outgoing energy changes through native angle–energy kinematics. Independent
energy-PDF shape changes remain on the discovered MT24 law. Native incident
interpolation transitions at adjacent out-of-window points remain documented.

This is a baseline-centered sensitivity prior, not evaluated nuclear-data
covariance. Symmetric inputs do not force symmetric responses or a mean
response exactly equal to baseline: nonlinear physical effects can cause a
real second-order shift.

## Transport noise and common-random-number control variates

±10% density-shape changes have much smaller effects than the previous large
one-sided tilts. Each sample therefore uses two native runs with identical
source histories/RNG seed: an unperturbed control and its perturbed counterpart.
The batch-by-batch difference supplies the SE of the difference, including the
covariance of the paired scores. Its estimator is:

`rate_estimate = high_statistics_baseline + (perturbed - matched_control)`.

The common baseline uses 400,000 particles × 20 batches = 8 million histories;
each member of a matched pair uses 40,000 × 20 = 800,000 histories. With 96
parameter points, this is 192 native scan runs plus baseline/restored-baseline
runs. The baseline common error shifts all response estimates together and is
reported separately; it is not added independently to every point or mistaken
for parameter-induced output variance. Its absolute statistical error must be
included when quoting the overall mean rate. `propagate_li6_pce.py` reports
mean-fit uncertainty with and without this common baseline error.

`scan.csv` retains raw perturbed/control rates, their difference/paired SE,
the common baseline estimate/SE, and the response used for PCE. Raw native
transport has not been replaced by a particle sampler. The control-variate
estimator can be negative in very weak-production mesh voxels because of
statistical noise; those values are preserved in arrays, not silently clipped.

The default PCE is linear (six terms) for these small fractional changes; the
quadratic model's extra terms were not resolved above training noise. The
initial rejected model/diagnostics are retained. Cached results can be refitted
without transport:

```bash
python -m sandbox.examples.fusion_li6_uq_scan \
  --refit-only sandbox/examples/fusion_li6_uq_symmetric10_output --degree 1 \
  --surrogate-draws 1000000
```

PCE validation uses independent parameter points and independent matched
history pairs. Because small-response spread can be comparable with transport
SE, raw R² alone is not an adequate test. The symmetric profile reports raw
R² but also requires mean squared standardized residual <=2, response spread
more than three times the fitted-response noise diagnostic, and relative RMSE
<=0.5%. Standardization includes held-out paired SE and training-fit prediction
SE. This is a noise-aware adequacy check, not proof that all polynomial
truncation bias is absent. The density is withheld if those checks fail.

The PDF describes parameter-induced expected `(n,t)` production *conditional on
the common baseline estimate*. It does not include counting noise or silently
fold the common baseline's Monte Carlo error into nuclear-data uncertainty.

## Checks and files

Tests cover antithetic design symmetry, exact centering in synthetic histogram
and linear-linear tables, discrete mass preservation, native ±10% readbacks,
unchanged out-of-window tables, unmodified MT105 XS, exact restoration, and
batch-estimate agreement with native tally means.

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg pytest -q \
  tests/unit_tests/test_li6_symmetric_uq.py tests/unit_tests/test_li6_uq_pce.py
```

The output files follow the existing scan: native 3D mesh data and errors,
mean/variation/sample-extreme beam images, scan CSV, independent PCE validation,
saved PCE, and total-production probability PDF/CDF. Earlier one-sided results
remain in their original directories for comparison.

## Larger ±20% scan

The `symmetric20` profile uses the same centered, linear construction with
coefficients uniform in [-0.2,+0.2]. XS multipliers range from 0.8 to 1.2;
the two energy modes share a 20% PDF budget, as do the two angular modes.
Normalization, discrete mass and native interpolation rules are unchanged.
Run the larger design separately to preserve the ±10% results:

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.fusion_li6_uq_scan \
  --profile symmetric20 --train-samples 128 --validation-samples 64 \
  --surrogate-draws 1000000 \
  --output sandbox/examples/fusion_li6_uq_symmetric20_output
```

This uses 192 parameter points (twice the earlier design), each with a matched
unperturbed transport control. Each transport member has 800,000 histories;
including the baseline and restoration runs gives 323.2 million histories.

## Separate effects on (n,t)

```bash
OMP_NUM_THREADS=2 MPLBACKEND=Agg python -m sandbox.examples.li6_separate_effects
```

This runs six native one-at-a-time cases and their matched controls: XS only,
energy shape only and angle shape only, each at −20% and +20%. Other classes
remain at baseline. The energy and angular comparisons move their respective
two shape coefficients together, preserving the shared 20% budget. They are
specific directions in parameter space, not independent prior averages.

`separate_effects/*_spatial.png` shows the baseline central voxel slice and
signed differences; `*_response.png` shows depth profiles and total-rate changes.
Error bars are one paired standard error and exclude the separately saved
common-baseline error. Finite-statistics spatial maps are not significance maps.
The associated `results.json` and `native_responses.npz` preserve raw totals,
paired differences, errors, parameters, seeds and spatial arrays.
