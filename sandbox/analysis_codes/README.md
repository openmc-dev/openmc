# Shared OpenMC nuclear-data analysis engine

For the detailed suite description, configuration reference and operating guide, see [description.md](description.md).

This directory contains the analysis implementation and its tests. It imports
OpenMC and standard scientific dependencies, but **does not import `Li_uq`,
`FLiBe_uq`, `FLiBe_analysis` or `nuclear_data_sensitivity`**. Those older suites
remain available for historical studies. The OpenMC runtime mutation APIs are
still supplied by this modified OpenMC installation.

A study is a directory with `config.py` and `model.py`. All results, checkpoints,
plots and nuclear-data provenance are written inside that study's configured
results directory. No geometry is defined in this engine.

From `/sandbox/openmc`, using the configured OpenMC Python environment:

```bash
export OPENMC_CROSS_SECTIONS=/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/opt/openmc-venv/bin/python -m sandbox.analysis_codes FLiBe_analysis --dry-run
/opt/openmc-venv/bin/python -m sandbox.analysis_codes FLiBe_analysis --run
/opt/openmc-venv/bin/python -m sandbox.analysis_codes FLiBe_analysis --resume
/opt/openmc-venv/bin/python -m sandbox.analysis_codes FLiBe_analysis --reanalyse
```

From any working directory, use the direct-script entry point:

```bash
/opt/openmc-venv/bin/python /sandbox/openmc/sandbox/analysis_codes/run.py /path/to/study --run
```

`--run` and `--resume` both reuse successful deterministic native checkpoints;
`--resume` is an explicit spelling of the same operation. `--dry-run` performs
native discovery and estimates histories, but does not run transport.
`--reanalyse` uses saved observations and performs no transport. Changing
transport inputs requires a new study folder or `--output /path/to/new/results`.
PCE settings may change for saved-data reanalysis. The old completed FLiBe
campaign is adopted only after matching its configuration and model builder;
its original transport provenance/checkpoints remain untouched. Partially
completed legacy campaigns cannot be resumed as new-engine checkpoints.

## Create another study

```bash
mkdir /sandbox/openmc/sandbox/studies/Li_analysis
cp /sandbox/openmc/sandbox/studies/FLiBe_analysis/model.py /sandbox/openmc/sandbox/studies/FLiBe_analysis/config.py /sandbox/openmc/sandbox/studies/Li_analysis/
# Edit these two files, then:
python -m sandbox.analysis_codes Li_analysis --dry-run
python -m sandbox.analysis_codes Li_analysis --run
```

Do not copy a completed study's results when changing its model. `model.py`
exports `build_model(particles, batches) -> openmc.Model`. It may import its
local settings with `from . import config as c`; directories with spaces or
outside this repository are supported. Define any geometry/material/source in
this function. Return a fresh model on each call. Avoid transport or library
initialisation at module import time. Use `openmc.reset_auto_ids()` before
manual repeated model construction.

`config.py` controls low-fidelity histories as `SCREEN_PARTICLES * SCREEN_BATCHES`
and high-fidelity histories as `HIGH_PARTICLES * HIGH_BATCHES`. `WORKERS` controls
Morris processes; `TRANSPORT_THREADS` controls threads in each native instance.
`HIGH_WORKERS` controls high-fidelity processes (default 1). Parallel workers use spawned, isolated native libraries and a shared atomic history budget; each worker resets and restores its native state for every evaluation.
`TOP_PER_FAMILY` selects up to n parameters independently in XS, energy and
angle; explicit `RESERVED_PARAMETERS` can add physics checks. Every retained
parameter is included in subsequent PCE/Sobol output.

### PCE backend

Set in the study's `config.py`:

```python
PCE_BACKEND = 'numpy'       # Existing weighted ridge/sparse implementation
# PCE_BACKEND = 'openturns' # Optional OpenTURNS weighted LARS/CLOO selection
```

Install OpenTURNS only if requested: `python -m pip install openturns`.
Both backends use orthonormal Legendre polynomials, declared hyperbolic and
interaction truncation, training-only CV and independent OpenMC validation.
OpenTURNS selects support with LARS/CLOO on an explicitly prewhitened design;
optional retained main effects and ridge are refitted on that support. This
is an intentionally specified workflow, not a claim that both regressions are
identical. Its predictions and coefficient Sobol indices are checked against
OpenTURNS' native chaos representation. Covariance estimates are conditional
on the selected support and omit model-selection/truncation uncertainty.

Results in `results/high/` reflect the last backend analysed; a separate copy
of each backend's fitted model, validation, Sobol results and distribution is
kept in `results/pce_backends/numpy/` or `openturns/`. Transport datasets are
shared. Switching backends does **not** require new OpenMC calculations.

### Tritium scoring and other models

The built-in response is the verified ENDF/B-VIII.0 partial-reaction sum:
Li6 MT105, Li7 MT52–82, Be9 MT700/701 and F19 MT105, each yielding one triton.
It scores only constituents present in the model. Cached MT205 and redundant
Be9 MT105 are not primary responses. Contributions are summed within each
batch before estimating standard errors. Data hashes, source ENDF metadata,
Q values and native full-grid baseline reconstruction are checked. Reference
ENDF files are cached in the **study results**, not another study's directory.
A fresh study can fetch the bounded NNDC archive prefix with pinned SHA checks.

For a Li-only model, leave the verified mapping intact or restrict
`PRODUCTION_CHANNELS` to its Li entries. New nuclides added to PRODUCTION_CHANNELS require established yield
mappings; other material nuclides need no production mapping. For a different or unsupported production response, define `responses(model) -> dict[str, Response]` in
its `model.py`, using `analysis_codes.Response`, and
`verify_response(output_directory)` to perform controlled response verification.
That callback must write `verification.json` with `status='passed'` and the
current `analysis_codes.fingerprint.implementation_sha()` as
`implementation_sha256`. Always include a scalar response named `total` first.
Model-specific response code belongs with the model; no engine change is needed.

## Workflow and modules

1. `context.py`, `config.py`: load/register arbitrary study models and responses.
2. `production.py`: verified tritium scoring and native controlled tests.
3. `discovery.py`, `thresholds.py`, `catalogue.py`: native inventory, remove
   thresholds **strictly above** source energy, plot all retained representations.
4. `registry.py`, `cross_sections.py`, `perturbations.py`: stable parameter IDs
   and representation-aware candidates; no writable native pointers.
5. `runner.py`, `budget.py`, `storage.py`: native execution, exact rollback,
   per-batch response reduction, atomic checkpoints and hard history budgets.
6. `morris.py`, `analysis.py`, `settings.py`: global trajectories, paired-batch
   errors, bootstrap intervals, convergence diagnostics and low-fidelity run.
7. `family_workflow.py`: per-family selection, high-fidelity designs and saved
   data analysis; `pce.py` / `openturns_pce.py`, `sobol.py`, `validation.py`:
   surrogate regression, coefficient sensitivity decomposition and held-out checks.
8. `native_plots.py`, `two_pass.py`: native perturbation galleries and complete
   Morris/Sobol figures. Unsupported correlated PDFs are documented, not invented.
9. `workflow.py`, `__main__.py`, `run.py`: orchestration and entry points;
   `fingerprint.py`: implementation provenance.

Independent coefficients in `[-PERTURBATION_LIMIT, +PERTURBATION_LIMIT]` are hypothetical priors, not covariance data.
PDF coefficients share the existing two-mode budget: a single mode has at most
half the configured local PDF bound, whereas both together stay within the full configured bound. XS coefficients
scale the native cross section directly. These are different physical operations.
Set any finite `0 < PERTURBATION_LIMIT < 1` in the study configuration; for example, `0.10` means ±10% XS coefficients and `0.30` means ±30%. Discovery, Morris, high-fidelity transport, native verification, plots and report labels all use this setting. A changed limit requires a fresh output directory; saved-data reanalysis rejects a range that differs from the saved parameter registry. Native
interpolation, thresholds, support and frames are preserved. Unsupported
mutations are listed. No conclusion of convergence is made solely from successful
execution; unresolved Morris/PCE/Sobol components remain labelled.

## Tests

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest analysis_codes/tests -q
```

## Additional configuration controls

Both study configs now expose the following settings. Geometry and source values
are consumed by that study's model builder; statistical settings are passed by
the shared workflow. Older study configs may omit new statistical keys and retain
the original defaults.

| Settings | Effect |
|---|---|
| `DENSITY_G_CM3`, `TEMPERATURE_K`, `LI6_ATOM_FRACTION` | Actual material density, temperature and Li isotopic composition; `None` uses natural Li. |
| `RADIUS_CM`, `WORLD_MARGIN_CM` | Target radius and vacuum-world margin. |
| `SOURCE_POSITION_CM` | Explicit source point in cm; `None` uses the automatic external -z position. |
| `SOURCE_STANDOFF_CM`, `SOURCE_DIRECTION`, `SOURCE_ENERGY_EV` | Automatic-position standoff, beam direction and monoenergetic source energy. |
| `MORRIS_MIN_TRAJECTORIES`, `MORRIS_RANK_CORRELATION`, `MORRIS_TOP_OVERLAP`, `MORRIS_MAX_RELATIVE_CI_WIDTH` | Minimum trajectory count and ranking/precision thresholds for the declared convergence diagnostic. |
| `MORRIS_BOOTSTRAPS`, `MORRIS_PROGRESS_BOOTSTRAPS`, `MORRIS_CONVERGENCE_BOOTSTRAPS` | Resamples for final results, progress results and convergence diagnostics. |
| `MORRIS_CONFIDENCE_LEVEL`, `MORRIS_NOISE_MULTIPLIER` | Bootstrap interval coverage and noise-resolution multiplier. |
| `MORRIS_BOOTSTRAP_SEED` | Bootstrap seed; `None` derives it from `SEED`. |
| `PCE_MIN_R2`, `PCE_MAX_NOISE_MSE`, `PCE_MIN_SIGNAL_TO_NOISE` | Validation requires either the R² criterion or noise-consistency criterion, and always requires resolved signal. Constant/unavailable responses remain rejected. |
| `NATIVE_PLOT_SEED` | Native diagnostic-sampler seed; `None` uses `SEED`. |
| `WORKERS`, `HIGH_WORKERS`, `TRANSPORT_THREADS` | Morris processes, high-fidelity processes and OpenMC threads per process. |

`MORRIS_TRAJECTORIES` still sets the planned trajectory count; changing convergence
criteria does not silently authorize extra transport histories. Increasing worker
counts does not increase histories per sample. Allow sufficient CPUs for workers
multiplied by threads.

**Li model correction:** its actual density and temperature now come from config.
The old builder replaced the configured material with a natural-Li material at
0.534 g/cm³ and did not assign its temperature. The current Li config's 1.94 g/cm³
and 900 K will therefore be used on a new run unless you edit them. Historical
results remain unchanged and must not be described using the corrected model.
Set the density you intend before launching new transport. New physical settings
require a new output directory; saved observations cannot be relabelled as a new
material/source calculation.


## Production QoI versus transport nuclides

`PRODUCTION_CHANNELS` defines which verified reaction rates are summed into the
tritium-production QoI. For example, Li6 MT105 plus Li7 MT52–82 measures lithium
production even when lead is present. Only the configured scored channels require
production-yield verification. No exclusion list is needed.

`TARGET_NUCLIDES` independently controls perturbation discovery. With `None`,
all material nuclides, including lead, remain eligible for safely supported
perturbations. Their transport changes can affect the configured lithium QoI.
The configured sum is not a claim about production from unlisted nuclides.
