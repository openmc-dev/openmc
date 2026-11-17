# Description of the OpenMC nuclear-data sensitivity code suite

This document describes the reusable analysis code in `analysis_codes/` and its
use with study directories such as `Li_analysis/` and `FLiBe_analysis/`.
It describes the implemented workflow; a feature mentioned as a limitation is
not implied to be available. Settings and results belong to a specific study.

## 1. Purpose and organisation

The suite measures how declared nuclear-data perturbations affect neutron-induced
tritium production. It uses native OpenMC transport for every physical evaluation,
Morris screening to rank parameters, a polynomial chaos expansion (PCE) to
approximate selected responses, and coefficient-based Sobol indices to attribute
response variance. Its input ranges are hypothetical sensitivity priors, not
an evaluated nuclear-data covariance model.

```text
/sandbox/openmc/
├── analysis_codes/             Shared implementation, documentation and tests
├── Li_analysis/
│   ├── model.py                Lithium model construction
│   ├── config.py               Physical and numerical settings
│   ├── README.md               Study notes, if supplied
│   └── results/                Data, checkpoints, plots and reports
└── FLiBe_analysis/
    ├── model.py                FLiBe model construction
    ├── config.py               Physical and numerical settings
    ├── README.md
    └── results/
```

Change the shared implementation in `analysis_codes/` when improving the engine.
Change a study's `model.py` and `config.py` when changing the problem being
analysed. Do not put study geometry into the engine. Historical suites such as
`Li_uq/` and `FLiBe_uq/` remain available but are not analysis-code dependencies
of the new engine. The engine still requires this modified OpenMC installation
and its native runtime mutation interfaces.

## 2. End-to-end workflow

1. **Load the study.** Import its local config and model builder, select its
   output directory and check requested settings against existing output metadata.
2. **Verify production scoring.** Check the supported evaluation-specific triton
   reaction mapping, nuclear-data hashes, ENDF metadata and baseline reconstruction.
3. **Discover native data.** Query loaded nuclides, reactions, products and native
   distributions. Distinguish inspectable data from safely mutable data.
4. **Filter by source energy.** Remove reactions whose native threshold is strictly
   greater than the maximum supported discrete source energy, including their
   associated energy and angle parameters. Write the exclusion list.
5. **Plot the retained native data.** Catalogue plots include supported/read-only
   representations with explicit limitations. They precede sensitivity selection.
6. **Check the budget.** Calculate planned Morris, high-fidelity and verification
   histories before starting the campaign. A dry run stops before transport.
7. **Run controlled native verification.** Perturb known production/scattering
   cross sections, check production sums and confirm exact baseline restoration.
8. **Run low-fidelity Morris screening.** Evaluate complete trajectories using
   the configured lower history count. Calculate elementary effects, intervals,
   noise diagnostics and ranking-convergence measures.
9. **Select inputs per family.** Keep up to `TOP_PER_FAMILY` inputs separately in
   cross section, energy and angle; add available configured reserved parameters.
10. **Generate new high-fidelity observations.** Use separate Latin-hypercube
    training and held-out validation designs for the selected parameter space.
11. **Fit and validate PCE.** Select polynomial order using training-only CV,
    then assess predictions against independent OpenMC validation observations.
12. **Calculate conditional Sobol indices.** Use orthonormal PCE coefficients,
    mark unresolved components and reject responses that fail validation.
13. **Generate response distributions.** Cheaply sample an accepted PCE and plot
    its response histogram, validation observations, baseline, mean and percentiles.
14. **Write reports and plots.** Preserve datasets, provenance and backend-specific
    analysis outputs inside the study's results directory.

The Morris observations are used for screening and selection. They are **not
appended to the current high-fidelity PCE training design**. PCE training uses
new simulations with the high-fidelity history settings.

## 3. Physical response and units

For the verified ENDF/B-VIII.0 evaluations, the built-in response sums:

| Constituent | Native partial reaction scores | Verified triton yield |
|---|---|---|
| Li6 | MT105 | One per reaction |
| Li7 | MT52 through MT82 | One per reaction |
| Be9 | MT700 and MT701 | One per reaction |
| F19 | MT105 | One per reaction |

Only constituents present in the model are scored. The engine verifies and
locks these mappings to the recorded evaluation/library identity. It does not
infer a general triton yield solely from an MT label, and it refuses to invent
mappings for new constituents explicitly requested in PRODUCTION_CHANNELS.

The primary response is **tritons per source neutron**. OpenMC tracklength
reaction-rate tallies include the appropriate target atom density; the suite
does not multiply an already macroscopic reaction tally by density again.
Constituent contributions are summed within each batch before calculating a
standard error, retaining their statistical covariance. The suite also reports
constituent contributions separately.

Cached MT205 is not the perturbed production response: it can remain unchanged
when partial cross sections are mutated. Redundant Be9 MT105 is also excluded
from the primary reaction sum in favour of its partial channels. Baseline MT205
may be used only as an independent diagnostic of the verified mapping.

These outputs are not tritons per second. For a steady source strength of
`S` neutrons/second, multiply a rate per source neutron by `S`. No radioactive
decay production or tritium transport is included.

A closed reaction channel can legitimately give exactly zero production. For
example, the loaded Li7 tritium channels begin above a 2 MeV source energy.
Zero-variance responses have undefined R²/noise diagnostics and are not meaningful
Sobol responses; this does not indicate failed neutron transport. Li7 elastic
scattering can still affect Li6 production indirectly through the neutron field.

## 4. Meaning of the perturbations

Let `L = PERTURBATION_LIMIT`. Independent coefficients have a uniform prior
between `-L` and `+L`, with `0 < L < 1` required.

### Cross sections

The standard whole-reaction operation is

```text
sigma_new(E) = sigma_baseline(E) * (1 + delta),  delta in [-L, L].
```

Thus `L=0.10` gives ±10% cross-section coefficients and `L=0.30` gives ±30%.
Native incident grids are retained, zero cross sections stay zero, and the
threshold is preserved. Native derived-data rebuilding and eligibility checks
are used; unsupported/dependent targets are not silently changed.

The underlying adapter also implements exponential scaling and disjoint
incident-energy masks through `Parameter`. The default study workflow does not
expose these as additional analysis methods or automatically create energy groups.

### Emission-energy distributions

Supported independent continuous-tabular neutron laws use bounded outgoing-energy
shape modes: a linear mode and a quadratic mode. The incident-energy grid and
outgoing-energy support remain native. The operation changes the PDF shape,
not the reaction cross section, multiplicity or a prescribed percentage of the
mean outgoing energy.

### Angular distributions

Supported independent tabular angular laws use `mu=cos(theta)` and P1/P2 modes:

```text
P1(mu) = mu
P2(mu) = (3*mu*mu - 1)/2
```

The existing strict shape operation centres each basis under the native PDF,
bounds it and shares the budget between the two declared modes. Schematically,

```text
p_new(x) = p_baseline(x) * [1 + (delta_1*f_1(x) + delta_2*f_2(x))/2]
```

with native normalization and CDF rebuilding. One mode at its endpoint changes
the local continuous PDF by at most `L/2`; both modes together remain within
`L`. The division by two remains even if screening excludes one mode, so a
parameter does not change meaning between Morris and PCE. Native normalization
roundoff is checked, and an inadmissible candidate is rejected.

Discrete probability masses, native interpolation, support and reference frames
are preserved. Histogram areas use histogram integration, not a generic
trapezoidal integral. The primary angular plots show density in mu; they are
not probability densities per degree.

Isotropic/analytic/correlated laws are not silently converted into independent
tabular PDFs. An inspectable law is not necessarily mutable. Unsupported native
representations and reasons are included in the catalogue.

## 5. Morris screening and selection

For `d` eligible parameters and `r` trajectories, the main Morris design needs
`r*(d+1)` OpenMC evaluations. The current low stage also has one baseline and
three timing-pilot observations. Each evaluation uses

```text
SCREEN_PARTICLES * SCREEN_BATCHES histories.
```

An elementary effect is a response difference divided by the signed input step
on the unit input cube. The implementation records signed mean `mu`, mean
absolute effect `mu_star`, variability `sigma`, bootstrap intervals and noise
resolution diagnostics. Raw magnitudes depend on the response units and prior
range; they are not Sobol variance fractions.

With `PAIRED_MORRIS=True`, vertices within a trajectory share a transport seed;
trajectories have different seeds. Effect standard errors use observed paired
batch differences, including covariance. Common random numbers are not assumed
to improve every comparison. With pairing disabled, independent-error handling
is used. Bootstrap resampling treats a whole trajectory as its resampling unit.

Convergence diagnostics compare ranking halves, top-ten overlap and interval
width, and require resolved leading effects and a minimum trajectory count.
`MORRIS_TRAJECTORIES` sets the planned count in this workflow. Changing the
convergence thresholds does not silently authorize extra trajectories. If the
criterion fails, the report must retain that limitation.

Selection is by `mu_star` within each family. A noisy large absolute effect can
enter the selected list: the code retains its unresolved flag and does not claim
that selection proves importance. Reserved physics checks can increase the
selected dimension beyond three times the per-family limit. If a family has
fewer than n eligible modes, all available modes are kept.

## 6. High-fidelity PCE, validation and Sobol

High-fidelity training and validation use new, reproducible Latin-hypercube
designs. Each evaluation uses `HIGH_PARTICLES * HIGH_BATCHES` histories. There
is also a high-fidelity baseline observation. Transport seeds are deterministic
functions of `SEED`, stage and sample index; validation is separate from training.

PCE uses orthonormal Legendre polynomials on independent uniform coordinates.
Hyperbolic truncation, maximum interaction order and a basis-size cap control
the candidate basis. `PCE_ORDERS` provides candidate degrees; training-only
cross-validation chooses among allowed candidates. The basis cap is not a
substitute for an adequate training design. Underdetermined unregularized dense
fits are rejected by the regression implementation.

Two selectable backends are available:

- `numpy`: the existing weighted ridge/sparse regression implementation.
- `openturns`: optional OpenTURNS LARS/CLOO support selection on a prewhitened
  design, followed by the configured support/main-effect and ridge handling.
  Native OpenTURNS predictions and Sobol values check coefficient mapping.

The regressions need not produce identical models. Transport-error weighting
helps avoid fitting Monte Carlo noise, but fitted coefficient covariance is
conditional on the selected support and does not capture all truncation or
model-selection uncertainty.

Independent validation reports R², RMSE, normalized RMSE, maximum error and
residuals relative to estimated transport/prediction uncertainty. Acceptance
requires either the configured R² criterion or noise-consistent residuals,
and always requires resolved parameter-induced signal. Constant or unavailable
responses remain rejected. A successful script exit alone does not validate a PCE.

For a validated orthonormal expansion, squared nonconstant coefficients give
variance contributions. First-order indices retain terms involving only one
input; total-order indices retain all terms involving that input. Family/group
indices are computed using polynomial supports, rather than by summing
potentially overlapping individual total-order indices. Weak components remain
marked unresolved. An additive fitted basis can give equal first/total indices;
that does not prove physical interactions are absent.

Sobol results are conditional on the selected parameter set, with excluded
inputs fixed at baseline. The current workflow does not establish that excluded
inputs have zero influence or that reduced-space indices describe full evaluated
nuclear-data uncertainty.

Response distributions are generated only for an accepted total-response PCE.
A large cheap surrogate sample produces the mean, standard deviation and 5th/95th
percentiles. Independent OpenMC validation points are displayed as observations,
not as a precise direct-Monte-Carlo probability density.

## 7. Configuration reference

Edit the study's `config.py`. The shared engine's `settings.py` supplies internal
backward-compatible defaults; it is not the normal place to configure a study.
The following settings are supported by the supplied study models and workflow.

| Setting | Meaning |
|---|---|
| `RADIUS_CM` | Target sphere radius in cm. |
| `DENSITY_G_CM3` | Actual target material density in g/cm³; used by the model builder. |
| `TEMPERATURE_K` | Actual material temperature in kelvin. |
| `SOURCE_ENERGY_EV` | Monoenergetic neutron source energy in eV; the source definition in model.py must use it. |
| `LI6_ATOM_FRACTION` | Li6 fraction within lithium. None means natural isotopic lithium; otherwise use a fraction between 0 and 1. |
| `TOP_PER_FAMILY` | Maximum ranked parameters selected independently in XS, energy and angle before adding reserved checks. |
| `MORRIS_TRAJECTORIES` | Planned number of complete global Morris trajectories. |
| `MORRIS_LEVELS` | Morris grid levels; the implemented design requires an even value of at least four. |
| `SCREEN_PARTICLES` | Source histories per batch during low-fidelity screening. |
| `SCREEN_BATCHES` | Number of low-fidelity batches; at least two are required. |
| `HIGH_PARTICLES` | Source histories per batch during high-fidelity training, validation and baseline runs. |
| `HIGH_BATCHES` | Number of high-fidelity batches; at least two are required. |
| `TRAINING_SAMPLES` | Number of new high-fidelity PCE training observations. |
| `VALIDATION_SAMPLES` | Number of independent held-out high-fidelity validation observations. |
| `SEED` | Base for reproducible transport and experimental-design seeds. Stages and sample indices add deterministic offsets. |
| `HISTORY_CAP` | Explicit total study history cap. The stage caps must also fit inside it. |
| `LOW_HISTORY_CAP` | Hard reservation budget for low-fidelity observations, including baseline and timing pilot. |
| `HIGH_HISTORY_CAP` | Hard reservation budget for high-fidelity baseline, training and validation. |
| `VERIFICATION_HISTORY_CAP` | Hard reservation budget for controlled native verification. |
| `WORKERS` | Number of spawned Morris transport processes. |
| `ANALYSIS_SEED` | Seed for training-only PCE cross-validation folds and cheap surrogate propagation. |
| `STUDY_NAME` | Study/model identifier used in metadata; use a single path component. |
| `TARGET_NUCLIDES` | Nuclides whose nuclear data are screened. None selects all material nuclides; this does not itself define which nuclides contribute to the response. |
| `TRANSPORT_THREADS` | OpenMC threads per native worker/session. |
| `RESULTS_DIRECTORY` | Results path relative to the study folder, unless absolute. CLI --output overrides it. |
| `PERTURBATION_LIMIT` | Fractional prior bound L, finite and strictly between 0 and 1. Changing it changes parameter identities and requires new transport data. |
| `PAIRED_MORRIS` | Use a common transport seed within each trajectory and paired batch errors. |
| `PCE_BACKEND` | numpy or openturns. OpenTURNS must be installed when requested. |
| `PCE_ORDERS` | Candidate polynomial degrees; selection uses training-only cross-validation. |
| `PCE_CV_FOLDS` | Number of training-only cross-validation folds. |
| `PCE_HYPERBOLIC_Q` | Hyperbolic truncation exponent, 0 < q <= 1. Smaller q discourages high-order interactions. |
| `PCE_MAX_INTERACTION` | Maximum number of distinct inputs involved in a candidate polynomial term. |
| `PCE_RIDGE` | Nonnegative regression regularization strength. |
| `PCE_SPARSE` | Enable support reduction/sparse fitting according to the selected backend. |
| `PCE_KEEP_MAIN_EFFECTS` | Retain univariate polynomial terms during sparse fitting so weak families are not removed solely relative to a dominant input. |
| `PCE_BASIS_LIMIT` | Upper bound on candidate polynomial terms; degrees exceeding the cap are excluded. |
| `PCE_DISTRIBUTION_SAMPLES` | Number of cheap surrogate evaluations for the response distribution; these are not OpenMC histories. |
| `PRODUCTION_CHANNELS` | Verified constituent-to-partial-MT mapping. The built-in adapter permits the established evaluation-specific mappings, not arbitrary unverified yield assumptions. |
| `RESPONSE_LABEL` | Text label for the plotted/reported response. It does not change tally physics or units. |
| `VERIFICATION_PARTICLES` | Histories per batch in native controlled verification. |
| `VERIFICATION_BATCHES` | Batches per native controlled verification case. |
| `RESERVED_PARAMETERS` | Additional available parameters identified by (nuclide, MT, family, mode), retained after per-family ranking. |
| `WORLD_MARGIN_CM` | Vacuum-world radius minus target radius, in cm. |
| `SOURCE_POSITION_CM` | Explicit (x,y,z) source point in cm. None selects the automatic point on the negative z side. |
| `SOURCE_STANDOFF_CM` | Distance outside the negative-z target surface for automatic source positioning; keep the point inside the simulated world. |
| `SOURCE_DIRECTION` | Beam direction vector used by OpenMC Monodirectional. Changing direction alone does not automatically reposition the source. |
| `MORRIS_MIN_TRAJECTORIES` | Minimum trajectory count for declaring the configured convergence criterion met. |
| `MORRIS_RANK_CORRELATION` | Minimum split-half ranking Spearman correlation for convergence. |
| `MORRIS_TOP_OVERLAP` | Minimum top-ten ranking overlap fraction for convergence. |
| `MORRIS_MAX_RELATIVE_CI_WIDTH` | Maximum median relative bootstrap interval width among leading effects for convergence. |
| `MORRIS_BOOTSTRAPS` | Number of trajectory bootstrap resamples for final Morris estimates. |
| `MORRIS_PROGRESS_BOOTSTRAPS` | Bootstrap resamples for intermediate Morris estimates. |
| `MORRIS_CONVERGENCE_BOOTSTRAPS` | Bootstrap resamples used within trajectory-convergence diagnostics. |
| `MORRIS_CONFIDENCE_LEVEL` | Bootstrap confidence coverage, strictly between 0 and 1. |
| `MORRIS_NOISE_MULTIPLIER` | Noise-resolution multiplier for Morris effect flags. |
| `MORRIS_BOOTSTRAP_SEED` | Bootstrap RNG seed. None derives it from SEED. |
| `PCE_MIN_R2` | Minimum held-out R² for the strong-fit acceptance branch. |
| `PCE_MAX_NOISE_MSE` | Maximum mean squared standardized residual for the noise-consistency acceptance branch. |
| `PCE_MIN_SIGNAL_TO_NOISE` | Minimum parameter-spread/coefficient-noise ratio required by both acceptance branches. |
| `NATIVE_PLOT_SEED` | Seed for native diagnostic-spectrum sampling in catalogue plots. None uses SEED. |
| `HIGH_WORKERS` | Number of spawned high-fidelity transport processes; 1 uses serial repeated runs. |

The supplied sphere models consume the geometry/source settings above. A custom
blanket builder must explicitly read the settings it needs; the engine cannot
apply a sphere radius or change arbitrary geometry by interpreting a config name.
Avoid retaining unused model-specific settings in a new study.

**Important historical Li model correction:** the previous lithium builder
created a configured material and then replaced it with a natural-Li material
at 0.534 g/cm³, without assigning the configured temperature. The current builder
uses `DENSITY_G_CM3`, `TEMPERATURE_K` and `LI6_ATOM_FRACTION`. Check the intended
values before launching a new study. Old observations retain their original
physical meaning and cannot be relabelled as runs of the corrected material.

## 8. What each module does

| File in analysis_codes | Responsibility |
|---|---|
| `__init__.py` | Exposes ModelSpec, Response and Parameter to user/model code. |
| `__main__.py` | Parses study path, mode and optional output override for python -m sandbox.analysis_codes. |
| `run.py` | Direct-script entry point, including use outside the repository working directory. |
| `context.py` | Loads a study's local model/config under an isolated module name, delegates construction and response hooks, and registers target nuclides/threads. |
| `workflow.py` | Coordinates stages, checks request provenance/budgets, performs per-family selection, invokes fitting and writes study reports/backend archives. |
| `config.py` | Defines ModelSpec and Response, including within-batch weighted tally reduction and optional volume scaling. This is not the editable study config.py. |
| `settings.py` | Internal dataclass carrying Morris execution/statistical settings between functions. |
| `production.py` | Verifies the supported triton mappings, caches pinned reference ENDF evidence, builds constituent/total tracklength tallies, and runs controlled native checks. |
| `discovery.py` | Enumerates loaded reactions and secondary neutron laws; tests native backup eligibility and reports unsupported representations. |
| `registry.py` | Parameter metadata, family/mode/interval semantics, coefficient bounds and stable hashed identities. |
| `thresholds.py` | Determines the supported discrete source-energy maximum, compares native thresholds and writes kept/removed lists. |
| `catalogue.py` | Saves complete native reaction/distribution inventories and baseline galleries for source-retained data. |
| `cross_sections.py` | Constructs nonnegative, threshold-preserving XS candidates on native grids; supports the lower-level linear/exponential and disjoint-mask operations. |
| `perturbations.py` | Owns Python baseline snapshots, constructs native-representation candidates, rebuilds PDF/CDF pairs using native integration and verifies readback/restoration. |
| `runner.py` | Initialises a native model, resets simulation/RNG state per sample, applies transactional mutations, scores batches, restores data and loads/saves deterministic observations. |
| `budget.py` | Shared SQLite history-reservation ledger and a budgeted Runner wrapper. Parallel workers share the stage cap. |
| `storage.py` | Deterministic sample identifiers, atomic per-sample HDF5 observations and JSON serialization helpers. |
| `morris.py` | Valid global trajectory design, elementary effects, bootstrap estimates, statistical resolution and ranking-convergence diagnostics. |
| `analysis.py` | Executes low-fidelity Morris, timing observations and spawned workers; writes screening tables and paired-noise diagnostics. |
| `family_workflow.py` | Selects top inputs per family, executes high-fidelity designs serially/in parallel, fits saved observations and writes validation/Sobol/distribution artifacts. |
| `pce.py` | Orthonormal Legendre basis construction, regression, cross-validation, predictions and portable PCE persistence. |
| `openturns_pce.py` | Optional OpenTURNS regression adapter with matching polynomial normalization, direct native checks and portable coefficient output. |
| `validation.py` | Held-out fit/error/noise metrics and eligibility checks, including defined handling of constant/zero responses. |
| `sobol.py` | First-, total-, pair- and group-support coefficient decompositions plus component-resolution diagnostics. The standard report focuses on first/total and families. |
| `native_plots.py` | Baseline/configured-endpoint native XS, energy and angle plots, actual readback checks, PDF normalization and representation diagnostics. |
| `two_pass.py` | Complete paginated Morris plots and the Morris/Sobol Markdown plot index; the filename is historical. |
| `fingerprint.py` | SHA fingerprint of shared analysis source for provenance and stale-verification checks. |
| `tests/` | Regression tests for settings propagation, model construction, Morris/PCE/Sobol mathematics, optional OpenTURNS and validation edge cases. |

Native C++ transport, C bindings and Python library mutation live in the existing
OpenMC implementation (`src/`, `include/openmc/`, `openmc/lib/`). They are a
required dependency, not duplicated in this suite. The plotted nuclear-data
PDFs/XS come through loaded runtime objects and native APIs. Reading reference
ENDF metadata for yield verification is separate from obtaining plotted angular
or emission-energy data.

## 9. Running the suite

### Environment

From the repository, using the environment configured for this modified OpenMC:

```bash
cd /sandbox/openmc
export OPENMC_CROSS_SECTIONS=/nuclear_data/endfb-viii.0-hdf5/cross_sections.xml
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

Use `/opt/openmc-venv/bin/python` in the following examples. If `python` in your
shell is already that environment, it is equivalent. The Python OpenMC package
and shared library must both contain the required runtime perturbation APIs.
Dependencies include NumPy, SciPy, pandas, Matplotlib and h5py, plus OpenMC's
normal Python dependencies. OpenTURNS is optional:

```bash
/opt/openmc-venv/bin/python -m pip install openturns
```

### Inspect cost without transport

```bash
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --dry-run --output sandbox/studies/Li_analysis/results_new
```

Discovery still initializes native data and can write provenance/cost metadata.
A fresh study may need reference ENDF files; the supported verifier can download
a bounded archive prefix and checks pinned integrity hashes. Dry-run is not an
OpenMC transport calculation.

### Run or resume

```bash
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --run --output sandbox/studies/Li_analysis/results_new
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --resume --output sandbox/studies/Li_analysis/results_new
```

`--run` and `--resume` use the same deterministic-checkpoint execution path.
Successful matching observations are reused. Failed or incompatible observations
are not treated as successful responses. Failed/reserved attempts can consume
a conservative budget reservation, so investigate failures rather than removing
the ledger to circumvent the configured cap.

For FLiBe, replace `Li_analysis` with `FLiBe_analysis` and choose a fresh output
if physical inputs/settings have changed. From another working directory:

```bash
/opt/openmc-venv/bin/python /sandbox/openmc/sandbox/analysis_codes/run.py /sandbox/openmc/sandbox/studies/Li_analysis --dry-run --output /sandbox/openmc/sandbox/studies/Li_analysis/results_new
```

### Refit/replot saved observations

```bash
/opt/openmc-venv/bin/python -m sandbox.analysis_codes Li_analysis --reanalyse --output sandbox/studies/Li_analysis/results_new
```

This performs no transport. Change `PCE_BACKEND` or fitting settings in config
and use this mode to analyse the same observations. It regenerates screening
figures and PCE/Sobol/distribution/report outputs; it does not recreate missing
transport observations or every native nuclear-data gallery. The saved range
must agree with the selected parameter registry. Do not change material/source
settings and then describe old data as new physical calculations.

The latest analysis is shown in `high/`; backend-specific copies are stored under
`pce_backends/numpy/` and `pce_backends/openturns/`. A later analysis with the same
backend can replace its previous derived artifacts. Copy that archive first if
you want to retain multiple hyperparameter variants.

### Create a new study

```bash
mkdir /sandbox/openmc/blanket_analysis
cp /sandbox/openmc/sandbox/studies/FLiBe_analysis/config.py /sandbox/openmc/sandbox/studies/FLiBe_analysis/model.py /sandbox/openmc/sandbox/studies/blanket_analysis/
```

Edit those copies. Do not copy the old results as the new study's observations.
The public model builder is:

```python
import openmc
from . import config as c

def build_model(particles=1000, batches=10):
    # Construct fresh materials, geometry and source here.
    model = openmc.Model(...)
    model.settings.run_mode = 'fixed source'
    model.settings.particles = particles
    model.settings.batches = batches
    # Assign source, temperature and other intended physical settings.
    return model
```

The ellipsis marks user-supplied model construction, not runnable geometry.
Return a fresh model each time. Do not run transport or initialize the native
library during module import. Preserve particle/batch arguments because the
engine supplies different settings for Morris, high-fidelity and verification.

A response requesting production channels from other nuclides requires a verified response definition.
Do not add arbitrary production MTs to the built-in mapping. Instead supply
`responses(model)` in model.py returning `dict[str, analysis_codes.Response]`
and `verify_response(output_directory)` for controlled verification. The latter
must write `verification.json` with `status='passed'` and the current
`analysis_codes.fingerprint.implementation_sha()` in `implementation_sha256`.
The standard analysis expects the scalar `total` response first. This callback
is a verification contract; a passing JSON file must reflect real controlled tests.

## 10. Outputs and where to look

```text
results_new/
├── analysis_request.json          Requested config/model/engine identity
├── analysis_plan.json             Planned counts and history budget
├── dry_run_plan.json              Cost report when dry-run was used
├── provenance/                    Mapping, data hashes and reference evidence
├── verification/                  Controlled response/restoration tests
├── reaction_registry.csv          Complete loaded reaction classification
├── retained_cross_sections.csv    Source-retained XS catalogue
├── retained_parameters.csv/.json  Eligible parameters, IDs, modes and bounds
├── removed_reactions.csv          Reactions above the source-energy maximum
├── removed_parameters.csv         Associated excluded perturbation modes
├── distribution_registry.csv      Native secondary representations/status
├── all_retained_plots/            Baseline data gallery and index.md
├── selected_parameters.csv        Per-family selection and reserved checks
├── parameter_decisions.csv        Full screening/selection decisions
├── low/
│   ├── morris_results.csv         Every analysed parameter and response
│   ├── morris_arrays.npz          Design, means, errors and paired batches
│   ├── morris_convergence.json    Ranking/precision convergence diagnostics
│   ├── paired_noise_diagnostics.json
│   ├── checkpoints/               Native observations and budget ledger
│   └── morris*.png/.pdf           Complete/paginated figures
├── high/
│   ├── pce_training.csv           Parameters, seeds, histories, means and SEs
│   ├── pce_validation.csv         Independent validation observations
│   ├── analysis_arrays.npz        Designs and response/error arrays
│   ├── pce_model.npz              Polynomial powers, coefficients/covariance
│   ├── pce_validation.json        CV choice and held-out diagnostics
│   ├── sobol_indices.csv          All selected first/total indices and flags
│   ├── sobol_families.csv          Support-based grouped family indices
│   ├── tritium_distribution.json  Accepted response-distribution summary
│   ├── checkpoints/               High-fidelity native observations/ledger
│   └── *.png/.pdf                 Validation, Sobol and distribution plots
├── nuclear_data_plots/            Selected native endpoint perturbation gallery
├── pce_backends/                  Backend-specific derived-analysis copies
├── sensitivity_plots.md           Linked Morris/Sobol plot gallery
└── summary.md                    Main concise study report
```

Some files depend on mode and success: a response-distribution figure/JSON is
not retained as a claimed valid result if the total PCE is rejected. A paired
noise file requires paired Morris. Older adopted studies can contain additional
historical files. Actual numeric values and configuration snapshots, not a fixed
parameter count in documentation, determine what was analysed.

Start with `summary.md`, then `sensitivity_plots.md`. Use `removed_reactions.csv`
for threshold exclusions and `distribution_registry.csv` for unsupported laws.
Use `pce_validation.json` before interpreting Sobol bars or a probability plot.
An asterisk/unresolved flag indicates insufficient component resolution, not
necessarily zero physical effect.

## 11. Costs, parallel execution and native restoration

The nominal low-stage histories are

```text
[4 + MORRIS_TRAJECTORIES*(d+1)] * SCREEN_PARTICLES * SCREEN_BATCHES.
```

High-stage histories are

```text
[1 + TRAINING_SAMPLES + VALIDATION_SAMPLES] * HIGH_PARTICLES * HIGH_BATCHES.
```

Add controlled verification cases at their configured particles/batches. These
counts are histories, not independent PCE coefficients or cheap surrogate samples.
Native retries/extra checks may consume additional reserved histories. Inspect
stage execution records and budget ledgers for actual consumption.

`WORKERS` and `HIGH_WORKERS` change concurrency, not the number of histories per
sample. Each process owns a native OpenMC library; the code spawns processes and
does not fork/share writable initialized native state. `TRANSPORT_THREADS` applies
inside a worker. Reserve approximately workers × threads CPUs for the active
stage and avoid additional BLAS oversubscription.

Within an evaluation, the runner resets simulation and tally state, sets the
recorded seed, constructs candidates from baseline snapshots, applies them through
the native perturbation context, verifies native readback, runs actual transport,
extracts batch responses and restores/verifies the baseline. A restoration failure
invalidates the runner rather than allowing contaminated subsequent samples.
No nuclear-data HDF5 library is rewritten for each trial.

## 12. Provenance, compatibility and limitations

Parameter identities include nuclide, MT, family, native product/distribution,
mode, incident interval, coefficient bounds and semantics. Checkpoints also carry
model and native-library provenance, settings and seeds. A changed model, transport
configuration or engine fingerprint can require a new output directory. Do not
force compatibility by deleting guard files or modifying cached response values.
Reanalysis is the intended path for changing fitting settings on existing data.

The completed historical FLiBe study has a special adoption path that checks its
original config and archived model builder before preserving/reusing completed
transport datasets. This is not a general migration path for partially completed
old engines. Source snapshots and original transport provenance remain archived.

Main limitations to retain when reporting conclusions:

- Eligible native laws are a subset of inspectable laws; conditional/correlated
  and unsupported analytic mutations are not silently approximated.
- XS mutation currently requires the framework's supported loaded-temperature
  representation; arbitrary multi-temperature transactional parameterization
  is not implemented.
- Strict threshold filtering is appropriate to the supported source/transport
  assumptions; unusual systems with meaningful neutron upscattering or other
  energy-gain mechanisms need a physics review of that screening assumption.
- The current runner's standard high-level figures focus on scalar total and
  constituent production; arbitrary fine spatial PCE reporting is not supplied
  merely by adding an OpenMC mesh tally.
- Reduced-space sensitivity is conditional; omitted parameters can still matter.
- Good total-response fit does not establish reliable weak-component rankings.
- Monte Carlo noise, PCE truncation and hypothetical input variability are
  distinct uncertainties and must not be conflated.
- Changing the source direction requires a compatible source position; no automatic
  ray/sphere intersection placement is inferred for an arbitrary geometry.

## 13. Tests and troubleshooting

```bash
cd /sandbox/openmc
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/openmc-venv/bin/python -m pytest analysis_codes/tests -q
```

The recent configuration work passed 30 automated tests and a bounded native
end-to-end check with nondefault Morris settings, ±12% perturbation coefficients,
two high-fidelity workers and exact restoration. These tests validate machinery,
not convergence of a large scientific campaign. Saved check logs are under
`Li_analysis/results/configuration_validation/` where present.

| Symptom | Meaning/action |
|---|---|
| Changed model/settings/engine guard | Use a new --output for new transport, or --reanalyse for compatible saved observations. Do not delete the provenance guard to mix data. |
| Configured limit differs from saved registry | The old transport data represent another input prior; regenerate data in a new output for the new range. |
| Budget exceeded | Review the dry-run counts, histories and stage caps; explicitly choose a sufficient budget or reduce the design. |
| Missing OpenTURNS | Install the optional dependency or select numpy. |
| Zero Li7 production below threshold | Physically closed channels can give zero response/variance; inspect thresholds before interpreting sensitivity. |
| No accepted PCE/distribution | Inspect held-out diagnostics, noise and training adequacy; a successful run is not automatic surrogate validation. |
| Unexpected material density/temperature | Inspect model.py and the exported actual model/provenance. The current supplied builders use the config values directly. |
| Fewer modes than requested per family | There may be fewer source-accessible safely mutable native laws; consult inventories and exclusions. |
| Source misses the target or starts outside the world | Check position, direction, radius, world margin and standoff together. |

For a new study, first check the physical inputs, run dry-run, use a deliberately
bounded test configuration, inspect verification/validation, and only then launch
the explicitly budgeted production campaign.


## Production QoI versus transport nuclides

`PRODUCTION_CHANNELS` defines which verified reaction rates are summed into the
tritium-production QoI. For example, Li6 MT105 plus Li7 MT52–82 measures lithium
production even when lead is present. Only the configured scored channels require
production-yield verification. No exclusion list is needed.

`TARGET_NUCLIDES` independently controls perturbation discovery. With `None`,
all material nuclides, including lead, remain eligible for safely supported
perturbations. Their transport changes can affect the configured lithium QoI.
The configured sum is not a claim about production from unlisted nuclides.
