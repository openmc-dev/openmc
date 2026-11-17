# Runtime nuclear-data perturbation audit

This note records the runtime mutation interfaces in this fork after the
project-specific analysis work was relocated under `sandbox/`. The OpenMC
runtime implementation remains in `openmc/lib/`, `include/openmc/` and
`src/`; regression tests remain under `tests/`.

## Supported edits

| Data family | Public Python interface | Native representation currently editable | Important constraints |
|---|---|---|---|
| Microscopic cross section | `openmc.lib.get_nuclide_xs`, `set_nuclide_xs`, `perturb_nuclide_xs`; transactional `nuclear_data_perturbation` | Selected non-redundant reaction partial cross section for supported nuclides/temperatures | Input uses the native incident-energy grid. The Python context retains the native threshold and zero prefix. Derived totals must be rebuilt after direct low-level writes. |
| Outgoing energy | `get_emission_distribution`, `set_emission_distribution`, native `sample_emission_energy` | `UncorrelatedAngleEnergy` with a `ContinuousTabular` energy component | Incident grid and table support/interpolation are retained; PDF and CDF are validated and rebuilt together. Correlated/analytic laws are read-only. |
| Angular | `get_angular_distribution`, `set_angular_distribution`, native `sample_angular_distribution` | Existing tabular-angle component of an uncorrelated law | Exact incident and cosine grids and interpolation are retained; cosine bounds are [-1,1]; PDF is nonnegative and normalized and CDF/integrals rebuilt. Isotropic laws are inspectable but not converted to tabular. |

The transactional context is implemented in
[`openmc/lib/perturbation.py`](../openmc/lib/perturbation.py). It snapshots
native data before applying edits, validates before/after mutation, rebuilds
derived cross sections, and restores backups in reverse order on normal exit
or exceptions. It refuses edits during a transport batch. Snapshots exposed to
Python are copies; callers do not receive writable pointers into native data.

## Trace through native code

- Python C bindings and cross-section getters/setters are in
  [`openmc/lib/core.py`](../openmc/lib/core.py) (`nuclide_rebuild_derived_xs`,
  `get_nuclide_xs`, and associated setters).
- Secondary-law Python inspection, mutation and sampling wrappers are in
  [`openmc/lib/secondary.py`](../openmc/lib/secondary.py).
- Shared native secondary resolver, validation, snapshot, setters and native
  sampling entry points are in [`src/secondary_data.cpp`](../src/secondary_data.cpp)
  and [`include/openmc/secondary_data.h`](../include/openmc/secondary_data.h).
- `Nuclide::rebuild_derived_xs` and `validate_derived_xs` are in
  [`src/nuclide.cpp`](../src/nuclide.cpp). They rebuild/check the cached total,
  absorption, fission, nu-fission and photon-production arrays from the
  non-redundant partial reactions for each temperature.
- Native continuous-tabular energy edits and PDF/CDF construction are in
  [`src/distribution_energy.cpp`](../src/distribution_energy.cpp).
- Batch-state guards and simulation reset behavior are in
  [`src/simulation.cpp`](../src/simulation.cpp).

Reaction-rate scoring and transport sample the reaction partial cross
sections, so a supported partial edit affects reaction selection/rates after
the derived arrays are rebuilt. A direct low-level partial cross-section
setter must be followed by `openmc.lib.nuclide_rebuild_derived_xs(index)`;
the high-level perturbation context does that itself. Fission, unresolved
resonance, multipole and redundant reaction data are not general-purpose
targets of the unified mutation context because they have dependent native
representations/caches.

### Secondary-law sampling semantics

For supported uncorrelated laws, energy and angle are sampled independently
from their respective components. `CorrelatedAngleEnergy` and `KalbachMann`
encode a joint/conditional energy-angle law; the independent angular sampler
rejects them, and this mutation API does not pretend their angle and energy
are independent. Analytic energy laws, correlated laws, and other unsupported
representations remain inspectable where implemented but reject mutation
explicitly. The API does not convert isotropic or correlated data to tabular
data.

## Lifecycle and restoration

Mutations require an initialized library and are blocked while OpenMC is in a
batch. The Python perturbation context initializes/resets simulation tally
results before a trial; repeated trial runs can therefore use the same
initialized library session with fresh tally results. It captures all
requested backups before applying candidates, and its `finally` path restores
already-applied edits even if a later edit fails. Cross-section restoration
rebuilds derived arrays after restoring the baseline partials. Secondary
restoration writes back deep native snapshots. Existing lifecycle tests cover
repeated edits, in-batch rejection, exception restoration, sampler response,
and baseline recovery.

One known boundary is cached MT205: this is preprocessed production data and
is not reconstructed by `rebuild_derived_xs`. UQ responses for runtime-edited
partial reactions must be tallied from verified reaction-specific channels,
not from cached MT205. See the separate `mt205_investigation/` report and the
corrected response workflow under `sandbox/analysis_codes/`.

## Verification performed

Commands run from the repository root:

```bash
cmake --build build -j2
ctest --test-dir build -R 'secondary_data|emission_spectrum' --output-on-failure
OMP_NUM_THREADS=2 python -m pytest -q \
  tests/unit_tests/test_runtime_uq.py \
  tests/unit_tests/test_lib_emission_spectrum.py \
  tests/unit_tests/test_lib_angular.py \
  tests/unit_tests/test_lib_angular_lifecycle.py \
  tests/unit_tests/test_lib_secondary.py \
  tests/unit_tests/test_lib_secondary_chain.py \
  tests/unit_tests/test_lib_secondary_watt.py
OMP_NUM_THREADS=2 python -m pytest -q \
  sandbox/analysis_codes/tests sandbox/tokomak_analysis/tests
```

Build passed; C++ tests passed 2/2; focused Python runtime tests passed
130/130; relocated analysis/tokamak tests passed 65/65. Relocated smoke checks
also exercised the Li-6 MT2 angular sampler (5,000 native samples; tabular
normalization 0.999999953) and the emission-energy transport example (the
perturbed run's sampled outgoing-energy mean differed from the control as
expected). A relocated tokamak sensitivity dry-run completed and resolved
the study/model paths successfully.

The runtime regression command

```bash
OMP_NUM_THREADS=2 python -m pytest -q \
  tests/regression_tests/runtime_nuclear_data_uq \
  tests/regression_tests/runtime_emission_spectrum
```

reported 3 passed and 1 reference-data mismatch. The configured
`OPENMC_CROSS_SECTIONS` points to ENDF/B-VIII.0 HDF5, while the failing
emission-spectrum regression reference was generated against the official
NNDC HDF5 dataset. The test's runtime transport assertion (same first batch,
then a changed outgoing-energy source bank) completed; only comparison with
the incompatible stored reference remains unverified. No data were downloaded
or changed to mask this mismatch.

## Relocated analysis layout and Git handling

Project-specific scripts, models, configs, analysis packages and report
sources now live under `sandbox/`. OpenMC's own examples and its unit and
regression test trees remain in their conventional locations. Analysis source
is intentionally not ignored; generated campaign results, local build
outputs, checkpoints, local cross-section indexes and nuclear-data files are
ignored by `.gitignore`. Git status shows these moved/new sandbox source files
as untracked because nothing was staged; they are available for review and can
be staged by the repository owner. No commits, pushes or remote changes were
made.
