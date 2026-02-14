# OpenMC Build & Testing Workflow for Agents

## Quick Reference

| Task | Command |
|------|---------|
| Clean build | `./build.sh` |
| Incremental build | `./build.sh --incremental` |
| Run all tests, check for regressions | `./check_tests.sh -q` |
| Smoke test (~10% of tests) | `./check_tests.sh -q --smoke` |
| Run a single test (quiet) | `./run_test.sh -q <test_id>` |
| Run a single test (verbose) | `./run_test.sh <test_id>` |
| Record new baseline | `./record_tests.sh` |

## Build

OpenMC is a C++/Python project. C++ changes require recompilation before testing.

```bash
# After modifying C++ files (src/ or include/), rebuild:
./build.sh --incremental

# If cmake configuration changed (CMakeLists.txt), do a clean build:
./build.sh
```

The build uses clang, no MPI, no DAGMC. Output goes to `~/openmc/build/`.

**Do not run `pip install`** — OpenMC is already installed in development mode and the executable/library etc is in the environment.

## Testing

### Checking for regressions

The file `passing_tests.txt` contains the baseline of tests that should pass.
To verify no regressions were introduced:

```bash
# Quiet mode — just progress counters and result (preferred for agents)
./check_tests.sh -q

# Smoke test — runs ~10% of test files, good for quick sanity checks
./check_tests.sh -q --smoke

# Verbose mode — also shows per-test PASSED/FAILED lines
./check_tests.sh
./check_tests.sh --smoke
```

The check script runs 8 test files in parallel (1 thread each) and stops at
the first regression. On failure, it saves full pytest output to
`/tmp/openmc_regression.txt` — read that file for failure details.

### Running a single test

When debugging a specific failure, use `run_test.sh`:

```bash
# Quiet mode (preferred for agents) — just pass/fail, details in /tmp/openmc_run_test.txt
./run_test.sh -q tests/regression_tests/random_ray_k_eff/test.py::test_random_ray_basic

# Verbose mode — full pytest output including OpenMC stdout
./run_test.sh tests/regression_tests/random_ray_k_eff/test.py::test_random_ray_basic

# Run all tests in a file:
./run_test.sh -q tests/regression_tests/random_ray_k_eff/test.py
```

On failure in quiet mode, if needed, read `/tmp/openmc_run_test.txt` for full details.

### Recording a new baseline

After intentionally changing test results (e.g., updating a feature that
changes numerical output), re-record the passing tests:

```bash
./record_tests.sh
```

This runs the full suite and writes `passing_tests.txt`. Other tests outside this file
may be failing in this environment due to differences in compiler/std libraries on this
system vs. the cloud CI (where ALL tests pass).

## Typical Workflow

1. **Edit code** in `src/` or `include/`
2. **Rebuild**: `./build.sh --incremental`
3. **Quick check**: `./check_tests.sh -q --smoke`
4. **Full check**: `./check_tests.sh -q`
5. **If a test regresses**: read `/tmp/openmc_regression.txt`, fix the issue, go to step 2
6. **If test results intentionally changed**: update reference files, run `./record_tests.sh`

## Notes

- Some tests have intra-file dependencies (e.g., `test_full` generates data
  that `test_depletion_results_to_material` reads). The check script handles
  this by running all tests in each file, not just the passing ones.
- There are sometimes test failures on the `develop` branch that
  are unrelated to our work (due to environmental differences on this machine vs. the cloud CI). These are excluded from `passing_tests.txt`.
