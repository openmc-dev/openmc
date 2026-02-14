# OpenMC Development Environment

Read `agent_build_and_testing_workflow.md` in this directory before starting work.
It contains the full build/test workflow, script usage, and conventions.

## Key Commands
- `./build.sh --incremental -q` — rebuild after C++ changes (quiet, errors in `/tmp/openmc_build.txt`)
- `./check_tests.sh -q --smoke` — quick regression check (~10% of tests)
- `./check_tests.sh -q` — full regression check (details in `/tmp/openmc_regression.txt`)
- `./run_test.sh -q <test_id>` — run a single test (details in `/tmp/openmc_run_test.txt`)
- `./record_tests.sh` — re-record baseline after intentional changes

## Workflow
1. Edit C++ or Python code
2. `./build.sh --incremental -q` — rebuild
3. `./check_tests.sh -q --smoke` — quick sanity check
4. `./check_tests.sh -q` — full regression check before committing

## Remotes
- `origin` — fork (git@github.com:jtramm/openmc.git) — push here
- `upstream` — official (https://github.com/openmc-dev/openmc.git) — pull from here

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work. Info in the
`agent_build_and_testing_workflow.md` supercedes anything in AGENTS.md.
