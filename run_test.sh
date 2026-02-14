#!/bin/bash
# Run a single test by name or path.
#
# Usage:
#   ./run_test.sh <test_id_or_file>          Verbose output
#   ./run_test.sh -q <test_id_or_file>       Quiet: pass/fail + details file only

if [[ -z "$1" ]]; then
    echo "Usage: ./run_test.sh [-q] <test_id_or_file>"
    echo ""
    echo "Examples:"
    echo "  ./run_test.sh tests/regression_tests/random_ray_k_eff/test.py::test_random_ray_basic"
    echo "  ./run_test.sh -q tests/regression_tests/random_ray_k_eff/test.py"
    exit 1
fi

QUIET=false
if [[ "$1" == "-q" ]]; then
    QUIET=true
    shift
fi

if [[ -z "$1" ]]; then
    echo "Error: no test specified"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if $QUIET; then
    OMP_NUM_THREADS=8 python -m pytest "$1" -v --tb=short > /tmp/openmc_run_test.txt 2>&1
    EXIT_CODE=$?
    if [[ $EXIT_CODE -eq 0 ]]; then
        grep -E "passed" /tmp/openmc_run_test.txt | tail -1
    else
        grep -E "(FAILED|ERROR|failed|error)" /tmp/openmc_run_test.txt | head -5
        echo ""
        echo "Details: /tmp/openmc_run_test.txt"
    fi
    exit $EXIT_CODE
else
    OMP_NUM_THREADS=8 python -m pytest "$1" -v --tb=short 2>&1
fi
