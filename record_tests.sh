#!/bin/bash
# Run the full regression test suite and record which tests pass.
# Saves passing test IDs to passing_tests.txt (one per line).
# Usage: ./record_tests.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_FILE="$SCRIPT_DIR/passing_tests.txt"
RAW_OUTPUT="/tmp/openmc_test_output.txt"

echo "=== Running full regression test suite ==="
echo "    (This takes ~3-5 minutes)"

cd "$SCRIPT_DIR"
OMP_NUM_THREADS=8 python -m pytest tests/regression_tests/ \
    --ignore=tests/regression_tests/cmfd_feed_ng \
    -v --tb=no 2>&1 | tee "$RAW_OUTPUT"

echo ""
echo "=== Extracting passing tests ==="

# Extract test IDs that PASSED from verbose output
grep " PASSED" "$RAW_OUTPUT" | sed 's/ PASSED.*//' | sed 's/^[[:space:]]*//' > "$OUTPUT_FILE"

TOTAL=$(wc -l < "$OUTPUT_FILE" | tr -d ' ')
echo "Recorded $TOTAL passing tests to: $OUTPUT_FILE"
echo ""

# Also record metadata at the top as comments
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
BRANCH=$(cd "$SCRIPT_DIR" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
COMMIT=$(cd "$SCRIPT_DIR" && git rev-parse --short HEAD 2>/dev/null || echo "unknown")

# Prepend metadata (use temp file to avoid clobbering)
TMPFILE=$(mktemp)
cat > "$TMPFILE" <<EOF
# OpenMC passing tests baseline
# Recorded: $TIMESTAMP
# Branch: $BRANCH
# Commit: $COMMIT
# Total passing: $TOTAL
EOF
cat "$OUTPUT_FILE" >> "$TMPFILE"
mv "$TMPFILE" "$OUTPUT_FILE"

echo "Done. Summary:"
tail -1 "$RAW_OUTPUT"
