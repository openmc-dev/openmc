#!/bin/bash
# Clean build OpenMC from scratch.
# Usage: ./build.sh [options]
#   --incremental  Skip cleaning and cmake (faster rebuilds)
#   -q             Quiet mode: only show errors

set -e

BUILD_DIR="$HOME/openmc/build"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Parse arguments
INCREMENTAL=false
QUIET=false
for arg in "$@"; do
    case "$arg" in
        --incremental) INCREMENTAL=true ;;
        -q)            QUIET=true ;;
    esac
done
mkdir -p $BUILD_DIR

if $QUIET; then
    LOG="/tmp/openmc_build.txt"

    if ! $INCREMENTAL; then
        rm -rf "$BUILD_DIR"/*
        cd "$BUILD_DIR"
        cmake \
            -DOPENMC_USE_MPI=off \
            -DCMAKE_INSTALL_PREFIX=./install \
            "$SCRIPT_DIR" > "$LOG" 2>&1
    fi

    cd "$BUILD_DIR"
    if make -j10 install >> "$LOG" 2>&1; then
        echo "=== Build succeeded ==="
    else
        echo "=== Build FAILED ==="
        echo ""
        # Show just the error lines
        grep -iE "(error:|fatal|undefined|cannot find)" "$LOG" | head -20
        echo ""
        echo "Full log: $LOG"
        exit 1
    fi
else
    if ! $INCREMENTAL; then
        echo "=== Cleaning build directory ==="
        rm -rf "$BUILD_DIR"/*
    fi

    cd "$BUILD_DIR"

    if ! $INCREMENTAL; then
        echo "=== Running cmake ==="
        cmake \
            -DOPENMC_USE_MPI=off \
            -DCMAKE_INSTALL_PREFIX=./install \
            "$SCRIPT_DIR" 2>&1
    fi

    echo "=== Building ==="
    make -j10 install 2>&1

    echo "=== Build complete ==="
fi
