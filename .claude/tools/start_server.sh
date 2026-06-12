#!/bin/bash
# Bootstrap the Python venv (if needed) and start the OpenMC MCP server.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CACHE_DIR="$(dirname "$SCRIPT_DIR")/cache"
VENV_DIR="$CACHE_DIR/.venv"
SENTINEL="$VENV_DIR/.installed"

if ! command -v python3 >/dev/null 2>&1; then
    echo "Error: python3 not found on PATH." >&2
    exit 1
fi

if ! python3 -c 'import sys; assert sys.version_info >= (3,12)' 2>/dev/null; then
    echo "Error: Python 3.12+ is required." >&2
    exit 1
fi

if [ ! -f "$SENTINEL" ]; then
    rm -rf "$VENV_DIR"
    mkdir -p "$CACHE_DIR"
    python3 -m venv "$VENV_DIR"

    if ! "$VENV_DIR/bin/pip" install -q -r "$SCRIPT_DIR/requirements.txt"; then
        echo "Error: pip install failed. Remove $VENV_DIR and retry." >&2
        rm -rf "$VENV_DIR"
        exit 1
    fi

    touch "$SENTINEL"
fi

exec "$VENV_DIR/bin/python" "$SCRIPT_DIR/openmc_mcp_server.py"
