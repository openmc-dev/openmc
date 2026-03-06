#!/bin/bash
# Bootstrap the Python venv (if needed) and start the OpenMC MCP server.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CACHE_DIR="$(dirname "$SCRIPT_DIR")/cache"
VENV_DIR="$CACHE_DIR/.venv"

if [ ! -d "$VENV_DIR" ]; then
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install -q -r "$SCRIPT_DIR/requirements.txt"
fi

exec "$VENV_DIR/bin/python" "$SCRIPT_DIR/openmc_mcp_server.py"
