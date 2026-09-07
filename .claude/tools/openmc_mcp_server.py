#!/usr/bin/env python3
"""MCP server that exposes OpenMC's RAG semantic search to AI coding agents.

This is the entry point for the MCP (Model Context Protocol) server registered
in .mcp.json at the repo root. When an MCP-capable agent (e.g. Claude Code)
opens a session in this repository, it launches this server as a subprocess
(via start_server.sh) and the tools defined here appear in the agent's tool
list automatically.

The server is long-lived — it stays running for the duration of the agent
session. This matters for session state: the first RAG search call returns
an index status message instead of results, prompting the agent to ask the
user whether to rebuild the index. That first-call flag resets each session.

Tools exposed:
  openmc_rag_search  — semantic search across the codebase and docs
  openmc_rag_rebuild — rebuild the RAG vector index

The actual search/indexing logic lives in the rag/ subdirectory (openmc_search.py,
indexer.py, chunker.py, embeddings.py). This file is just the MCP interface
layer and session state management.
"""

from mcp.server.fastmcp import FastMCP
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# MCP communicates over stdin/stdout with JSON-RPC framing. Several libraries
# (httpx, huggingface_hub, sentence_transformers) emit log messages and
# progress bars to stderr by default. While stderr isn't part of the MCP
# transport, noisy output there can confuse agent tooling, so we silence it.
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("huggingface_hub").setLevel(logging.ERROR)
logging.getLogger("sentence_transformers").setLevel(logging.WARNING)

# Path constants. This file lives at .claude/tools/openmc_mcp_server.py,
# so parents[2] is the OpenMC repo root.
OPENMC_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = OPENMC_ROOT / ".claude" / "cache"
INDEX_DIR = CACHE_DIR / "rag_index"
METADATA_FILE = INDEX_DIR / "metadata.json"

# The RAG modules (openmc_search, indexer, etc.) live in .claude/tools/rag/.
# We add that directory to sys.path so we can import them directly.
TOOLS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS_DIR / "rag"))

mcp = FastMCP("openmc-code-tools")

# First-call flag: the first openmc_rag_search call of each session returns
# index status info instead of search results, so the agent can ask the user
# whether to rebuild. This resets when the server process restarts (i.e. each
# new agent session).
_rag_first_call = True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_current_branch():
    """Get the current git branch name."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, cwd=str(OPENMC_ROOT),
        )
        if result.returncode != 0 or not result.stdout.strip():
            return "unknown"
        return result.stdout.strip()
    except Exception:
        return "unknown"


def _get_index_metadata():
    """Read index build metadata, or None if unavailable."""
    if not METADATA_FILE.exists():
        return None
    try:
        return json.loads(METADATA_FILE.read_text())
    except Exception:
        return None


def _save_index_metadata():
    """Save index build metadata alongside the index."""
    metadata = {
        "built_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "branch": _get_current_branch(),
    }
    METADATA_FILE.write_text(json.dumps(metadata, indent=2))


def _check_index_first_call():
    """On the first RAG call of the session, return a status message for the
    agent to relay to the user.  Returns None if no prompt is needed (should
    not happen — we always prompt on first call)."""
    current_branch = _get_current_branch()

    if not INDEX_DIR.exists():
        return (
            "No RAG index found. Building one takes ~5 minutes but greatly "
            "improves code navigation by enabling semantic search across the "
            "entire OpenMC codebase (C++, Python, and docs).\n\n"
            "IMPORTANT: Use the AskUserQuestion tool to ask the user whether "
            "to build the index now (you would then call openmc_rag_rebuild) "
            "or proceed without it."
        )

    meta = _get_index_metadata()
    if meta:
        built_at = meta.get("built_at", "unknown time")
        built_branch = meta.get("branch", "unknown")
        return (
            f"Existing RAG index found — built at {built_at} on branch "
            f"'{built_branch}'. Current branch is '{current_branch}'.\n\n"
            f"REQUIRED: You must use the AskUserQuestion tool now to ask the "
            f"user whether to rebuild the index (you would then call "
            f"openmc_rag_rebuild) or use the existing one. Do not skip this "
            f"step — the user may have uncommitted changes. Do not decide "
            f"on their behalf."
        )

    return (
        f"RAG index found but has no build metadata. "
        f"Current branch is '{current_branch}'.\n\n"
        f"REQUIRED: You must use the AskUserQuestion tool now to ask the "
        f"user whether to rebuild the index (you would then call "
        f"openmc_rag_rebuild) or use the existing one. Do not skip this "
        f"step. Do not decide on their behalf."
    )


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@mcp.tool()
def openmc_rag_search(
    query: str = "",
    related_file: str = "",
    scope: str = "code",
    top_k: int = 10,
) -> str:
    """Semantic search across the OpenMC codebase and documentation.

    Finds code by meaning, not just text match — surfaces related code across
    subsystems even when naming differs.  Use for discovery and exploration
    before reaching for grep.  Covers C++, Python, and RST docs.

    Args:
        query: Search query (e.g. "particle weight adjustment variance reduction")
        related_file: Instead of a text query, find code related to this file
        scope: "code" (default), "docs", or "all"
        top_k: Number of results to return (default 10)
    """
    global _rag_first_call

    # First call of the session — prompt the agent to check with the user
    if _rag_first_call:
        _rag_first_call = False
        status = _check_index_first_call()
        if status:
            return status

    # No index available
    if not INDEX_DIR.exists():
        return (
            "No RAG index available. Call openmc_rag_rebuild() to build one "
            "(takes ~5 minutes)."
        )

    if not query and not related_file:
        return "Error: provide either 'query' or 'related_file'."

    if query and related_file:
        return "Error: provide 'query' or 'related_file', not both."

    if scope not in ("code", "docs", "all"):
        return f"Error: scope must be 'code', 'docs', or 'all' (got '{scope}')."

    if top_k < 1:
        return f"Error: top_k must be at least 1 (got {top_k})."

    try:
        from openmc_search import (
            get_db_and_embedder, search_table, format_results, search_related,
        )

        db, embedder = get_db_and_embedder()

        if related_file:
            results = search_related(db, embedder, related_file, top_k)
            return format_results(results, f"Code related to {related_file}")
        elif scope == "all":
            code_results = search_table(db, embedder, "code", query, top_k)
            doc_results = search_table(db, embedder, "docs", query, top_k)
            return (format_results(code_results, "Code") + "\n"
                    + format_results(doc_results, "Documentation"))
        elif scope == "docs":
            results = search_table(db, embedder, "docs", query, top_k)
            return format_results(results, "Documentation")
        else:
            results = search_table(db, embedder, "code", query, top_k)
            return format_results(results, "Code")
    except Exception as e:
        return f"Error during search: {e}"


@mcp.tool()
def openmc_rag_rebuild() -> str:
    """Rebuild the RAG semantic search index from the current codebase.

    Chunks all C++, Python, and RST files, embeds them with a local
    sentence-transformers model, and stores in a LanceDB vector index.
    Takes ~5 minutes on 10 CPU cores.  Call this after pulling new code
    or switching branches.
    """
    global _rag_first_call
    _rag_first_call = False  # no need to prompt after an explicit rebuild

    try:
        import io
        from indexer import build_index

        old_stdout = sys.stdout
        sys.stdout = captured = io.StringIO()
        try:
            build_index()
        finally:
            sys.stdout = old_stdout

        _save_index_metadata()

        branch = _get_current_branch()
        build_output = captured.getvalue()
        return (
            f"Index rebuilt successfully on branch '{branch}'.\n\n"
            f"{build_output}"
        )
    except Exception as e:
        return f"Error rebuilding index: {e}"


if __name__ == "__main__":
    mcp.run()
