#!/usr/bin/env python3
"""MCP server providing OpenMC code navigation tools.

Exposes two tools:
  - openmc_rag_search:    Semantic search across the codebase and docs
  - openmc_rag_rebuild:   Rebuild the RAG vector index
"""

from mcp.server.fastmcp import FastMCP
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# Suppress noisy logging from httpx and huggingface_hub before any imports
# that trigger HTTP requests.  The MCP transport uses stderr, so stray log
# lines there would corrupt the protocol stream.
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("huggingface_hub").setLevel(logging.ERROR)
logging.getLogger("sentence_transformers").setLevel(logging.WARNING)


OPENMC_ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = OPENMC_ROOT / ".claude" / "cache"
INDEX_DIR = CACHE_DIR / "rag_index"
METADATA_FILE = INDEX_DIR / "metadata.json"

# Add tool subdirectories to path for imports
TOOLS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS_DIR / "rag"))

mcp = FastMCP("openmc-code-tools")

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
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
