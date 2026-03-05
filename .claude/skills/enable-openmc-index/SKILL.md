---
name: enable-openmc-index
description: Enable the OpenMC codebase index for this session. Provides semantic code search and structural repo mapping. Run this when investigating, modifying, or debugging OpenMC code.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase index for this session. This gives you two tools:
1. **Semantic search** (`openmc_search.py`) - Find related code across the codebase by concept
2. **Structural map** (`openmc_map.py`) - See condensed code structure of files and their neighbors

## Step 1: Ensure the virtual environment exists

Check if `.claude/cache/.venv/` exists. If not, create it and install dependencies:

```bash
if [ ! -d .claude/cache/.venv ]; then
    python3 -m venv .claude/cache/.venv
    .claude/cache/.venv/bin/pip install -r .claude/tools/requirements.txt -q
    echo "INSTALLED"
else
    echo "VENV_EXISTS"
fi
```

## Step 2: Ensure the RAG index exists

Check if `.claude/cache/rag_index/` exists. If not, build it:

```bash
if [ ! -d .claude/cache/rag_index ]; then
    echo "Building RAG index for the first time (this takes ~3 minutes)..."
    HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
    echo "INDEX_BUILT"
else
    echo "INDEX_EXISTS"
fi
```

Note: The repo map tool (`openmc_map.py`) does NOT need a pre-built index - it generates maps on the fly using tree-sitter. Only the RAG search needs the vector index.

## Step 3: Learn the tool APIs

Run `--help` for both tools to see their full APIs:

```bash
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help
.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help
```

Read and internalize the output so you know all available options.

## Step 4: Confirm activation

Tell the user the OpenMC index is active and briefly describe the two tools:
- **Semantic search**: Find related code by concept (e.g., "particle seed initialization")
- **Structural map**: See condensed code structure around specific files

## How to use the tools after activation

**Typical workflow:**

1. Use `openmc_search.py` to discover which files are relevant to your task
2. Use `openmc_map.py` on those files to understand their structure and neighbors
3. Use Read/Grep to dive into the specific code you need to change

**When to use semantic search** (`openmc_search.py`):
- Investigating how a change might affect other parts of the codebase
- Finding code that does something conceptually similar but with different naming
- Discovering cross-cutting concerns across run modes

**When to use the repo map** (`openmc_map.py`):
- Understanding the structure of unfamiliar files before modifying them
- Seeing what classes/methods neighbor the code you're working on
- Getting a condensed overview of a subsystem (pass multiple files)

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help` (structural map)"
