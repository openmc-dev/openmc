---
name: enable-openmc-index
description: Enable the OpenMC codebase index for this session. Provides a structural repo map and semantic code search. Run this when investigating, modifying, or debugging OpenMC code.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase index for this session. This gives you:
1. A structural repo map showing the most important classes, functions, and their relationships
2. Semantic search across all source code, tests, and documentation

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

## Step 2: Ensure the index exists

Check if `.claude/cache/rag_index/` exists. If not, build it:

```bash
if [ ! -d .claude/cache/rag_index ]; then
    echo "Building index for the first time (this takes ~3 minutes)..."
    .claude/cache/.venv/bin/python .claude/tools/repomap/generate_repomap.py
    HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
    echo "INDEX_BUILT"
else
    echo "INDEX_EXISTS"
fi
```

## Step 3: Load the repo map

Read the file `.claude/cache/repomap.md` and internalize the codebase structure.

## Step 4: Confirm activation

Tell the user the OpenMC index is active and briefly describe what's available:
- The repo map is loaded (structural overview)
- Semantic search is ready via: `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py "query"`

## Using semantic search after activation

For the rest of this session, before modifying unfamiliar code or when investigating how a change might affect other parts of the codebase, search for related code:

```bash
# Search source code
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py "your query here"

# Search documentation
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py "your query" --docs

# Search both code and docs
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py "your query" --all

# Find code related to a specific file
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --related src/somefile.cpp
```

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:
"The OpenMC search index is available. Read .claude/cache/repomap.md for a structural overview. Use `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py 'query'` for semantic search."
