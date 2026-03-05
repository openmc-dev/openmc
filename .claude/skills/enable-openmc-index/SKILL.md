---
name: enable-openmc-index
description: Enable the OpenMC codebase tools for this session. Provides semantic code search, structural repo mapping, and LSP-based C++ code navigation.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase tools for this session:

1. **`openmc_search.py`** — RAG semantic search. Embeds your query and searches a vector index. Good for finding conceptually related code even when naming differs. Covers C++, Python, and docs.
2. **`openmc_map.py`** — Structural repo map via aider/tree-sitter. Given focus files you already have open, shows condensed code skeletons of their *neighbors* — other files that share symbols with them. Useful for seeing surrounding context. Caveat: matches identifiers by name only, so common names like `push_back` or `__init__` create false connections in the neighbor ranking.
3. **`openmc_lsp.py`** — LSP navigation via clangd. Compiler-accurate `definition`, `references`, and `related` commands for C++. Zero false edges. Requires clangd and compile_commands.json.

## Step 1: Ensure the virtual environment exists

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

The semantic search tool needs a pre-built vector index. The other two tools work without it.

```bash
if [ ! -d .claude/cache/rag_index ]; then
    echo "Building RAG index for the first time (this takes ~3 minutes)..."
    HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
    echo "INDEX_BUILT"
else
    echo "INDEX_EXISTS"
fi
```

## Step 3: Learn the tool APIs

Run `--help` for each tool to learn their full APIs:

```bash
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help
.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help
.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help
```

Read and internalize the output so you know all available options.

## Step 4: Confirm activation

Tell the user the tools are active. Mention if clangd or compile_commands.json
are missing (the LSP tool will report this itself if invoked without them).

## When to use each tool

- **`openmc_search.py`**: "What code is conceptually related to X?" — best for broad discovery, cross-cutting concerns, searching docs and Python code
- **`openmc_lsp.py`**: "Where is this C++ symbol defined, who calls it, and what files are connected to this one?" — returns file:line locations with compiler accuracy, zero false positives
- **`openmc_map.py`**: "What other code surrounds the files I'm working on?" — shows condensed signatures of neighboring files. Be aware its neighbor ranking has false edges from common method names

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help` (structural map)
- `.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help` (LSP navigation, C++ only)"
