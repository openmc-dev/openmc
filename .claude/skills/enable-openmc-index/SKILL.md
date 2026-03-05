---
name: enable-openmc-index
description: Enable the OpenMC codebase tools for this session. Provides semantic code search, structural repo mapping, and LSP-based C++ code navigation.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase tools for this session:

1. **Semantic search** (`openmc_search.py`) - Find related code by concept across C++, Python, and docs
2. **Structural map** (`openmc_map.py`) - Condensed code structure of files and their neighbors (C++ and Python)
3. **LSP navigation** (`openmc_lsp.py`) - Compiler-accurate definition, references, and related-file discovery (C++ only, requires clangd and compile_commands.json)

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

- **`openmc_search.py`**: Finding code by concept, discovering cross-cutting concerns, searching docs
- **`openmc_lsp.py`**: Go-to-definition, find-references, discovering which C++ files are connected by real typed references
- **`openmc_map.py`**: Seeing condensed code structure and class/function signatures of files you're about to modify

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help` (structural map)
- `.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help` (LSP navigation, C++ only)"
