---
name: enable-openmc-index
description: Enable the OpenMC codebase tools for this session. Provides semantic code search, structural repo mapping, and LSP-based C++ code navigation.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase tools for this session:

1. **`openmc_search.py`** — RAG semantic search. Chunks code at function/class boundaries, embeds with sentence-transformers, searches a LanceDB vector index. Returns code previews with file paths and line numbers. Covers C++, Python, and docs.
2. **`openmc_map.py`** — Structural repo map via aider/tree-sitter. Builds a cross-file reference graph, ranks files with PageRank relative to your focus files, then shows the top-ranked files as condensed code skeletons fitted to a token budget. Focus files are excluded (assumes you already have them). Caveat: the graph matches identifiers by name only — common names like `push_back` or `__init__` create false edges in the ranking.
3. **`openmc_lsp.py`** — LSP navigation via clangd. Talks to the C++ compiler frontend for symbol resolution. `definition`, `references`, `symbols`, and `related` commands with compiler accuracy — zero false edges. Requires clangd and compile_commands.json.

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

- **`openmc_search.py`**: "What code is conceptually related to X?" — broad discovery by meaning, cross-cutting concerns, Python and docs
- **`openmc_lsp.py`**: "Where is this C++ symbol defined, who calls it, and what files are truly connected to this one?" — compiler-accurate file:line locations, zero false positives
- **`openmc_map.py`**: "Show me the code structure of files neighboring my focus files" — PageRank-ranked code skeletons fitted to a token budget. Neighbor ranking is noisy for common identifiers; use `openmc_lsp.py related` for accurate C++ file connections

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help` (structural map)
- `.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help` (LSP navigation, C++ only)"
