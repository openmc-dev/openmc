---
name: enable-openmc-index
description: Enable the OpenMC codebase tools for this session. Provides semantic code search and LSP-based C++ code navigation.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase tools for this session:

1. **`openmc_search.py`** — RAG semantic search. Chunks code into overlapping fixed-size windows, embeds locally on CPU with the all-MiniLM-L6-v2 model (22M parameters, no GPU or API key needed), and searches a LanceDB vector index. Returns code previews with file paths and line numbers. Covers C++, Python, and docs.
2. **`openmc_lsp.py`** — LSP (Language Server Protocol) navigation via clangd. Talks to the C++ compiler frontend for symbol resolution. `definition`, `references`, `symbols`, and `related` commands with compiler accuracy — zero false edges. Requires clangd and that OpenMC has been built (for `build/compile_commands.json`).

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

The semantic search tool needs a pre-built vector index. The LSP tool works without it.

```bash
if [ ! -d .claude/cache/rag_index ]; then
    echo "Building RAG index for the first time (takes ~5 minutes on 10 CPU cores)..."
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
.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help
```

Read and internalize the output so you know all available options.

## Step 4: Demonstrate the RAG tool to yourself

Run this search and read the results carefully:

```bash
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py "particle weight adjustment variance reduction"
```

Notice how the results span multiple files across different subsystems (weight
windows, random ray, variance reduction) — files you would never find with
`grep "particle weight"` because they use different terminology for related
concepts. This is the key advantage: **the RAG search finds code by meaning,
not by exact text match.** When you `grep` for `init_particle`, you miss
`initialize_seeds` which does the same thing with a different name. The RAG
search finds both.

Your trained instinct is to reach for `grep` and `Read` — those are excellent
tools for precise symbol lookup, but they are blind to semantic relationships.
Use each tool for what it's best at:
- **RAG search**: discovery, exploration, "what else in the codebase relates to
  this concept?", checking what a change might affect across subsystems
- **grep/Read**: precise symbol tracing, "every line that writes to `variable_x`"

Don't force RAG searches for exact symbol lookups, and don't rely on grep alone
for broad exploration.

## Step 5: Demonstrate the LSP tool to yourself

Run this references query and read the results carefully:

```bash
.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py references src/tallies/tally.cpp:835
```

Line 835 is `Tally::reset()`. The LSP tool uses the C++ compiler frontend to
resolve this — it returns only references to **this specific** `reset()`, not
the other 3 classes that also define `void reset()` (Timer, ParticleData,
SharedArray). Compare with `grep 'reset()'` which returns 62 mixed hits across
20 files including vendor code. The LSP tool gives you the exact 10 files that
call or reference `Tally::reset()`, with line numbers — zero false positives.

This is why the LSP tool exists: **`grep` matches text, LSP resolves types.**
When a common method name like `reset`, `get`, `size`, or `create` is used by
multiple classes, `grep` gives you a haystack. LSP gives you the needle.

Use each tool for what it's best at:
- **LSP**: "who calls *this specific* C++ method?" — type-accurate references
- **grep**: "every line containing this unique string" — fast exact text match

## Step 6: Confirm activation

Tell the user the tools are active. If the LSP tool can't find clangd or
`build/compile_commands.json`, it will report this itself when invoked.

## When to use each tool

- **`openmc_search.py`**: "What code is conceptually related to X?" — broad discovery by meaning, cross-cutting concerns, Python and docs. **Use this before grep when exploring unfamiliar code or checking what a change might affect.**
- **`openmc_lsp.py`**: "Where is this C++ symbol defined, who calls it, and what files are truly connected to this one?" — compiler-accurate file:line locations, zero false positives

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help` (LSP navigation, C++ only)"
