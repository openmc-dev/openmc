---
name: enable-openmc-index
description: Enable the OpenMC codebase index for this session. Provides semantic code search, structural repo mapping, and LSP-based code navigation. Run this when investigating, modifying, or debugging OpenMC code.
allowed-tools: Bash(*), Read
---

# Enable OpenMC Index

Set up (if needed) and activate the OpenMC codebase index for this session. This gives you three tools:
1. **Semantic search** (`openmc_search.py`) - Find related code across the codebase by concept
2. **Structural map** (`openmc_map.py`) - See condensed code structure of files and their neighbors
3. **LSP navigation** (`openmc_lsp.py`) - Compiler-accurate go-to-definition, find-references, and related-file discovery for C++ code (requires clangd and compile_commands.json)

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

## Step 3: Check LSP tool prerequisites

The LSP tool (`openmc_lsp.py`) requires clangd and `compile_commands.json`. Check if they're available:

```bash
if command -v clangd &>/dev/null || compgen -c clangd- 2>/dev/null | head -1 | grep -q .; then
    echo "CLANGD_AVAILABLE"
else
    echo "CLANGD_MISSING"
fi
if [ -f build/compile_commands.json ]; then
    echo "COMPILE_COMMANDS_AVAILABLE"
else
    echo "COMPILE_COMMANDS_MISSING"
fi
```

If clangd is missing, tell the user: "The LSP tool needs clangd. Install with `apt-get install clangd` (or `clangd-15`, `clangd-16`, etc.)."
If compile_commands.json is missing, tell the user: "The LSP tool needs compile_commands.json. Generate with `cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`."

## Step 4: Learn the tool APIs

Run `--help` for all tools to see their full APIs:

```bash
.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help
.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help
.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help
```

Read and internalize the output so you know all available options.

## Step 5: Confirm activation

Tell the user the OpenMC index is active and briefly describe the tools:
- **Semantic search**: Find related code by concept (e.g., "particle seed initialization")
- **Structural map**: See condensed code structure around specific files
- **LSP navigation** (if available): Compiler-accurate definition/references/related-files for C++

## How to use the tools after activation

**Typical workflow:**

1. Use `openmc_search.py` to discover which files are relevant to your task
2. Use `openmc_lsp.py related` to find files connected by real typed references (C++ only)
3. Use `openmc_map.py` on those files to see their code structure
4. Use Read/Grep to dive into the specific code you need to change

**When to use semantic search** (`openmc_search.py`):
- Investigating how a change might affect other parts of the codebase
- Finding code that does something conceptually similar but with different naming
- Discovering cross-cutting concerns across run modes
- Searching Python code and documentation

**When to use LSP navigation** (`openmc_lsp.py`):
- Finding exactly where a C++ symbol is defined (`definition`)
- Finding all callers of a C++ function (`references`)
- Discovering which C++ files are connected by real typed references (`related`)
- Best for C++ — uses the compiler's own type system, zero false edges

**When to use the repo map** (`openmc_map.py`):
- Understanding the structure of unfamiliar files before modifying them
- Seeing what classes/methods neighbor the code you're working on
- Getting a condensed overview of a subsystem (pass multiple files)
- Works for both C++ and Python files

## Subagent guidance

When spawning subagents that will investigate or modify OpenMC code, include in their prompt:

"The OpenMC code index is available. Run `--help` on these tools to see their full API:
- `.claude/cache/.venv/bin/python .claude/tools/rag/openmc_search.py --help` (semantic search)
- `.claude/cache/.venv/bin/python .claude/tools/repomap/openmc_map.py --help` (structural map)
- `.claude/cache/.venv/bin/python .claude/tools/lsp/openmc_lsp.py --help` (LSP navigation, C++ only)"
