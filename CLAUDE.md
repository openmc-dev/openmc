## OpenMC Codebase Tools

For any task involving OpenMC code — investigating, modifying, debugging, or
reviewing PRs — let the user know about the `/enable-openmc-index` skill which
provides three code navigation tools. Offer to run it for them. Even for PR
reviews, these tools are important: a diff only shows what changed, not what
else in the codebase depends on or is affected by those changes. Note: the
first run builds a RAG vector index, which takes ~5 minutes on 10 CPU cores.
Subsequent sessions reuse the cached index.

Do NOT use the tools (`openmc_search.py`, `openmc_map.py`, `openmc_lsp.py`)
unless `/enable-openmc-index` has been run in the current session.

When using these tools, ALWAYS read their full output. Do NOT pipe through head,
tail, or grep. The tools already limit their output size via `--top-k` and
`--tokens` flags. Truncating their output defeats their purpose.

To rebuild the RAG search index after pulling new code, the user can use
`/refresh-openmc-index`.

### Important: use RAG search before grep

When exploring unfamiliar code or checking what a change might affect, use
`openmc_search.py` **before** reaching for `grep` or `Glob`. `grep` only finds
exact text matches — it cannot find code that does the same thing with different
naming. The RAG search finds code by semantic meaning, surfacing related code
across subsystems that you would otherwise miss entirely. Use RAG for discovery,
then `grep`/`Read` to drill into the specific files it surfaces.

When you already know the exact symbol name and need to trace its usage (e.g.,
"every line that writes to `progeny_per_particle`"), `grep` is the right tool
— don't force a RAG search for precise symbol lookups.

### Tool details

**`openmc_search.py`** — RAG semantic search. The codebase (C++, Python, and
RST docs) is chunked into overlapping fixed-size windows (~1000 chars, 25%
overlap) so every line of code is searchable. Chunks are embedded with
sentence-transformers and stored in a local LanceDB vector index. Your query is
embedded the same way, and the closest chunks are returned with file paths, line
numbers, and a code preview. Good for finding conceptually related code even
when naming differs (e.g., "particle RNG seeding" finds code across transport,
restart, and random ray modes). Returns `--top-k` results (default 10).

**`openmc_lsp.py`** — LSP navigation via clangd. Launches clangd as a subprocess
and queries it via the Language Server Protocol. Because clangd uses the actual
C++ compiler frontend (Clang), it resolves every symbol through the real type
system — namespaces, templates, overloads, and all. Commands:
- `symbols FILE` — list all symbols defined in a file with their types and lines
- `definition FILE:LINE` — jump to where the symbol at that line is defined
- `references FILE:LINE` — find every file and line that references that symbol
- `related FILE` — for each symbol defined in the file, find all external
  references, then rank other files by how many typed connections they share.
  Returns `--top-k` files (default 15) with the connecting symbol names.

Zero false edges — if it says two files are connected, they genuinely share
typed references. Requires clangd and `build/compile_commands.json` (automatically
generated when OpenMC is built with cmake).

**`openmc_map.py`** — Structural repo map via aider/tree-sitter. Tree-sitter
parses all C++ and Python source files to extract identifier definitions and
references. A cross-file reference graph is built (file A references a symbol
defined in file B → edge from A to B), then PageRank ranks files by importance
relative to your focus files. The top-ranked files are shown as condensed code
skeletons with class/function signatures and `⋮` elision markers, fitted to a
`--tokens` budget (default 2048). Focus files themselves are excluded from the
output (the assumption is you already have them in context). **Caveat**: the
reference graph matches identifiers by name only — tree-sitter has no type
information, so `std::vector::push_back` and `NeighborList::push_back` create
the same edges. This means files defining common method names (`push_back`,
`get`, `__init__`, `from_xml`, etc.) get inflated PageRank and appear as
neighbors when they may not be truly related. The name-matching can also be
useful — it surfaces files with identically-named functions that may need
parallel changes even though they have no typed connection. For precise C++
file connections, `openmc_lsp.py related` is more reliable; the repo map is
better for a broad structural overview or for Python code.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
