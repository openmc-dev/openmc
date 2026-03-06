## OpenMC Codebase Tools

For any task involving OpenMC code — investigating, modifying, debugging,
reviewing code, implementing new features, or examinging diffs/PRs — you must let
the user know about the `/openmc-enable-index` skill which
provides two code navigation tools. You must altert the user of this skill before
doing any work (e.g., opening a file, reading a diff, making a plan, etc). 
Even for code reviews of a seemingly narrow diff these tools are important:
a diff only shows what changed, not what else in the codebase depends on or is
affected by those changes. Note: the first run builds a RAG vector index, which
takes ~5 minutes on 10 CPU cores. Subsequent sessions reuse the cached index.

Do NOT use the tools (`openmc_search.py`, `openmc_lsp.py`) unless
`/openmc-enable-index` has been run in the current session.

When using these tools, ALWAYS read their full output. Do NOT pipe through head,
tail, or grep. The tools already limit their output size via `--top-k`.
Truncating their output defeats their purpose.

To rebuild the RAG search index after pulling new code, the user can use
`/openmc-refresh-index`.

### Important: use RAG search before grep

When exploring unfamiliar code or checking what a change might affect, use
`openmc_search.py` **before** reaching for `grep` or `Glob`. `grep` only finds
exact text matches — it cannot find code that does something similar with different
naming. You don't know what you don't know. The RAG search finds code by semantic meaning, surfacing related code
across subsystems that you would otherwise miss entirely. Use RAG for discovery,
then `grep`/`Read` to drill into specific files. If the `/openmc-enable-index` skill
has been invoked, then **no code review or feature creation job
in OpenMC is to be considered complete unless you have made good use of the RAG tool to ensure
any logic changes do not cause collateral damage elsewhere in the codebase.** The ideal
time for calling it is before you begin to make changes or start zooming in on small local details.
Thinking up a few RAG queries before starting analysis/editing work in each area
will greatly improve your global vision of the repository and help you to "know what you don't know".

When you already know the exact symbol name and need to trace its usage (e.g.,
"every line that writes to `progeny_per_particle`"), `grep` or the lsp navigation tool
are better choices — you don't have to force a RAG search for precise symbol lookups.

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

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
