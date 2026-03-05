## OpenMC Codebase Tools

If the user asks you to investigate, modify, or debug OpenMC code, let them know
about the `/enable-openmc-index` skill which provides three code navigation tools.
Offer to run it for them.

Do NOT use the tools (`openmc_search.py`, `openmc_map.py`, `openmc_lsp.py`)
unless `/enable-openmc-index` has been run in the current session.

When using these tools, ALWAYS read their full output. Do NOT pipe through head,
tail, or grep. The tools already limit their output size via `--top-k` and
`--tokens` flags. Truncating their output defeats their purpose.

### Tool summary

**`openmc_search.py`** — RAG semantic search. Embeds your query and searches a
vector index of the codebase. Good for finding conceptually related code even
when naming differs (e.g., "particle RNG seeding" finds code across transport,
restart, and random ray modes). Covers C++, Python, and docs.

**`openmc_lsp.py`** — LSP navigation via clangd. Uses the C++ compiler's own
type system to resolve symbols. `definition`, `references`, and `related`
commands give compiler-accurate results with zero false edges. Use this when
you need to know exactly which files reference a C++ symbol or are connected
to a C++ file. Requires clangd and `build/compile_commands.json`.

**`openmc_map.py`** — Structural repo map via aider/tree-sitter. Given focus files
you already have open, shows condensed code skeletons (class/function signatures)
of their *neighbors* — other files that share symbols with them. The focus files
themselves are excluded (you already have them in context). **Caveat**: the
underlying graph matches identifiers by name only (tree-sitter has no type
information), so files defining common method names like `push_back`, `get`,
`__init__`, or `from_xml` may appear as neighbors when they are not truly
related. For finding which files are genuinely connected to a C++ file, prefer
`openmc_lsp.py related`.

To rebuild the RAG search index after pulling new code, the user can use
`/refresh-openmc-index`.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
