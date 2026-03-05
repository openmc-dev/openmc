## OpenMC Codebase Tools

If the user asks you to investigate, modify, or debug OpenMC code, let them know
about the `/enable-openmc-index` skill which provides semantic code search,
structural repo mapping, and LSP-based C++ code navigation. Offer to run it.

Do NOT use the tools (`openmc_search.py`, `openmc_map.py`, `openmc_lsp.py`)
unless `/enable-openmc-index` has been run in the current session.

When using these tools, ALWAYS read their full output. Do NOT pipe through head,
tail, or grep. The tools already limit their output size via `--top-k` and
`--tokens` flags. Truncating their output defeats their purpose.

To rebuild the RAG search index after pulling new code, the user can use
`/refresh-openmc-index`.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
