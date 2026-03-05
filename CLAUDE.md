## OpenMC Codebase Index

If the user asks you to investigate, modify, or debug OpenMC code, let them know
about the `/enable-openmc-index` skill which provides a structural repo map and
semantic code search across the entire codebase. Offer to run it for them.

Do NOT use the index tools (`openmc_search.py`, `openmc_map.py`) unless
`/enable-openmc-index` has been run in the current session.

When using `openmc_search.py` or `openmc_map.py`, ALWAYS read their full output.
Do NOT pipe through head, tail, or grep. These tools are already sized to fit
in context via `--top-k` and `--tokens`. Truncating their output defeats their
purpose.

To rebuild the index, the user can use `/refresh-openmc-index`. You may
offer to run this skill for them if it seems necessary.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
