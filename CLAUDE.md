## OpenMC Codebase Tools

Two MCP tools are available for navigating the OpenMC codebase. They are
registered in `.mcp.json` and appear automatically in every session.

### Tool overview

**`openmc_rag_search`** — Semantic search across the codebase (C++, Python, RST
docs). Finds code by meaning, not just text match. Surfaces related code across
subsystems even when naming differs (e.g., "particle RNG seeding" finds code
across transport, restart, and random ray modes — files you would never find
with `grep "particle seed"`). The index uses a small 22M-param embedding model
(384-dim). Phrase-level natural-language queries work much better than single
keywords or symbol names.

**`openmc_rag_rebuild`** — Rebuild the RAG vector index. Call after pulling new
code or switching branches. The first RAG search of each session will report
the index status and ask whether to rebuild — you can also call this explicitly.

**First-call behavior:** The first `openmc_rag_search` call of each session
returns an index status message instead of search results. When this happens,
you MUST use the `AskUserQuestion` tool to present the rebuild/use-existing
choice to the user. Do not ask conversationally — always use the widget.
Do not skip this step even if the index looks current — the user may have
uncommitted changes that warrant a rebuild.

### Important: use RAG search before grep

When exploring unfamiliar code or checking what a change might affect, use
`openmc_rag_search` **before** reaching for `grep` or `Glob`. `grep` only finds
exact text matches — it cannot find code that does something similar with
different naming. You don't know what you don't know. The RAG search finds code
by semantic meaning, surfacing related code across subsystems that you would
otherwise miss entirely. Use RAG for discovery, then `grep`/`Read` to drill into
specific files.

**No code review or feature creation job in OpenMC is to be considered complete
unless you have made good use of `openmc_rag_search` to ensure any logic changes
do not cause collateral damage elsewhere in the codebase.** The ideal time for
calling it is before you begin to make changes or start zooming in on small local
details. Thinking up a few RAG queries before starting analysis/editing work in
each area will greatly improve your global vision of the repository and help you
to "know what you don't know".

When you already know the exact symbol name and need to trace its usage (e.g.,
"every line that writes to `progeny_per_particle`"), `grep` is the better
choice — you don't have to force a RAG search for precise symbol lookups.

### When to use each tool

- **`openmc_rag_search`**: "What code is conceptually related to X?" — broad
  discovery by meaning, cross-cutting concerns, Python and docs. **Use this
  before grep when exploring unfamiliar code or checking what a change might
  affect.**
- **`grep`/`Glob`/`Read`**: Precise text match, unique string lookup, reading
  specific files. Best when you know the exact symbol name.

### Why global awareness matters

An agent reviewed a large OpenMC PR using only diff, grep, and Read. It found
1 of 11 serious bugs. Here is its own post-mortem:

> **I treated the diff as a closed system.** I verified internal consistency of
> the changed code obsessively, but never built a global understanding of how
> the changed code fits into the wider codebase. The diff altered assumptions
> that code elsewhere silently relied on — but I couldn't see that because I
> never looked beyond the diff. I couldn't see the forest for the trees.
>
> **Why I resisted RAG:** Overconfidence. My internal model was "I can see the
> diff, I understand the data structures, I can trace the logic." The diff felt
> self-contained. RAG felt like it would return noisy results about tangentially
> related code. But in a codebase this large, changes in one subsystem can
> quietly break invariants that distant subsystems depend on — and you need
> global awareness to foresee that.
>
> **In the post-mortem**, I re-ran the RAG queries I should have run during the
> review. They directly surfaced the files containing the bugs I missed — files
> I never thought to open because they weren't in the diff.

The takeaway: **use RAG throughout your work to maintain global awareness.**
Before diving into details, ask "what else in this codebase depends on the
behavior being changed?" As you explore each area, keep querying to build your
mental map of affected subsystems. The diff tells you *what* changed; RAG tells
you *what else cares*.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
