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

### Why RAG matters

OpenMC is large enough that changes in one subsystem can silently break
invariants that distant subsystems depend on — and those distant files often
use different naming, so grep won't find them. The RAG search finds code by
meaning, surfacing files you wouldn't have thought to open.

An agent reviewed a large OpenMC PR without RAG. It found 1 of 11 serious
bugs. Its post-mortem:

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

The takeaway: when reviewing or modifying code, ask yourself "what else in this
codebase might depend on the behavior I'm changing?" If you aren't sure, that's
a good time for a RAG query. It won't replace the grep-based investigation you
should already be doing — but it can surface files you wouldn't have thought to
open.

## Additional OpenMC info

Read the FULL `AGENTS.md` in this directory also before starting work.
