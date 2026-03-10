## OpenMC Codebase Tools

Read the FULL `AGENTS.md` in this directory before starting work. It contains
project context, coding conventions, and documentation of the RAG search tools
registered in `.mcp.json`.

### Claude Code-specific: first-call behavior

The first `openmc_rag_search` call of each session returns an index status
message instead of search results. When this happens, you MUST use the
`AskUserQuestion` tool to present the rebuild/use-existing choice to the user.
Do not ask conversationally — always use the widget. Do not skip this step even
if the index looks current — the user may have uncommitted changes that warrant
a rebuild.
