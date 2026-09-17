.. _devguide_agentic_tools:

===========================
Agentic Development Tools
===========================

OpenMC ships a set of tools designed for AI coding agents (such as
`Claude Code`_) that agents can use to navigate and understand the codebase.

.. _Claude Code: https://claude.ai/code

Motivation
----------

Agentic tools like Claude Code are skilled at using grep to navigate and
understand large code bases. However, grep can only find exact text matches —
it cannot discover code that is *conceptually* related but uses different
naming. Without a "global view" of the codebase that a human developer will
build up over time, the agent is generally blind to any file it hasn't
tokenized fully. While it can grep to see who else calls a function, it
remains blind if other areas might be related but not share identical naming
conventions.

This problem is mitigated somewhat by using a model with a longer context
window. OpenMC has somewhere around ~1 million tokens of C++ and ~1 million
tokens of python. While Claude Code in early 2026 only has a context window
of 200k tokens, beta versions have extended context windows of 1M tokens,
and it's not unreasonable to assume that models may be available in the near
future that greatly exceed these limits.

However, even assuming the entire repository can be fit within a context
window, there are several downsides to doing this.
`Model performance degrades significantly as context size increases`_.
Benchmark results are
greatly improved if the model has less garbage to pick through. Additionally, API usage
is typically billed as tokens in/out per turn. As the context file
grows these costs become much larger. As such, there is still significant
motivation to solving the above problem, so as to ensure only relevant
information is drawn into context so as to maximize model performance and
minimize costs.

Setup
-----

The tools are registered as an `MCP (Model Context Protocol)`_ server in
``.mcp.json`` at the repository root. AI agents that support MCP (such as
Claude Code) discover them automatically on session start. The underlying
Python scripts can also be run directly from the command line.

All tools run entirely locally — no API keys or external service accounts are
required. Python dependencies are installed automatically into an isolated
virtual environment at ``.claude/cache/.venv/`` on first use.

.. _Model performance degrades significantly as context size increases: https://www.anthropic.com/news/claude-opus-4-6
.. _MCP (Model Context Protocol): https://modelcontextprotocol.io

RAG Semantic Search
-------------------

The RAG (Retrieval-Augmented Generation) semantic search addresses this
problem — it finds code by meaning, not just text match, surfacing related code
across subsystems that ``grep`` would miss entirely. Two MCP tools are provided:

- **openmc_rag_search** — Given a natural-language query, returns the most
  relevant code chunks with file paths, line numbers, and a preview. Can search
  code, documentation, or both. Can also find code related to a given file.
- **openmc_rag_rebuild** — Rebuilds the search index. Should be called after
  pulling new code or switching branches.

How it works
^^^^^^^^^^^^

The search pipeline runs entirely on your local CPU:

1. **Chunking.** All C++, Python, and RST files are split into overlapping
   fixed-size windows (~1000 characters, 25% overlap). This ensures every line
   of code appears in at least one chunk and most lines appear in two.

2. **Embedding.** Each chunk is embedded into a 384-dimensional vector using
   the `all-MiniLM-L6-v2`_ sentence-transformer model (22 million parameters).
   This model runs on CPU with no GPU required. No API key is needed — the
   model weights are downloaded once from Hugging Face and cached locally.

3. **Indexing.** The vectors are stored in a local LanceDB_ database on disk.
   Building the full index takes approximately 5 minutes on a machine with
   10 CPU cores. The index is stored in ``.claude/cache/rag_index/`` and
   persists across sessions.

4. **Searching.** Your query is embedded using the same model, and the closest
   chunks are retrieved by vector similarity. Results include the file path,
   line range, file type, similarity distance, and a text preview.

.. _all-MiniLM-L6-v2: https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2
.. _LanceDB: https://lancedb.com

Requirements
^^^^^^^^^^^^

No system dependencies beyond **Python 3.12+** with ``pip``. An internet
connection is required on first use to download the Python packages and
embedding model weights; subsequent runs are fully offline. The Python packages
(``sentence-transformers``, ``lancedb``) and their dependencies (including
PyTorch, ~2GB) are installed automatically into an isolated virtual environment
on first use.
