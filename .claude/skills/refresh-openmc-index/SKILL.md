---
name: refresh-openmc-index
description: Rebuild the OpenMC RAG search index. Use after pulling new code or switching branches. The repo map tool does not need refreshing (it generates maps on the fly).
allowed-tools: Bash(*)
---

# Refresh OpenMC Index

Rebuild the RAG vector index from scratch. Note: only the semantic search index needs refreshing. The repo map tool (`openmc_map.py`) always works on the current code.

## Step 1: Ensure venv exists

```bash
if [ ! -d .claude/cache/.venv ]; then
    python3 -m venv .claude/cache/.venv
    .claude/cache/.venv/bin/pip install -r .claude/tools/requirements.txt -q
fi
```

## Step 2: Rebuild the RAG index

```bash
echo "Rebuilding RAG index (this takes ~3 minutes)..."
HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
```

## Step 3: Confirm

Tell the user the index has been refreshed and is ready to use.
