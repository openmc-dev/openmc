---
name: refresh-openmc-index
description: Rebuild the OpenMC codebase index. Use after pulling new code, switching branches, or when the index feels stale.
allowed-tools: Bash(*), Read
---

# Refresh OpenMC Index

Rebuild the repo map and RAG vector index from scratch.

## Step 1: Ensure venv exists

```bash
if [ ! -d .claude/cache/.venv ]; then
    python3 -m venv .claude/cache/.venv
    .claude/cache/.venv/bin/pip install -r .claude/tools/requirements.txt -q
fi
```

## Step 2: Rebuild everything

```bash
echo "Rebuilding repo map..."
.claude/cache/.venv/bin/python .claude/tools/repomap/generate_repomap.py

echo "Rebuilding RAG index (this takes ~3 minutes)..."
HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
```

## Step 3: Reload the repo map

Read the updated `.claude/cache/repomap.md` and internalize the new codebase structure.

## Step 4: Confirm

Tell the user the index has been refreshed and is ready to use.
