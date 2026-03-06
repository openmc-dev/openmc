---
name: refresh-openmc-index
description: Rebuild the OpenMC RAG search index. Use after pulling new code or switching branches. The LSP tool does not need refreshing.
allowed-tools: Bash(*)
---

# Refresh OpenMC Index

Rebuild the RAG semantic search index from scratch. Only this index needs
refreshing — the LSP tool always works on the current code.

## Step 1: Ensure venv exists

```bash
if [ ! -d .claude/cache/.venv ]; then
    python3 -m venv .claude/cache/.venv
    .claude/cache/.venv/bin/pip install -r .claude/tools/requirements.txt -q
fi
```

## Step 2: Rebuild the RAG index

```bash
echo "Rebuilding RAG index (takes ~5 minutes on 10 CPU cores)..."
HF_HUB_DISABLE_TELEMETRY=1 .claude/cache/.venv/bin/python .claude/tools/rag/indexer.py
```

## Step 3: Confirm

Tell the user the index has been refreshed and is ready to use.
