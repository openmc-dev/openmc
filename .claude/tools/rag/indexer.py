#!/usr/bin/env python3
"""Build the RAG vector index for the OpenMC codebase.

This is the index-building half of the RAG pipeline. All operations are local
once the embedding model has been downloaded and cached (see embeddings.py for
details on model download, caching, and network behavior). It walks the repo,
chunks every
C++/Python/RST file (via chunker.py), embeds all chunks into 384-dim vectors
(via embeddings.py), and stores them in a local LanceDB database on disk. The
result is a .claude/cache/rag_index/ directory containing two tables — "code"
and "docs" — that openmc_search.py queries at search time.

Building the full index takes ~5 minutes on a 10-core machine. The bottleneck
is the embedding step (running all chunks through the MiniLM model on CPU).

Can be run standalone:  python indexer.py
Or called programmatically:  from indexer import build_index; build_index()
The MCP server (openmc_mcp_server.py) uses the latter when the agent calls
openmc_rag_rebuild.
"""

import lancedb
import sys
import time
from pathlib import Path

# This file lives at .claude/tools/rag/indexer.py. The sys.path insert lets
# us import sibling modules (embeddings, chunker) when run as a standalone
# script. When imported from the MCP server, the server has already done this.
TOOLS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS_DIR / "rag"))

from embeddings import EmbeddingProvider
from chunker import chunk_file


OPENMC_ROOT = Path(__file__).resolve().parents[3]
CACHE_DIR = OPENMC_ROOT / ".claude" / "cache"
INDEX_DIR = CACHE_DIR / "rag_index"

CODE_PATTERNS = [
    "src/**/*.cpp",
    "include/openmc/**/*.h",
    "openmc/**/*.py",
    "tests/**/*.py",
    "examples/**/*.py",
]

DOC_PATTERNS = [
    "docs/**/*.rst",
]


def collect_chunks(patterns, openmc_root):
    """Collect all chunks from files matching the given patterns."""
    chunks = []
    for pattern in patterns:
        for filepath in sorted(openmc_root.glob(pattern)):
            if "__pycache__" in str(filepath):
                continue
            file_chunks = chunk_file(filepath, openmc_root)
            chunks.extend(file_chunks)
    return chunks


def build_index():
    """Build or rebuild the complete vector index."""
    start = time.time()

    # Collect all chunks
    print("Collecting code chunks...")
    code_chunks = collect_chunks(CODE_PATTERNS, OPENMC_ROOT)
    print(f"  {len(code_chunks)} code chunks")

    print("Collecting doc chunks...")
    doc_chunks = collect_chunks(DOC_PATTERNS, OPENMC_ROOT)
    print(f"  {len(doc_chunks)} doc chunks")

    all_chunks = code_chunks + doc_chunks
    if not all_chunks:
        print("ERROR: No chunks collected!", file=sys.stderr)
        sys.exit(1)

    # Create embeddings
    all_texts = [c["text"] for c in all_chunks]
    print("Creating embedding provider...")
    embedder = EmbeddingProvider()
    print(f"  dim={embedder.dim}")

    print("Embedding chunks...")
    all_embeddings = embedder.embed(all_texts)

    # Build LanceDB tables
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    db = lancedb.connect(str(INDEX_DIR))

    # Separate code vs doc records by index (code_chunks come first in all_chunks)
    n_code = len(code_chunks)
    code_records = []
    doc_records = []
    for i, (chunk, emb) in enumerate(zip(all_chunks, all_embeddings)):
        record = {
            "text": chunk["text"],
            "filepath": chunk["filepath"],
            "kind": chunk["kind"],
            "symbol": chunk.get("symbol", ""),
            "start_line": chunk.get("start_line", 0),
            "end_line": chunk.get("end_line", 0),
            "vector": emb,
        }
        if i < n_code:
            code_records.append(record)
        else:
            doc_records.append(record)

    # Create tables (drop existing)
    result = db.table_names() if hasattr(db, "table_names") else db.list_tables()
    existing = result.tables if hasattr(result, "tables") else list(result)
    for table_name in ("code", "docs"):
        if table_name in existing:
            db.drop_table(table_name)

    if code_records:
        db.create_table("code", code_records)
        print(f"  Created 'code' table: {len(code_records)} rows")

    if doc_records:
        db.create_table("docs", doc_records)
        print(f"  Created 'docs' table: {len(doc_records)} rows")

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s")


if __name__ == "__main__":
    build_index()
