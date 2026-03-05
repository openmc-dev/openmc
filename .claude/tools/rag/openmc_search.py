#!/usr/bin/env python3
"""Semantic search across the OpenMC codebase and documentation.

Usage:
    openmc_search.py "query"                    # Search code (default)
    openmc_search.py "query" --docs             # Search documentation
    openmc_search.py "query" --all              # Search both code and docs
    openmc_search.py --related src/particle.cpp # Find related code
    openmc_search.py "query" --top-k 20         # Return more results

Examples:
    openmc_search.py "particle random number seed initialization"
    openmc_search.py "how to define tallies" --docs
    openmc_search.py --related src/simulation.cpp
"""

import argparse
import sys
from pathlib import Path

# Add tools dir to path
TOOLS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS_DIR / "rag"))

OPENMC_ROOT = Path(__file__).resolve().parents[3]
CACHE_DIR = OPENMC_ROOT / ".claude" / "cache"
INDEX_DIR = CACHE_DIR / "rag_index"


def get_db_and_embedder():
    """Load the LanceDB database and embedding provider."""
    import lancedb
    from embeddings import EmbeddingProvider

    if not INDEX_DIR.exists():
        print("ERROR: No index found. Run /enable-openmc-index first.",
              file=sys.stderr)
        sys.exit(1)

    db = lancedb.connect(str(INDEX_DIR))

    embedder = EmbeddingProvider()
    return db, embedder


def search_table(db, embedder, table_name, query, top_k):
    """Search a LanceDB table with a text query."""
    if table_name not in db.list_tables().tables:
        print(f"Table '{table_name}' not found in index.", file=sys.stderr)
        return []

    table = db.open_table(table_name)
    query_vec = embedder.embed_query(query)
    results = table.search(query_vec).limit(top_k).to_list()
    return results


def format_results(results, label=""):
    """Format search results for display."""
    if not results:
        return "No results found.\n"

    output = []
    if label:
        output.append(f"=== {label} ===\n")

    for i, r in enumerate(results, 1):
        filepath = r["filepath"]
        start = r["start_line"]
        end = r["end_line"]
        kind = r["kind"]
        symbol = r["symbol"]
        dist = r.get("_distance", 0)

        header = f"[{i}] {filepath}:{start}-{end} ({kind}"
        if symbol:
            header += f": {symbol}"
        header += f", dist={dist:.3f})"
        output.append(header)

        # Show text preview (first 500 chars)
        text = r["text"][:500]
        if len(r["text"]) > 500:
            text += "\n  ..."
        # Indent the text
        for line in text.split("\n"):
            output.append(f"  {line}")
        output.append("")

    return "\n".join(output)


def search_related(db, embedder, filepath, top_k):
    """Find code related to a given file."""
    if "code" not in db.list_tables().tables:
        print("No 'code' table in index.", file=sys.stderr)
        return []

    table = db.open_table("code")

    # Normalize filepath
    fp = filepath
    if filepath.startswith("/"):
        try:
            fp = str(Path(filepath).relative_to(OPENMC_ROOT))
        except ValueError:
            pass

    # Get chunks from target file
    try:
        target_chunks = table.search().where(
            f"filepath = '{fp}'"
        ).limit(50).to_list()
    except Exception:
        # LanceDB where clause might not work in all versions
        # Fall back to fetching all and filtering
        all_data = table.to_pandas()
        target_rows = all_data[all_data["filepath"] == fp]
        if target_rows.empty:
            print(f"No chunks found for '{fp}'", file=sys.stderr)
            return []
        target_chunks = target_rows.head(50).to_dict("records")

    if not target_chunks:
        print(f"No chunks found for '{fp}'", file=sys.stderr)
        return []

    # Combine top chunks as the query
    combined_text = " ".join(c["text"][:200] for c in target_chunks[:5])
    query_vec = embedder.embed_query(combined_text)

    # Search excluding the source file
    results = table.search(query_vec).limit(top_k + 10).to_list()
    # Filter out same file
    results = [r for r in results if r["filepath"] != fp][:top_k]
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Semantic search across OpenMC codebase and docs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  %(prog)s "particle random number seed initialization"
  %(prog)s "how to define tallies" --docs
  %(prog)s "weight window variance reduction" --all
  %(prog)s "where is cross section data loaded" --top-k 15
  %(prog)s --related src/simulation.cpp
  %(prog)s --related src/particle_restart.cpp --top-k 5""",
    )
    parser.add_argument("query", nargs="?", help="Search query")
    parser.add_argument("--docs", action="store_true",
                        help="Search documentation instead of code")
    parser.add_argument("--all", action="store_true",
                        help="Search both code and documentation")
    parser.add_argument("--related", metavar="FILE",
                        help="Find code related to a given file")
    parser.add_argument("--top-k", type=int, default=10,
                        help="Number of results (default: 10)")
    args = parser.parse_args()

    if not args.query and not args.related:
        parser.print_help()
        sys.exit(1)

    db, embedder = get_db_and_embedder()

    if args.related:
        results = search_related(db, embedder, args.related, args.top_k)
        print(format_results(results, f"Code related to {args.related}"))
    elif args.all:
        code_results = search_table(db, embedder, "code", args.query, args.top_k)
        doc_results = search_table(db, embedder, "docs", args.query, args.top_k)
        print(format_results(code_results, "Code"))
        print(format_results(doc_results, "Documentation"))
    elif args.docs:
        results = search_table(db, embedder, "docs", args.query, args.top_k)
        print(format_results(results, "Documentation"))
    else:
        results = search_table(db, embedder, "code", args.query, args.top_k)
        print(format_results(results, "Code"))


if __name__ == "__main__":
    main()
