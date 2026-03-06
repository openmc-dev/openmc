"""Chunk OpenMC source files and documentation for RAG indexing.

Uses fixed-size overlapping windows so every line of code is searchable.
Window size is tuned to fit within the MiniLM embedding model's 256-token
context (~1000 chars). 25% overlap ensures most content appears in at least
two chunks.
"""

from pathlib import Path

# ~256 tokens for MiniLM. 1 token ≈ 4 chars for code.
WINDOW_CHARS = 1000
# 25% overlap — most lines appear in at least 2 chunks
STRIDE_CHARS = 750
MIN_CHUNK_CHARS = 50

SUPPORTED_EXTENSIONS = {".cpp", ".h", ".py", ".rst"}


def chunk_file(filepath, openmc_root):
    """Chunk a single file into overlapping fixed-size windows."""
    filepath = Path(filepath)
    if filepath.suffix not in SUPPORTED_EXTENSIONS:
        return []

    rel = str(filepath.relative_to(openmc_root))
    try:
        content = filepath.read_text(errors="replace")
    except Exception:
        return []

    if len(content) < MIN_CHUNK_CHARS:
        return []

    kind = _file_kind(filepath)

    # Build a char-offset → line-number map
    line_starts = []
    offset = 0
    for line in lines:
        line_starts.append(offset)
        offset += len(line) + 1  # +1 for newline

    chunks = []
    start = 0
    while start < len(content):
        end = min(start + WINDOW_CHARS, len(content))

        # Snap end to a line boundary to avoid splitting mid-line
        if end < len(content):
            newline_pos = content.rfind("\n", start, end)
            if newline_pos > start:
                end = newline_pos + 1

        text = content[start:end].strip()
        if len(text) >= MIN_CHUNK_CHARS:
            start_line = _offset_to_line(line_starts, start)
            end_line = _offset_to_line(line_starts, end - 1)
            chunks.append({
                "text": text,
                "filepath": rel,
                "kind": kind,
                "symbol": "",
                "start_line": start_line,
                "end_line": end_line,
            })

        start += STRIDE_CHARS

    return chunks


def _file_kind(filepath):
    """Map file extension to a kind label."""
    ext = filepath.suffix
    if ext in (".cpp", ".h"):
        return "cpp"
    elif ext == ".py":
        return "py"
    elif ext == ".rst":
        return "doc"
    return "other"


def _offset_to_line(line_starts, offset):
    """Convert a character offset to a 1-based line number."""
    # Binary search for the line containing this offset
    lo, hi = 0, len(line_starts) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if line_starts[mid] <= offset:
            lo = mid
        else:
            hi = mid - 1
    return lo + 1  # 1-based
