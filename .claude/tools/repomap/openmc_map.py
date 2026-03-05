#!/usr/bin/env python3
"""Generate a focused repo map around specific OpenMC files.

Uses aider's RepoMap to produce a condensed structural overview of the
codebase, ranked by relevance to the files you're currently working on.

Usage:
    openmc_map.py src/particle.cpp                    # Map around one file
    openmc_map.py src/simulation.cpp src/source.cpp   # Map around multiple files
    openmc_map.py --tokens 4096                       # Larger map (default: 2048)
    openmc_map.py                                     # Map of the whole repo (top-ranked files)

Examples:
    openmc_map.py src/particle_restart.cpp src/random_lcg.cpp
    openmc_map.py openmc/deplete/coupled_operator.py --tokens 4096
    openmc_map.py include/openmc/cell.h include/openmc/surface.h
"""

import argparse
import glob
import os
import sys
from pathlib import Path

OPENMC_ROOT = Path(__file__).resolve().parents[3]

# File patterns to include in the map
FILE_PATTERNS = [
    "src/**/*.cpp",
    "include/openmc/**/*.h",
    "openmc/**/*.py",
]


class TokenCounter:
    """Simple token counter that doesn't need an API model."""

    def token_count(self, text):
        # Rough approximation: ~4 chars per token for code
        return len(text) // 4


class FakeModel:
    """Minimal model stand-in for aider's RepoMap token counting."""

    def __init__(self):
        self._token_counter = TokenCounter()

    def token_count(self, text):
        return self._token_counter.token_count(text)


def get_all_files():
    """Collect all source files matching our patterns."""
    files = []
    for pattern in FILE_PATTERNS:
        for fp in sorted(OPENMC_ROOT.glob(pattern)):
            if "__pycache__" in str(fp):
                continue
            files.append(str(fp.relative_to(OPENMC_ROOT)))
    return files


def generate_map(focus_files=None, map_tokens=2048):
    """Generate a repo map, optionally focused on specific files.

    Args:
        focus_files: List of file paths to focus on. If None, generates
                     a general overview of the most important files.
        map_tokens: Approximate token budget for the map.

    Returns:
        The repo map as a string.
    """
    from aider.io import InputOutput
    from aider.repomap import RepoMap

    os.chdir(OPENMC_ROOT)

    io = InputOutput(yes=True)
    model = FakeModel()

    rm = RepoMap(
        map_tokens=map_tokens,
        root=str(OPENMC_ROOT),
        io=io,
        main_model=model,
    )

    all_files = get_all_files()

    # Normalize focus files to relative paths
    chat_fnames = []
    if focus_files:
        for f in focus_files:
            # Handle both absolute and relative paths
            fp = Path(f)
            if fp.is_absolute():
                try:
                    fp = fp.relative_to(OPENMC_ROOT)
                except ValueError:
                    pass
            rel = str(fp)
            if rel in all_files:
                chat_fnames.append(rel)
            else:
                # Try to find a match
                matches = [af for af in all_files if rel in af]
                if matches:
                    chat_fnames.extend(matches)
                else:
                    print(f"Warning: '{f}' not found in indexed files",
                          file=sys.stderr)

    # other_fnames = files NOT in chat_fnames
    other_fnames = [f for f in all_files if f not in chat_fnames]

    repo_map = rm.get_repo_map(chat_fnames, other_fnames)

    if not repo_map:
        return "No map generated. Try with different files or a larger --tokens budget."

    return repo_map


def main():
    parser = argparse.ArgumentParser(
        description="Generate a focused structural map of OpenMC code",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  %(prog)s src/particle_restart.cpp src/random_lcg.cpp
  %(prog)s openmc/deplete/coupled_operator.py --tokens 4096
  %(prog)s include/openmc/cell.h include/openmc/surface.h
  %(prog)s                                    # overview of whole repo""",
    )
    parser.add_argument(
        "files", nargs="*",
        help="Files to focus the map on (shows their structure and neighbors)")
    parser.add_argument(
        "--tokens", type=int, default=2048,
        help="Approximate token budget for the map (default: 2048)")

    args = parser.parse_args()

    # Suppress aider's scanning output
    repo_map = generate_map(
        focus_files=args.files if args.files else None,
        map_tokens=args.tokens,
    )
    print(repo_map)


if __name__ == "__main__":
    main()
