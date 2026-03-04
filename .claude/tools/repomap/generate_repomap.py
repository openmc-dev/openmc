#!/usr/bin/env python3
"""Generate a structural repo map of the OpenMC codebase.

Uses tree-sitter to parse C++ and Python files, extracts class/function
signatures, builds a cross-file reference graph, applies PageRank-like
ranking, and outputs a concise markdown map grouped by subsystem.

Output:
  .claude/cache/repomap.md        - Concise map (~170 lines) for agent context
  .claude/cache/repomap_full.json - Full symbol data for other tools
"""

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import tree_sitter_cpp as tscpp
import tree_sitter_python as tspy
from tree_sitter import Language, Parser

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

OPENMC_ROOT = Path(__file__).resolve().parents[3]  # .claude/tools/repomap -> repo root
CACHE_DIR = OPENMC_ROOT / ".claude" / "cache"

CPP_PATTERNS = ["src/**/*.cpp", "include/openmc/**/*.h"]
PY_PATTERNS = ["openmc/**/*.py"]

# Map file path keywords to logical subsystems
SUBSYSTEM_RULES = [
    # (path substring, subsystem name)
    ("random_ray", "Random Ray Solver"),
    ("tallies/", "Tallies & Filters"),
    ("deplete", "Depletion"),
    ("mgxs", "Multi-Group Cross Sections"),
    ("data/", "Nuclear Data"),
    ("lattice", "Geometry"),
    ("universe", "Geometry"),
    ("surface", "Geometry"),
    ("cell", "Geometry"),
    ("geometry", "Geometry"),
    ("dagmc", "Geometry (DAGMC)"),
    ("mesh", "Mesh"),
    ("material", "Materials"),
    ("nuclide", "Nuclear Data"),
    ("cross_section", "Nuclear Data"),
    ("thermal", "Nuclear Data"),
    ("wmp", "Nuclear Data"),
    ("particle", "Particle Transport"),
    ("physics", "Particle Transport"),
    ("collision", "Particle Transport"),
    ("photon", "Particle Transport"),
    ("bremsstrahlung", "Particle Transport"),
    ("secondary_", "Particle Transport"),
    ("reaction", "Particle Transport"),
    ("source", "Sources & Distributions"),
    ("distribution", "Sources & Distributions"),
    ("eigenvalue", "Eigenvalue Solver"),
    ("cmfd", "CMFD Acceleration"),
    ("weight_window", "Variance Reduction"),
    ("tally", "Tallies & Filters"),
    ("filter", "Tallies & Filters"),
    ("trigger", "Tallies & Filters"),
    ("plot", "Plotting"),
    ("track", "Plotting"),
    ("volume_calc", "Volume Calculation"),
    ("random_lcg", "Random Number Generation"),
    ("random_dist", "Random Number Generation"),
    ("settings", "Settings & Configuration"),
    ("simulation", "Simulation Control"),
    ("initialize", "Simulation Control"),
    ("finalize", "Simulation Control"),
    ("state_point", "I/O & Serialization"),
    ("summary", "I/O & Serialization"),
    ("hdf5", "I/O & Serialization"),
    ("xml_interface", "I/O & Serialization"),
    ("output", "I/O & Serialization"),
    ("bank", "Particle Banking"),
    ("event", "Event-Based Transport"),
    ("error", "Utilities"),
    ("string_utils", "Utilities"),
    ("math_functions", "Utilities"),
    ("file_utils", "Utilities"),
    ("timer", "Utilities"),
    ("memory", "Utilities"),
    ("position", "Utilities"),
    ("model", "Model Builder (Python)"),
    ("stats", "Statistics (Python)"),
    ("lib/", "C API Bindings (Python)"),
]

MAX_REPOMAP_LINES = 160

# ---------------------------------------------------------------------------
# Tree-sitter setup
# ---------------------------------------------------------------------------

CPP_LANG = Language(tscpp.language())
PY_LANG = Language(tspy.language())


def make_parser(lang):
    p = Parser(lang)
    return p


cpp_parser = make_parser(CPP_LANG)
py_parser = make_parser(PY_LANG)

# ---------------------------------------------------------------------------
# Symbol extraction
# ---------------------------------------------------------------------------


def extract_cpp_symbols(filepath, content):
    """Extract class/struct/function definitions from C++ code."""
    symbols = []
    tree = cpp_parser.parse(content.encode())

    def visit(node, namespace=""):
        if node.type == "namespace_definition":
            ns_name = ""
            for child in node.children:
                if child.type == "namespace_identifier":
                    ns_name = child.text.decode()
                    break
            body = None
            for child in node.children:
                if child.type == "declaration_list":
                    body = child
                    break
            if body:
                prefix = f"{namespace}{ns_name}::" if ns_name else namespace
                for child in body.children:
                    visit(child, prefix)
            return

        if node.type in ("class_specifier", "struct_specifier"):
            name_node = node.child_by_field_name("name")
            if name_node:
                name = name_node.text.decode()
                full_name = f"{namespace}{name}"
                # Extract method signatures
                methods = []
                body = node.child_by_field_name("body")
                if body:
                    for child in body.children:
                        if child.type == "function_definition":
                            sig = _cpp_func_signature(child)
                            if sig:
                                methods.append(sig)
                        elif child.type == "declaration":
                            # Could be a method declaration
                            sig = _cpp_decl_signature(child)
                            if sig:
                                methods.append(sig)
                kind = "class" if node.type == "class_specifier" else "struct"
                symbols.append({
                    "name": full_name,
                    "kind": kind,
                    "signature": f"{kind} {full_name}",
                    "methods": methods[:10],  # Cap to avoid bloat
                    "file": str(filepath),
                    "line": node.start_point[0] + 1,
                })

        elif node.type == "function_definition":
            sig = _cpp_func_signature(node)
            if sig:
                symbols.append({
                    "name": sig.split("(")[0].split()[-1] if "(" in sig else sig,
                    "kind": "function",
                    "signature": sig,
                    "methods": [],
                    "file": str(filepath),
                    "line": node.start_point[0] + 1,
                })

        # Visit children for top-level traversal
        for child in node.children:
            visit(child, namespace)

    for child in tree.root_node.children:
        visit(child)

    return symbols


def _cpp_func_signature(node):
    """Extract a concise function signature from a function_definition node."""
    declarator = node.child_by_field_name("declarator")
    if not declarator:
        return None
    # Get return type
    ret_type = ""
    for child in node.children:
        if child == declarator:
            break
        if child.type not in ("comment", "attribute_declaration"):
            ret_type += child.text.decode() + " "
    ret_type = ret_type.strip()
    decl_text = declarator.text.decode()
    # Truncate long signatures
    sig = f"{ret_type} {decl_text}".strip()
    if len(sig) > 120:
        sig = sig[:117] + "..."
    return sig


def _cpp_decl_signature(node):
    """Extract signature from a declaration that might be a method decl."""
    text = node.text.decode().strip()
    if "(" in text and ";" in text:
        sig = text.rstrip(";").strip()
        if len(sig) > 120:
            sig = sig[:117] + "..."
        return sig
    return None


def extract_py_symbols(filepath, content):
    """Extract class/function definitions from Python code."""
    symbols = []
    tree = py_parser.parse(content.encode())

    for node in tree.root_node.children:
        if node.type == "class_definition":
            name_node = node.child_by_field_name("name")
            if not name_node:
                continue
            name = name_node.text.decode()
            # Get superclasses
            superclass = ""
            for child in node.children:
                if child.type == "argument_list":
                    superclass = child.text.decode()
                    break
            # Get method names
            methods = []
            body = node.child_by_field_name("body")
            if body:
                for child in body.children:
                    if child.type == "function_definition":
                        mname = child.child_by_field_name("name")
                        if mname:
                            mtext = mname.text.decode()
                            if not mtext.startswith("_") or mtext in (
                                "__init__", "__repr__", "__iter__"
                            ):
                                params = child.child_by_field_name("parameters")
                                psig = params.text.decode() if params else "()"
                                methods.append(f"{mtext}{psig}")

            sig = f"class {name}{superclass}" if superclass else f"class {name}"
            symbols.append({
                "name": name,
                "kind": "class",
                "signature": sig,
                "methods": methods[:10],
                "file": str(filepath),
                "line": node.start_point[0] + 1,
            })

        elif node.type == "function_definition":
            name_node = node.child_by_field_name("name")
            if not name_node:
                continue
            name = name_node.text.decode()
            if name.startswith("_") and name != "__init__":
                continue
            params = node.child_by_field_name("parameters")
            psig = params.text.decode() if params else "()"
            sig = f"def {name}{psig}"
            if len(sig) > 120:
                sig = sig[:117] + "..."
            symbols.append({
                "name": name,
                "kind": "function",
                "signature": sig,
                "methods": [],
                "file": str(filepath),
                "line": node.start_point[0] + 1,
            })

    return symbols


# ---------------------------------------------------------------------------
# Reference graph and ranking
# ---------------------------------------------------------------------------


def build_reference_graph(all_symbols, file_contents):
    """Build a graph of cross-file symbol references.

    Returns a dict: symbol_name -> number of other files that reference it.
    Only counts classes, structs, and non-trivial functions.
    """
    # Filter out trivial/common names that would create noise
    TRIVIAL_NAMES = {
        "name", "type", "end", "begin", "size", "empty", "get", "set",
        "data", "value", "index", "clear", "push", "pop", "front", "back",
        "format", "write", "read", "to_string", "operator", "iterator",
        "const_iterator", "surface", "run", "reset", "init",
    }

    # Collect all symbol names and their source files
    symbol_files = {}  # name -> set of files where defined
    for sym in all_symbols:
        name = sym["name"]
        # Skip trivial names and very short names
        if name.lower() in TRIVIAL_NAMES or len(name) < 4:
            continue
        # Skip pure accessor patterns
        if sym["kind"] == "function" and sym["signature"]:
            sig_lower = sym["signature"].lower()
            if any(p in sig_lower for p in [
                "const {", "() const", "& name()", "type() const override"
            ]):
                # Only skip if it's a simple accessor (short signature)
                if len(sym["signature"]) < 60 and "(" in sig_lower:
                    parts = sig_lower.split("(")[0].split()
                    if parts and parts[-1] in TRIVIAL_NAMES:
                        continue
        symbol_files.setdefault(name, set()).add(sym["file"])

    # Count how many OTHER files reference each symbol
    ref_counts = defaultdict(int)
    for filepath, content in file_contents.items():
        for sym_name, def_files in symbol_files.items():
            if filepath not in def_files and sym_name in content:
                ref_counts[sym_name] += 1

    return ref_counts


def rank_symbols(all_symbols, ref_counts):
    """Rank symbols by cross-file reference count (simplified PageRank).

    Boost classes/structs since they represent key abstractions.
    """
    for sym in all_symbols:
        base_score = ref_counts.get(sym["name"], 0)
        # Boost classes/structs - they're more informative than individual functions
        if sym["kind"] in ("class", "struct"):
            base_score = int(base_score * 1.5) + 2
        sym["score"] = base_score
    return sorted(all_symbols, key=lambda s: (-s["score"], s["name"]))


# ---------------------------------------------------------------------------
# Subsystem categorization
# ---------------------------------------------------------------------------


def categorize_file(filepath):
    """Map a file path to a logical subsystem."""
    rel = str(filepath).replace("\\", "/").lower()
    for keyword, subsystem in SUBSYSTEM_RULES:
        if keyword in rel:
            return subsystem
    return "Other"


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------


def generate_repomap_md(ranked_symbols, max_lines=MAX_REPOMAP_LINES):
    """Generate concise markdown repo map."""
    # Deduplicate: keep highest-scored version of each name per subsystem
    seen = set()
    deduped = []
    for sym in ranked_symbols:
        subsystem = categorize_file(sym["file"])
        key = (subsystem, sym["name"])
        if key not in seen:
            seen.add(key)
            deduped.append(sym)

    # Group by subsystem
    groups = defaultdict(list)
    for sym in deduped:
        subsystem = categorize_file(sym["file"])
        if subsystem == "Other" and sym["score"] < 3:
            continue  # Skip low-value "Other" symbols
        groups[subsystem].append(sym)

    # Sort groups by max score in group, drop "Other" to the end
    sorted_groups = sorted(
        groups.items(),
        key=lambda g: (g[0] == "Other", -max(s["score"] for s in g[1]) if g[1] else 0),
    )

    lines = [
        "# OpenMC Repo Map",
        "",
        "Auto-generated structural overview. Top symbols ranked by cross-file usage.",
        "",
    ]

    for group_name, syms in sorted_groups:
        if len(lines) >= max_lines - 2:
            break

        lines.append(f"## {group_name}")

        # Show top symbols in this group
        shown = 0
        for sym in syms:
            if shown >= 5 or len(lines) >= max_lines - 1:
                break
            rel_file = os.path.relpath(sym["file"], OPENMC_ROOT)
            # Flatten signature to single line
            sig_flat = " ".join(sym["signature"].split())
            if len(sig_flat) > 80:
                sig_flat = sig_flat[:77] + "..."
            lines.append(f"- `{sig_flat}` ({rel_file}:{sym['line']})")

            # Show key methods for classes (max 2, single line each)
            if sym["kind"] in ("class", "struct") and sym["methods"]:
                for method in sym["methods"][:2]:
                    # Flatten to single line
                    method_flat = " ".join(method.split())
                    if len(method_flat) > 70:
                        method_flat = method_flat[:67] + "..."
                    lines.append(f"  - `{method_flat}`")
                    if len(lines) >= max_lines - 1:
                        break

            shown += 1

        lines.append("")

    return "\n".join(lines[:max_lines])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print(f"OpenMC root: {OPENMC_ROOT}")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Collect all source files
    file_contents = {}
    cpp_files = []
    py_files = []

    for pattern in CPP_PATTERNS:
        for fp in sorted(OPENMC_ROOT.glob(pattern)):
            try:
                content = fp.read_text(errors="replace")
                rel = str(fp.relative_to(OPENMC_ROOT))
                file_contents[rel] = content
                cpp_files.append((rel, content))
            except Exception as e:
                print(f"  Warning: could not read {fp}: {e}", file=sys.stderr)

    for pattern in PY_PATTERNS:
        for fp in sorted(OPENMC_ROOT.glob(pattern)):
            if "__pycache__" in str(fp):
                continue
            try:
                content = fp.read_text(errors="replace")
                rel = str(fp.relative_to(OPENMC_ROOT))
                file_contents[rel] = content
                py_files.append((rel, content))
            except Exception as e:
                print(f"  Warning: could not read {fp}: {e}", file=sys.stderr)

    print(f"Found {len(cpp_files)} C++/H files, {len(py_files)} Python files")

    # Extract symbols
    all_symbols = []
    for rel, content in cpp_files:
        syms = extract_cpp_symbols(rel, content)
        all_symbols.extend(syms)

    for rel, content in py_files:
        syms = extract_py_symbols(rel, content)
        all_symbols.extend(syms)

    print(f"Extracted {len(all_symbols)} symbols")

    # Build reference graph and rank
    ref_counts = build_reference_graph(all_symbols, file_contents)
    ranked = rank_symbols(all_symbols, ref_counts)

    # Generate outputs
    repomap_md = generate_repomap_md(ranked)
    repomap_path = CACHE_DIR / "repomap.md"
    repomap_path.write_text(repomap_md)
    print(f"Wrote {repomap_path} ({len(repomap_md.splitlines())} lines)")

    # Full JSON for other tools
    json_path = CACHE_DIR / "repomap_full.json"
    json_data = []
    for sym in ranked:
        json_data.append({
            "name": sym["name"],
            "kind": sym["kind"],
            "signature": sym["signature"],
            "methods": sym["methods"],
            "file": sym["file"],
            "line": sym["line"],
            "score": sym["score"],
        })
    json_path.write_text(json.dumps(json_data, indent=2))
    print(f"Wrote {json_path} ({len(json_data)} symbols)")


if __name__ == "__main__":
    main()
