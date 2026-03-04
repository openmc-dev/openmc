"""Chunk OpenMC source files and documentation for RAG indexing.

Code files are chunked at the function/class level using tree-sitter.
RST documentation is chunked by section headers.
"""

import re
from pathlib import Path

import tree_sitter_cpp as tscpp
import tree_sitter_python as tspy
from tree_sitter import Language, Parser

CPP_LANG = Language(tscpp.language())
PY_LANG = Language(tspy.language())

cpp_parser = Parser(CPP_LANG)
py_parser = Parser(PY_LANG)

MAX_CHUNK_CHARS = 1500
MIN_CHUNK_CHARS = 50


def chunk_file(filepath, openmc_root):
    """Chunk a single file based on its extension."""
    filepath = Path(filepath)
    rel = str(filepath.relative_to(openmc_root))
    try:
        content = filepath.read_text(errors="replace")
    except Exception:
        return []

    if filepath.suffix in (".cpp", ".h"):
        return _chunk_cpp(rel, content)
    elif filepath.suffix == ".py":
        return _chunk_python(rel, content)
    elif filepath.suffix == ".rst":
        return _chunk_rst(rel, content)
    return []


def _chunk_cpp(rel_path, content):
    """Extract function and class-level chunks from C++ code."""
    tree = cpp_parser.parse(content.encode())
    chunks = []
    used_ranges = []

    def _extract_node(node, kind_override=None):
        text = content[node.start_byte:node.end_byte]
        if len(text) < MIN_CHUNK_CHARS:
            return
        # Extract symbol name
        name = _get_node_name(node)
        kind = kind_override or node.type
        for sub in _split_if_large(text):
            chunks.append({
                "text": sub,
                "filepath": rel_path,
                "kind": kind,
                "symbol": name or "",
                "start_line": node.start_point[0] + 1,
                "end_line": node.end_point[0] + 1,
            })
        used_ranges.append((node.start_byte, node.end_byte))

    def _visit(node):
        if node.type in (
            "function_definition", "class_specifier",
            "struct_specifier", "enum_specifier",
        ):
            _extract_node(node)
        elif node.type == "namespace_definition":
            # Visit children inside namespaces
            for child in node.children:
                _visit(child)
        elif node.type == "declaration_list":
            for child in node.children:
                _visit(child)
        else:
            for child in node.children:
                if child.type in (
                    "function_definition", "class_specifier",
                    "struct_specifier", "namespace_definition",
                ):
                    _visit(child)

    for child in tree.root_node.children:
        _visit(child)

    # Add file header (includes, forward declarations) as a separate chunk
    header_lines = []
    for line in content.split("\n")[:50]:
        if line.strip().startswith("#include") or line.strip().startswith("namespace") \
           or line.strip().startswith("//") or line.strip().startswith("using") \
           or line.strip() == "":
            header_lines.append(line)
        else:
            break
    header = "\n".join(header_lines).strip()
    if len(header) >= MIN_CHUNK_CHARS:
        chunks.append({
            "text": header,
            "filepath": rel_path,
            "kind": "file_header",
            "symbol": Path(rel_path).name,
            "start_line": 1,
            "end_line": len(header_lines),
        })

    return chunks


def _chunk_python(rel_path, content):
    """Extract function and class-level chunks from Python code."""
    tree = py_parser.parse(content.encode())
    chunks = []

    for node in tree.root_node.children:
        if node.type in ("class_definition", "function_definition"):
            text = content[node.start_byte:node.end_byte]
            if len(text) < MIN_CHUNK_CHARS:
                continue
            name_node = node.child_by_field_name("name")
            name = name_node.text.decode() if name_node else ""
            for sub in _split_if_large(text):
                chunks.append({
                    "text": sub,
                    "filepath": rel_path,
                    "kind": node.type.replace("_definition", ""),
                    "symbol": name,
                    "start_line": node.start_point[0] + 1,
                    "end_line": node.end_point[0] + 1,
                })

    # Module-level docstring + imports as header
    header_lines = []
    for line in content.split("\n")[:40]:
        stripped = line.strip()
        if stripped.startswith(("import ", "from ", "#", '"""', "'''", "")) \
           or stripped == "":
            header_lines.append(line)
        elif stripped.startswith(("def ", "class ")):
            break
        else:
            header_lines.append(line)
    header = "\n".join(header_lines).strip()
    if len(header) >= MIN_CHUNK_CHARS:
        chunks.append({
            "text": header,
            "filepath": rel_path,
            "kind": "file_header",
            "symbol": Path(rel_path).name,
            "start_line": 1,
            "end_line": len(header_lines),
        })

    return chunks


def _chunk_rst(rel_path, content):
    """Chunk RST documentation by section headers."""
    # RST sections are indicated by underlines of =, -, ~, ^, etc.
    section_pattern = re.compile(
        r'^(.+)\n([=\-~^"+]+)\s*$', re.MULTILINE
    )
    chunks = []

    # Find all section positions
    positions = [0]
    for m in section_pattern.finditer(content):
        # The section title starts at the beginning of the title line
        positions.append(m.start())
    positions.append(len(content))

    for i in range(len(positions) - 1):
        section = content[positions[i]:positions[i + 1]].strip()
        if len(section) < MIN_CHUNK_CHARS:
            continue
        # Extract title
        title_match = section_pattern.match(section)
        title = title_match.group(1).strip() if title_match else ""
        start_line = content[:positions[i]].count("\n") + 1
        end_line = content[:positions[i + 1]].count("\n") + 1
        for sub in _split_if_large(section):
            chunks.append({
                "text": sub,
                "filepath": rel_path,
                "kind": "doc_section",
                "symbol": title,
                "start_line": start_line,
                "end_line": end_line,
            })

    return chunks


def _get_node_name(node):
    """Extract the name from a tree-sitter node."""
    name_node = node.child_by_field_name("name")
    if name_node:
        return name_node.text.decode()
    # For function_definition, check declarator
    decl = node.child_by_field_name("declarator")
    if decl:
        # Walk down to find the identifier
        while decl.type not in ("identifier", "qualified_identifier",
                                "field_identifier", "destructor_name"):
            found = False
            for child in decl.children:
                if child.type in ("function_declarator", "identifier",
                                  "qualified_identifier", "field_identifier",
                                  "destructor_name", "template_function"):
                    decl = child
                    found = True
                    break
            if not found:
                break
        return decl.text.decode()
    return ""


def _split_if_large(text, max_chars=MAX_CHUNK_CHARS):
    """Split text into chunks if it exceeds max_chars."""
    if len(text) <= max_chars:
        return [text]
    # Split on line boundaries
    lines = text.split("\n")
    chunks = []
    current = []
    current_len = 0
    for line in lines:
        if current_len + len(line) + 1 > max_chars and current:
            chunks.append("\n".join(current))
            current = [line]
            current_len = len(line)
        else:
            current.append(line)
            current_len += len(line) + 1
    if current:
        chunks.append("\n".join(current))
    return chunks
