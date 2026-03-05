#!/usr/bin/env python3
"""LSP-based code navigation for OpenMC using clangd.

Uses the Language Server Protocol to provide compiler-accurate symbol
resolution, go-to-definition, find-references, and related-file discovery.
Unlike tree-sitter-based tools, this resolves symbols through the actual
C++ type system — no false edges from name collisions.

Requires:
  - clangd (apt-get install clangd, or clangd-15/clangd-16/etc.)
  - compile_commands.json (cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON)

Usage:
    openmc_lsp.py symbols src/simulation.cpp
    openmc_lsp.py definition src/simulation.cpp:132
    openmc_lsp.py references src/simulation.cpp:132
    openmc_lsp.py related src/simulation.cpp
    openmc_lsp.py related src/simulation.cpp --top-k 20

Examples:
    openmc_lsp.py symbols src/particle.cpp
    openmc_lsp.py definition src/simulation.cpp:132    # where is write_message defined?
    openmc_lsp.py references include/openmc/error.h:55 # who calls write_message?
    openmc_lsp.py related src/simulation.cpp            # files connected by real references
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
from collections import Counter, defaultdict
from pathlib import Path

OPENMC_ROOT = Path(__file__).resolve().parents[3]

# Symbol kind names (LSP spec)
SYMBOL_KINDS = {
    1: "File", 2: "Module", 3: "Namespace", 4: "Package", 5: "Class",
    6: "Method", 7: "Property", 8: "Field", 9: "Constructor", 10: "Enum",
    11: "Interface", 12: "Function", 13: "Variable", 14: "Constant",
    15: "String", 16: "Number", 17: "Boolean", 18: "Array", 19: "Object",
    20: "Key", 21: "Null", 22: "EnumMember", 23: "Struct", 24: "Event",
    25: "Operator", 26: "TypeParameter",
}


class ClangdClient:
    """Minimal LSP client that talks to clangd via JSON-RPC over stdin/stdout."""

    def __init__(self, compile_commands_dir=None):
        clangd = self._find_clangd()
        if not clangd:
            print("ERROR: clangd not found. Install with: apt-get install clangd",
                  file=sys.stderr)
            sys.exit(1)

        if not compile_commands_dir:
            compile_commands_dir = self._find_compile_commands()
        if not compile_commands_dir:
            print("ERROR: compile_commands.json not found. Generate with:\n"
                  "  cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
                  file=sys.stderr)
            sys.exit(1)

        args = [clangd, '--compile-commands-dir=' + str(compile_commands_dir)]
        self.proc = subprocess.Popen(
            args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        self._id = 0
        self._opened_files = set()
        self._initialize()

    def _find_clangd(self):
        """Find clangd binary, trying common names."""
        for name in ['clangd', 'clangd-15', 'clangd-16', 'clangd-17', 'clangd-18']:
            path = shutil.which(name)
            if path:
                return path
        return None

    def _find_compile_commands(self):
        """Find compile_commands.json in common locations."""
        for d in [OPENMC_ROOT / 'build', OPENMC_ROOT]:
            if (d / 'compile_commands.json').exists():
                return str(d)
        return None

    def _initialize(self):
        """Send LSP initialize/initialized handshake."""
        self.request("initialize", {
            "processId": os.getpid(),
            "rootUri": OPENMC_ROOT.as_uri(),
            "capabilities": {}
        })
        self.notify("initialized")

    def _next_id(self):
        self._id += 1
        return self._id

    def _send(self, msg_dict):
        body = json.dumps(msg_dict)
        encoded = body.encode('utf-8')
        header = f"Content-Length: {len(encoded)}\r\n\r\n"
        self.proc.stdin.write(header.encode('ascii') + encoded)
        self.proc.stdin.flush()

    def request(self, method, params=None):
        """Send a request and wait for the response."""
        rid = self._next_id()
        self._send({"jsonrpc": "2.0", "id": rid, "method": method,
                     "params": params or {}})
        while True:
            msg = self._read_msg()
            if msg.get('id') == rid:
                if 'error' in msg:
                    return None
                return msg.get('result')

    def notify(self, method, params=None):
        """Send a notification (no response expected)."""
        self._send({"jsonrpc": "2.0", "method": method,
                     "params": params or {}})

    def _read_msg(self):
        headers = {}
        while True:
            line = self.proc.stdout.readline()
            if not line:
                raise EOFError("clangd process terminated")
            line = line.decode('utf-8').strip()
            if not line:
                break
            k, v = line.split(': ', 1)
            headers[k] = v
        length = int(headers['Content-Length'])
        body = self.proc.stdout.read(length)
        return json.loads(body)

    def open_file(self, filepath):
        """Open a file in clangd and wait for it to be indexed."""
        fpath = Path(filepath)
        if not fpath.is_absolute():
            fpath = OPENMC_ROOT / fpath
        uri = fpath.as_uri()
        if uri in self._opened_files:
            return uri
        text = fpath.read_text()
        self.notify("textDocument/didOpen", {
            "textDocument": {
                "uri": uri, "languageId": "cpp", "version": 1, "text": text
            }
        })
        self._opened_files.add(uri)
        # Give clangd time to parse. First file takes longer (preamble build).
        wait = 8 if len(self._opened_files) == 1 else 3
        time.sleep(wait)
        return uri

    def get_symbols(self, filepath):
        """Get all symbols defined in a file."""
        uri = self.open_file(filepath)
        result = self.request("textDocument/documentSymbol", {
            "textDocument": {"uri": uri}
        })
        return result or []

    def get_definition(self, filepath, line, character):
        """Get definition location for symbol at position."""
        uri = self.open_file(filepath)
        result = self.request("textDocument/definition", {
            "textDocument": {"uri": uri},
            "position": {"line": line, "character": character}
        })
        return result or []

    def get_references(self, filepath, line, character,
                       include_declaration=True):
        """Get all references to symbol at position."""
        uri = self.open_file(filepath)
        result = self.request("textDocument/references", {
            "textDocument": {"uri": uri},
            "position": {"line": line, "character": character},
            "context": {"includeDeclaration": include_declaration}
        })
        return result or []

    def close(self):
        """Shutdown clangd cleanly."""
        try:
            self.request("shutdown")
            self.notify("exit")
            self.proc.wait(timeout=5)
        except Exception:
            self.proc.kill()


def uri_to_relpath(uri):
    """Convert file:// URI to path relative to OPENMC_ROOT."""
    path = urllib.parse.unquote(uri.replace('file://', ''))
    try:
        return str(Path(path).relative_to(OPENMC_ROOT))
    except ValueError:
        return path


def is_project_file(relpath):
    """Check if a path is an OpenMC project file (not system/vendor)."""
    if relpath.startswith('/'):
        return False  # absolute path = system header
    if relpath.startswith('vendor/'):
        return False
    return True


def get_symbol_range(sym):
    """Extract start line/character from either SymbolInformation or DocumentSymbol."""
    # DocumentSymbol format: has 'range' and 'selectionRange' at top level
    if 'selectionRange' in sym:
        return sym['selectionRange']['start']
    # DocumentSymbol without selectionRange
    if 'range' in sym and isinstance(sym['range'], dict) and 'start' in sym['range']:
        return sym['range']['start']
    # SymbolInformation format: has 'location.range'
    if 'location' in sym:
        return sym['location']['range']['start']
    return {'line': 0, 'character': 0}


def flatten_symbols(symbols, depth=0):
    """Flatten nested document symbols into a flat list with depth info."""
    result = []
    for s in symbols:
        result.append((s, depth))
        children = s.get('children', [])
        if children:
            result.extend(flatten_symbols(children, depth + 1))
    return result


def cmd_symbols(client, filepath):
    """List all symbols defined in a file."""
    symbols = client.get_symbols(filepath)
    flat = flatten_symbols(symbols)
    for sym, depth in flat:
        kind_name = SYMBOL_KINDS.get(sym['kind'], f"kind={sym['kind']}")
        start = get_symbol_range(sym)
        line = start['line']
        indent = "  " * depth
        print(f"{indent}{kind_name}: {sym['name']} (line {line + 1})")


def cmd_definition(client, filepath, line, character=None):
    """Find the definition of a symbol."""
    if character is None:
        # Find first non-whitespace identifier on the line
        fpath = Path(filepath)
        if not fpath.is_absolute():
            fpath = OPENMC_ROOT / fpath
        lines = fpath.read_text().split('\n')
        if line - 1 < len(lines):
            text = lines[line - 1]
            # Skip leading whitespace
            character = len(text) - len(text.lstrip())

    result = client.get_definition(filepath, line - 1, character)
    if not result:
        print("No definition found.")
        return

    if isinstance(result, dict):
        result = [result]
    for loc in result:
        rel = uri_to_relpath(loc['uri'])
        ln = loc['range']['start']['line'] + 1
        print(f"  {rel}:{ln}")


def cmd_references(client, filepath, line, character=None):
    """Find all references to a symbol."""
    if character is None:
        fpath = Path(filepath)
        if not fpath.is_absolute():
            fpath = OPENMC_ROOT / fpath
        lines = fpath.read_text().split('\n')
        if line - 1 < len(lines):
            text = lines[line - 1]
            character = len(text) - len(text.lstrip())

    result = client.get_references(filepath, line - 1, character)
    if not result:
        print("No references found.")
        return

    # Group by file
    by_file = defaultdict(list)
    for loc in result:
        rel = uri_to_relpath(loc['uri'])
        ln = loc['range']['start']['line'] + 1
        by_file[rel].append(ln)

    print(f"{len(result)} references across {len(by_file)} files:\n")
    for fpath, lines_list in sorted(by_file.items()):
        lines_str = ", ".join(str(l) for l in sorted(lines_list))
        print(f"  {fpath}:{lines_str}")


def cmd_related(client, filepath, top_k=15):
    """Find files related to a given file through real typed references.

    For each symbol defined in the target file, finds all files that
    reference it. Returns files ranked by connection count.
    """
    symbols = client.get_symbols(filepath)
    flat = flatten_symbols(symbols)

    # Filter to meaningful symbols (functions, classes, methods, variables)
    interesting_kinds = {5, 6, 8, 12, 13, 23}  # Class, Method, Field, Function, Variable, Struct
    interesting = [(s, d) for s, d in flat if s['kind'] in interesting_kinds]

    if not interesting:
        print("No interesting symbols found in file.")
        return

    target_rel = filepath
    if Path(filepath).is_absolute():
        target_rel = str(Path(filepath).relative_to(OPENMC_ROOT))

    file_connections = Counter()  # file -> number of symbols referencing it
    symbol_details = defaultdict(set)  # file -> set of symbol names

    print(f"Analyzing {len(interesting)} symbols in {target_rel}...\n",
          file=sys.stderr)

    # Read the file so we can find exact symbol name positions
    fpath_obj = Path(filepath)
    if not fpath_obj.is_absolute():
        fpath_obj = OPENMC_ROOT / fpath_obj
    file_lines = fpath_obj.read_text().split('\n')

    for sym, depth in interesting:
        start = get_symbol_range(sym)
        line = start['line']
        char = start['character']

        # The range start may point to the type, not the symbol name.
        # Find the actual symbol name position within the line.
        if line < len(file_lines):
            name_col = file_lines[line].find(sym['name'], char)
            if name_col >= 0:
                char = name_col

        refs = client.get_references(filepath, line, char,
                                     include_declaration=False)
        if not refs:
            continue

        for loc in refs:
            rel = uri_to_relpath(loc['uri'])
            if rel == target_rel:
                continue
            if not is_project_file(rel):
                continue
            file_connections[rel] += 1
            symbol_details[rel].add(sym['name'])

    if not file_connections:
        print("No external references found.")
        return

    print(f"Files related to {target_rel} "
          f"(ranked by typed reference count):\n")
    for fpath, count in file_connections.most_common(top_k):
        syms = sorted(symbol_details[fpath])
        sym_preview = ", ".join(syms[:5])
        if len(syms) > 5:
            sym_preview += f", ... (+{len(syms)-5} more)"
        print(f"  [{count:3d} refs] {fpath}")
        print(f"           via: {sym_preview}")


def parse_file_location(location):
    """Parse 'filepath:line' into (filepath, line) or (filepath, None)."""
    # Handle filepath:line format
    parts = location.rsplit(':', 1)
    if len(parts) == 2:
        try:
            line = int(parts[1])
            return parts[0], line
        except ValueError:
            pass
    return location, None


def main():
    parser = argparse.ArgumentParser(
        description="LSP-based code navigation for OpenMC (via clangd)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  %(prog)s symbols src/simulation.cpp
  %(prog)s definition src/simulation.cpp:132
  %(prog)s references src/simulation.cpp:132
  %(prog)s related src/simulation.cpp
  %(prog)s related src/simulation.cpp --top-k 20""",
    )
    parser.add_argument("command",
                        choices=["symbols", "definition", "references",
                                 "related"],
                        help="Command to run")
    parser.add_argument("location",
                        help="File path, or file:line for definition/references")
    parser.add_argument("--top-k", type=int, default=15,
                        help="Number of related files to show (default: 15)")
    parser.add_argument("--compile-commands-dir", type=str, default=None,
                        help="Directory containing compile_commands.json")

    args = parser.parse_args()

    filepath, line = parse_file_location(args.location)

    # Validate file exists
    fpath = Path(filepath)
    if not fpath.is_absolute():
        fpath = OPENMC_ROOT / fpath
    if not fpath.exists():
        print(f"ERROR: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)

    client = ClangdClient(compile_commands_dir=args.compile_commands_dir)
    try:
        if args.command == "symbols":
            cmd_symbols(client, filepath)
        elif args.command == "definition":
            if line is None:
                print("ERROR: definition requires file:line format",
                      file=sys.stderr)
                sys.exit(1)
            cmd_definition(client, filepath, line)
        elif args.command == "references":
            if line is None:
                print("ERROR: references requires file:line format",
                      file=sys.stderr)
                sys.exit(1)
            cmd_references(client, filepath, line)
        elif args.command == "related":
            cmd_related(client, filepath, top_k=args.top_k)
    finally:
        client.close()


if __name__ == "__main__":
    main()
