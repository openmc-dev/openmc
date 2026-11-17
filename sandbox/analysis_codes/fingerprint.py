"""Fingerprint the shared analysis implementation for provenance."""
from pathlib import Path
import hashlib


def implementation_sha():
    root = Path(__file__).resolve().parent
    return hashlib.sha256(b''.join(p.name.encode()+p.read_bytes()
        for p in sorted(root.glob('*.py')))).hexdigest()
