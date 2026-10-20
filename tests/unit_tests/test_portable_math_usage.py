"""Guard against transcendental calls that bypass the portable math layer.

Transport results are only reproducible across platforms when every
transcendental function in the C++ sources goes through ``openmc/math.h``.
This test scans the C++ sources for direct ``std::`` calls to the functions
that header provides, so a new call site that would reintroduce a dependence
on the system libm fails here instead of surfacing later as a platform
difference in the regression tests.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIRS = [ROOT / 'src', ROOT / 'include' / 'openmc']

# Allowed to call the system libm: third-party code bundled verbatim and the
# header that maps the openmc:: names to a backend.
EXEMPT = [ROOT / 'src' / 'external', ROOT / 'include' / 'openmc' / 'external',
          ROOT / 'include' / 'openmc' / 'math.h']

FUNCTIONS = ('exp', 'expm1', 'log', 'log1p', 'log2', 'pow', 'sin', 'cos',
             'tan', 'asin', 'acos', 'atan', 'atan2', 'sinh', 'cosh', 'erf',
             'erfc', 'lgamma', 'tgamma', 'cbrt', 'hypot')

# A direct std:: call, an unqualified call (which resolves to the C library
# function through <cmath>), or a using-declaration that would make std::
# names visible. Calls qualified with openmc:: and member calls are fine.
PATTERN = re.compile(r'(?:(?<![\w.>:])|(?<=std::))(' + '|'.join(FUNCTIONS)
                     + r')\s*\(')
DECLARATION = re.compile(r'^\s*(inline\s+|static\s+|extern\s+)?'
                         r'(double|float|Tensor<T>)\s+('
                         + '|'.join(FUNCTIONS) + r')\s*\(')
USING = re.compile(r'\busing\s+(namespace\s+std\s*;|std::('
                   + '|'.join(FUNCTIONS) + r')\s*;)')


def _is_exempt(path):
    return any(path == e or e in path.parents for e in EXEMPT)


def _sources():
    for d in SOURCE_DIRS:
        for suffix in ('*.cpp', '*.h'):
            for path in sorted(d.rglob(suffix)):
                if not _is_exempt(path):
                    yield path


def test_transcendentals_use_portable_layer():
    if not all(d.is_dir() for d in SOURCE_DIRS):
        pytest.skip('C++ sources are not available in this checkout')

    offenders = []
    for path in _sources():
        in_block = False
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            code = line
            if in_block:
                end = code.find('*/')
                if end < 0:
                    continue
                code = code[end + 2:]
                in_block = False
            code = re.sub(r'/\*.*?\*/', '', code)
            start = code.find('/*')
            if start >= 0:
                code, in_block = code[:start], True
            code = code.split('//')[0]
            if USING.search(code) or (
                    PATTERN.search(code) and not DECLARATION.match(code)):
                offenders.append(f'{path.relative_to(ROOT)}:{lineno}: '
                                 f'{line.strip()}')

    assert not offenders, (
        'Transcendental functions must be called through openmc::math '
        '(include "openmc/math.h") so results are reproducible across '
        'platforms:\n' + '\n'.join(offenders))
