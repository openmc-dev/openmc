#!/usr/bin/env python3
"""Accuracy report for the transcendental functions behind openmc/math.h.

Reads the output of tools/dev/portable_math_sample.cpp and compares both the
openmc:: results (the CORE-MATH backend) and the system libm results against
correctly rounded references computed with MPFR through gmpy2. For every function the
report gives the maximum and mean error in units in the last place (ulp), the
fraction of correctly rounded results, and any disagreement on NaN or infinite
results.

Usage:

    python tools/dev/portable_math_report.py samples.txt [--json report.json]
        [--markdown report.md] [--kat tests/cpp_unit_tests/portable_math_kat.h]

The --kat option writes the special values (all of them for unary functions,
a strided subset of the pairs for binary ones) and the first ordinary sampled
inputs of each function, together with the result bits, as a C++ header of
known-answer vectors. Every vector is checked to be correctly rounded first. The
C++ unit test in tests/cpp_unit_tests/test_portable_math.cpp checks those bits
on every platform, so a platform that produced different bits would fail.

Requires gmpy2 (pip install gmpy2).
"""

import argparse
import json
import math
import struct
import sys
from collections import defaultdict

try:
    import gmpy2
    from gmpy2 import mpfr
except ImportError:  # pragma: no cover
    sys.exit("gmpy2 is required: pip install gmpy2")

# A context that reproduces IEEE 754 binary64 arithmetic, including subnormal
# results and overflow to infinity, with round-to-nearest-even.
CTX = gmpy2.ieee(64)
CTX.trap_invalid = False
CTX.trap_divzero = False
CTX.trap_overflow = False
CTX.trap_underflow = False
CTX.trap_inexact = False
CTX.trap_erange = False
gmpy2.set_context(CTX)


def lgamma_ref(x):
    value, _sign = gmpy2.lgamma(x)
    return value


REFERENCE = {
    'exp': gmpy2.exp,
    'expm1': gmpy2.expm1,
    'log': gmpy2.log,
    'log1p': gmpy2.log1p,
    'log2': gmpy2.log2,
    'sin': gmpy2.sin,
    'cos': gmpy2.cos,
    'tan': gmpy2.tan,
    'asin': gmpy2.asin,
    'acos': gmpy2.acos,
    'atan': gmpy2.atan,
    'sinh': gmpy2.sinh,
    'cosh': gmpy2.cosh,
    'erf': gmpy2.erf,
    'erfc': gmpy2.erfc,
    'lgamma': lgamma_ref,
    'tgamma': gmpy2.gamma,
    'cbrt': gmpy2.cbrt,
    'pow': lambda x, y: x ** y,
    'atan2': gmpy2.atan2,
    'hypot': gmpy2.hypot,
}

BINARY = {'pow', 'atan2', 'hypot'}


def bits_to_double(word):
    return struct.unpack('<d', struct.pack('<Q', int(word, 16)))[0]


def double_to_mpfr(x):
    # mpfr(float) is exact for every finite double, including subnormals
    return mpfr(x)


def mpfr_to_double(v):
    # The ieee(64) context has already rounded v to a double
    return float(v)


def ordered(x):
    """Map a finite double to an integer so that consecutive doubles differ
    by one and the ordering matches the real ordering."""
    u = struct.unpack('<q', struct.pack('<d', x))[0]
    return u if u >= 0 else -(u & 0x7FFFFFFFFFFFFFFF)


def ulp_distance(a, b):
    """Return the number of doubles between a and b, or None when the two are
    not comparable (one is NaN or infinite and the other is not)."""
    if math.isnan(a) or math.isnan(b):
        return 0 if (math.isnan(a) and math.isnan(b)) else None
    if math.isinf(a) or math.isinf(b):
        return 0 if a == b else None
    return abs(ordered(a) - ordered(b))


class Stats:
    def __init__(self):
        self.n = 0
        self.max_ulp = 0
        self.sum_ulp = 0
        self.exact = 0
        self.over_one = 0
        self.incomparable = 0
        self.worst = None

    def add(self, inputs, value, reference):
        self.n += 1
        d = ulp_distance(value, reference)
        if d is None:
            self.incomparable += 1
            if self.worst is None or self.worst[0] != 'incomparable':
                self.worst = ('incomparable', inputs, value, reference)
            return
        self.sum_ulp += d
        if d == 0:
            self.exact += 1
        if d > 1:
            self.over_one += 1
        if d > self.max_ulp:
            self.max_ulp = d
            self.worst = (d, inputs, value, reference)

    def as_dict(self):
        return {
            'points': self.n,
            'max_ulp': self.max_ulp,
            'mean_ulp': self.sum_ulp / self.n if self.n else 0.0,
            'correctly_rounded_fraction': self.exact / self.n if self.n else 0.0,
            'over_one_ulp': self.over_one,
            'nan_or_inf_disagreements': self.incomparable,
            'worst': None if self.worst is None else {
                'ulp': self.worst[0],
                'inputs': [repr(v) for v in self.worst[1]],
                'value': repr(self.worst[2]),
                'reference': repr(self.worst[3]),
            },
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('samples', help="sample file, '-' for stdin, or a JSON report written earlier with --json")
    parser.add_argument('--json')
    parser.add_argument('--markdown')
    parser.add_argument('--rst', help='write the table as a reStructuredText list-table')
    parser.add_argument('--kat', help='write known-answer C++ header')
    parser.add_argument('--kat-count', type=int, default=96,
                        help='ordinary sampled inputs per function in the '
                             'known-answer header')
    parser.add_argument('--kat-special-count', type=int, default=128,
                        help='special values or special-value pairs per '
                             'function in the known-answer header')
    args = parser.parse_args()

    ours = defaultdict(Stats)
    system = defaultdict(Stats)
    kat = defaultdict(list)
    kat_special = defaultdict(list)
    seen = defaultdict(int)
    specials = None
    order = []

    if args.samples.endswith('.json'):
        with open(args.samples) as fh:
            report = json.load(fh)
        order = list(report)
        if args.kat:
            sys.exit('--kat needs a sample file, not a JSON report')
        write_outputs(args, report, order, kat, kat_special)
        return

    fh = sys.stdin if args.samples == '-' else open(args.samples)
    with fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == '#':
                # Header written by the sampler: '# specials <count>'
                if len(parts) >= 3 and parts[1] == 'specials':
                    specials = [None] * int(parts[2])
                continue
            if specials is None:
                sys.exit('sample file lacks the "# specials" header')
            name = parts[0]
            if name not in REFERENCE:
                sys.exit(f'unknown function in samples: {name}')
            if name not in order:
                order.append(name)
            if name in BINARY:
                x, y, r_ours, r_sys = (bits_to_double(p) for p in parts[1:5])
                inputs = (x, y)
                ref = REFERENCE[name](double_to_mpfr(x), double_to_mpfr(y))
            else:
                x, r_ours, r_sys = (bits_to_double(p) for p in parts[1:4])
                inputs = (x,)
                ref = REFERENCE[name](double_to_mpfr(x))
            ref = mpfr_to_double(ref)
            ours[name].add(inputs, r_ours, ref)
            system[name].add(inputs, r_sys, ref)
            if args.kat:
                if ulp_distance(r_ours, ref) != 0:
                    sys.exit(f'{name}{inputs}: result {r_ours!r} is not the '
                             f'correctly rounded {ref!r}; refusing to write '
                             'known-answer vectors')
                seen[name] += 1
                entry = (parts[1:-2], parts[-2])
                n_special = (len(specials) ** 2 if name in BINARY
                             else len(specials))
                if seen[name] <= n_special:
                    # Keep every special value for unary functions and a
                    # strided subset of the special pairs for binary ones
                    stride = max(1, n_special // args.kat_special_count)
                    if (seen[name] - 1) % stride == 0 and \
                            len(kat_special[name]) < args.kat_special_count:
                        kat_special[name].append(entry)
                elif len(kat[name]) < args.kat_count:
                    kat[name].append(entry)

    report = {name: {'portable': ours[name].as_dict(),
                     'system_libm': system[name].as_dict()} for name in order}
    write_outputs(args, report, order, kat, kat_special)


def write_outputs(args, report, order, kat, kat_special):
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(report, fh, indent=2)

    lines = ['| Function | Points | Portable max ulp | Portable mean ulp | '
             'Portable correctly rounded | System max ulp | System mean ulp | '
             'System correctly rounded | NaN/Inf disagreements (portable/system) |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for name in order:
        p = report[name]['portable']
        s = report[name]['system_libm']
        lines.append(
            f"| {name} | {p['points']} | {p['max_ulp']} | {p['mean_ulp']:.4f} | "
            f"{100 * p['correctly_rounded_fraction']:.3f}% | {s['max_ulp']} | "
            f"{s['mean_ulp']:.4f} | {100 * s['correctly_rounded_fraction']:.3f}% | "
            f"{p['nan_or_inf_disagreements']}/{s['nan_or_inf_disagreements']} |")
    table = '\n'.join(lines)
    print(table)
    for name in order:
        w = report[name]['portable']['worst']
        if w and w['ulp'] not in (0, 1):
            print(f"worst portable {name}: {w['ulp']} ulp at {w['inputs']} -> "
                  f"{w['value']} (reference {w['reference']})")
    if args.markdown:
        with open(args.markdown, 'w') as fh:
            fh.write(table + '\n')

    if args.rst:
        with open(args.rst, 'w') as fh:
            fh.write('.. Generated by tools/dev/portable_math_report.py --rst; '
                     'do not edit.\n\n')
            fh.write('.. list-table::\n   :header-rows: 1\n\n')
            fh.write('   * - Function\n     - Points\n'
                     '     - Portable max ulp\n     - Portable mean ulp\n'
                     '     - Portable correctly rounded\n'
                     '     - System max ulp\n     - System mean ulp\n'
                     '     - System correctly rounded\n')
            for name in order:
                p = report[name]['portable']
                s = report[name]['system_libm']
                fh.write(f"   * - {name}\n     - {p['points']}\n"
                         f"     - {p['max_ulp']}\n     - {p['mean_ulp']:.4f}\n"
                         f"     - {100 * p['correctly_rounded_fraction']:.3f}%\n"
                         f"     - {s['max_ulp']}\n     - {s['mean_ulp']:.4f}\n"
                         f"     - {100 * s['correctly_rounded_fraction']:.3f}%\n")

    if args.kat:
        with open(args.kat, 'w') as fh:
            fh.write('// Generated by tools/dev/portable_math_report.py --kat.\n'
                     '// Known-answer vectors: input bits and the bits the\n'
                     '// portable functions returned on the reference platform.\n'
                     '// Do not edit by hand.\n\n'
                     '#ifndef OPENMC_PORTABLE_MATH_KAT_H\n'
                     '#define OPENMC_PORTABLE_MATH_KAT_H\n\n'
                     '#include <cstdint>\n\n'
                     'namespace portable_math_kat {\n\n'
                     'struct Unary {\n  const char* name;\n  uint64_t x;\n'
                     '  uint64_t result;\n};\n\n'
                     'struct Binary {\n  const char* name;\n  uint64_t x;\n'
                     '  uint64_t y;\n  uint64_t result;\n};\n\n')
            fh.write('const Unary unary[] = {\n')
            for name in order:
                if name in BINARY:
                    continue
                for inputs, result in kat_special[name] + kat[name]:
                    fh.write(f'  {{"{name}", 0x{inputs[0]}ull, 0x{result}ull}},\n')
            fh.write('};\n\nconst Binary binary[] = {\n')
            for name in order:
                if name not in BINARY:
                    continue
                for inputs, result in kat_special[name] + kat[name]:
                    # Wrapped at 80 columns as clang-format lays it out
                    head = f'  {{"{name}", 0x{inputs[0]}ull, 0x{inputs[1]}ull,'
                    tail = f'0x{result}ull}},'
                    if len(head) + 1 + len(tail) <= 80:
                        fh.write(f'{head} {tail}\n')
                    else:
                        fh.write(f'{head}\n    {tail}\n')
            fh.write('};\n\n} // namespace portable_math_kat\n\n'
                     '#endif // OPENMC_PORTABLE_MATH_KAT_H\n')


if __name__ == '__main__':
    main()
