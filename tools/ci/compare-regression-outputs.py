#!/usr/bin/env python3
"""Compare regression test results produced on different platforms.

Each platform runs the regression suite with ``--junitxml`` and leaves a
``results_error.dat`` next to every test whose results differ from the
committed reference. This script takes one directory per platform (the layout
produced by the reproducibility workflow: ``<dir>/junit.xml`` and
``<dir>/tests/regression_tests/<test>/results_error.dat``) and checks:

1. the set of tests that fail or error is the same on every platform;
2. for every test with a ``results_error.dat``, the output is byte-identical
   on every platform;
3. failures that leave no output to compare are reported. A test whose input
   files differ from the reference (generated on the Python side) or whose
   reference comparison does not keep an output file is reported as such; an
   exception or a crash fails the comparison, since it would otherwise go
   unseen.

Together these show that the platforms agree with each other even where they
disagree with references generated with a particular system math library.

Usage: compare-regression-outputs.py --json summary.json DIR [DIR ...]
"""

import argparse
import hashlib
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

INPUT_MISMATCH = 'Input files are broken'


def collect_outputs(root):
    files = {}
    for path in sorted(Path(root).rglob('results_error.dat')):
        test = str(path.relative_to(root).parent)
        files[test] = hashlib.sha256(path.read_bytes()).hexdigest()
    return files


def collect_junit(root):
    """Return {test id: message} for every failed or errored test case."""
    path = Path(root) / 'junit.xml'
    if not path.is_file():
        sys.exit(f'{path} not found: run pytest with --junitxml')
    failed = {}
    for case in ET.parse(path).getroot().iter('testcase'):
        for kind in ('failure', 'error'):
            node = case.find(kind)
            if node is not None:
                test_id = f"{case.get('classname')}::{case.get('name')}"
                failed[test_id] = (kind, node.get('message', '') or '')
    return failed


def test_dir(test_id):
    """Map a junit test id to the test directory under tests/regression_tests."""
    classname = test_id.split('::')[0]
    parts = classname.split('.')
    if parts[-1] == 'test':
        parts = parts[:-1]
    return '/'.join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('dirs', nargs='+')
    parser.add_argument('--json')
    args = parser.parse_args()
    if len(args.dirs) < 2:
        sys.exit('at least two platform directories are required')

    outputs = {d: collect_outputs(d) for d in args.dirs}
    junit = {d: collect_junit(d) for d in args.dirs}

    problems = []
    rows = []

    # 1. the same tests fail everywhere
    all_failed = sorted(set().union(*(set(j) for j in junit.values())))
    for test_id in all_failed:
        on = [d for d in args.dirs if test_id in junit[d]]
        if len(on) != len(args.dirs):
            messages = {d: junit[d][test_id][1][:120] for d in on}
            if all(INPUT_MISMATCH in m for m in messages.values()):
                rows.append({'test': test_id,
                             'status': 'input files differ from the reference '
                                       f'on {len(on)} of {len(args.dirs)} '
                                       'platforms (Python-side difference)'})
                continue
            problems.append(f'{test_id}: fails on {len(on)} of '
                            f'{len(args.dirs)} platforms: {messages}')
            rows.append({'test': test_id, 'status': 'platform-dependent'})

    # 2. outputs of tests that differ from the reference are identical
    tests = sorted(set().union(*(set(f) for f in outputs.values())))
    for test in tests:
        hashes = [outputs[d].get(test) for d in args.dirs]
        if any(h is None for h in hashes):
            problems.append(f'{test}: differs from the reference on only some '
                            'platforms')
            status = 'platform-dependent'
        elif len(set(hashes)) == 1:
            status = 'identical on all platforms'
        else:
            problems.append(f'{test}: outputs differ between platforms')
            status = 'differs between platforms'
        rows.append({'test': test, 'status': status,
                     'sha256': dict(zip(args.dirs, hashes))})

    # 3. failures with nothing to compare
    with_output = {test_dir(t) for t in tests}
    for test_id in all_failed:
        if all(test_id in junit[d] for d in args.dirs) and \
                not any(test_dir(test_id).startswith(w) for w in with_output):
            kind, message = junit[args.dirs[0]][test_id]
            if INPUT_MISMATCH in message:
                status = 'input files differ from the reference on all platforms'
            elif kind == 'failure' and 'AssertionError' in message:
                # A result comparison that keeps no file (for example the
                # depletion tests compare HDF5 references directly)
                status = 'fails the reference comparison on all platforms ' \
                         '(no output file to compare)'
            else:
                status = f'{kind} on all platforms: {message[:120]}'
                problems.append(f'{test_id}: {status}')
            rows.append({'test': test_id, 'status': status})

    print(f'{len(all_failed)} tests failed or errored, {len(tests)} of them '
          'with output files')
    for row in rows:
        print(f"  {row['test'][:70]:70s} {row['status']}")

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump({'platforms': args.dirs, 'tests': rows}, fh, indent=2)

    if problems:
        print('\n'.join(problems))
        return 1
    print(f'all {len(args.dirs)} platforms agree')
    return 0


if __name__ == '__main__':
    sys.exit(main())
