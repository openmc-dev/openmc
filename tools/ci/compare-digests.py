#!/usr/bin/env python3
"""Compare track and tally digests produced by tools/ci/track-digest.py on
different platforms and fail when any of them differ.

Usage: compare-digests.py DIGEST.json [DIGEST.json ...]
"""

import json
import sys


def main(paths):
    digests = []
    for path in paths:
        with open(path) as fh:
            digests.append((path, json.load(fh)))

    if len(digests) < 2:
        sys.exit('at least two digest files are required')

    print(f"{'file':40s} {'platform':40s} {'states':>8s} tracks sha256 / "
          "statepoint sha256")
    for path, d in digests:
        print(f"{path:40s} {d['platform'][:40]:40s} "
              f"{d['tracks']['states']:8d} {d['tracks']['sha256'][:16]} / "
              f"{d['statepoint']['sha256'][:16]}")

    reference = digests[0][1]
    mismatches = []
    for path, d in digests[1:]:
        for key in ('tracks', 'statepoint'):
            if d[key]['sha256'] != reference[key]['sha256']:
                mismatches.append(f'{key} digest of {path} differs from '
                                  f'{digests[0][0]}')
        if d['tracks']['states'] != reference['tracks']['states']:
            mismatches.append(f'state count of {path} differs from '
                              f'{digests[0][0]}')

    if mismatches:
        print('\n'.join(mismatches))
        return 1
    print(f'all {len(digests)} digests are identical')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
