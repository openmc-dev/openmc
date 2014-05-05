#!/usr/bin/env python2

import sys

if len(sys.argv) != 3:
    print("Must supply element and atom/b-cm")
    sys.exit(1)

element = sys.argv[1]
ao = float(sys.argv[2])

for line in open('abundances_modified.txt', 'r'):
    words = line.split()
    if words[1] == element:
        print('<nuclide name="{0}-{1}" ao="{2:10.4e}" />'.format(
              element, words[2], float(words[3])*ao))
