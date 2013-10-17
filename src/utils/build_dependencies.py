#!/usr/bin/env python2

import glob
import re

dependencies = {}

for src in glob.iglob('*.F90'):
    module = src.strip('.F90')
    deps = set()
    d = re.findall(r'\n\s*use\s+(\w+)',
                   open(src, 'r').read())
    for name in d:
        if name in ['mpi', 'hdf5', 'h5lt', 'petscsys', 'petscmat', 'petscksp',
                    'petscsnes', 'petscvec', 'omp_lib', 'fox_dom']:
            continue
        deps.add(name)
    if deps:
        dependencies[module] = sorted(list(deps))


for module in sorted(dependencies.keys()):
    for dep in dependencies[module]:
        print("{0}.o: {1}.o".format(module, dep))
    print('')
