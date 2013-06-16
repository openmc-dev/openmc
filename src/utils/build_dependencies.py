#!/usr/bin/env python2

import glob
import re

dependencies = {}

for src in glob.iglob('*.F90'):
    module = src.strip('.F90')
    deps = set()
    d = re.findall(r'\n\s*use\s+(\w+)',
                   open(src,'r').read())
    for name in d:
        if name in ['mpi','hdf5','h5lt']:
            continue
        if name.startswith('xml_data_'):
            name = name.replace('xml_data_', 'templates/')
        deps.add(name)
    if deps:
        dependencies[module] = sorted(list(deps))


for module in dependencies.keys():
    for dep in dependencies[module]:
        print("{0}.o: {1}.o".format(module, dep))
    print('')
