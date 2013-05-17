#!/usr/bin/env python

import glob

dependencies = {}

for src in glob.iglob('*.F90'):
    module = src.strip('.F90')
    d = set()
    for line in open(src, 'r'):
        words = line.split()
        if words and words[0].lower() == 'use':
            name = words[1].strip(',')
            if name in ['mpi','hdf5','h5lt']:
                continue
            if name.startswith('xml_data_'):
                name = name.replace('xml_data_', 'templates/')
            d.add(name)
    if d:
        d = list(d)
        d.sort()
        dependencies[module] = d


keys = dependencies.keys()
keys.sort()
for module in keys:
    for dep in dependencies[module]:
        print("{0}.o: {1}.o".format(module, dep))
    print('')
