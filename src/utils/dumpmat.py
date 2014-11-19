#!/usr/bin/env python2

from __future__ import division, print_function

import struct
import sys

################################################################################
def parse_options():
    """Process command line arguments"""

    from optparse import OptionParser
    usage = r"""%prog [options] <material file>"""
    p = OptionParser(usage=usage)
    parsed = p.parse_args()
    if not parsed[1]:
        p.print_help()
        return parsed
    return parsed

################################################################################
def main(file_, o):

    print(file_)
    fh, hdf5 = open_file(file_)
    n_nuc, n_inst = get_header(fh, hdf5)
    print('n_nuclides: %s' % n_nuc)
    print('n_instances: %s' % n_inst)
    comps = get_comps(file_)
    for c in comps:
        #print(c)
        print(c[0]/sum(c),c[1]/sum(c))

################################################################################
def open_file(filename):

    if filename.endswith('.h5'):
        import h5py
        return h5py.File(filename, 'r'), True
    else:
        return open(filename, 'rb'), False

################################################################################
def get_comps(file_):
    
    fh, hdf5 = open_file(file_)
    n_nuc, n_inst = get_header(fh, hdf5)
    if not hdf5:
        # the rest of the first record is padding
        dummy = get_long(fh, n_nuc-2)
    comps = []
    if hdf5:
        allcomps = get_double(fh, path='comps')
        comps = chunks(allcomps, n_nuc)
    else:
        for i in range(n_inst):
            vals = get_double(fh, n_nuc)
            comps.append(vals)
    fh.close()
    return comps

################################################################################
def get_n_nucs(file_):
    fh, hdf5 = open_file(file_)
    n_nuc, n_inst = get_header(fh, hdf5)
    fh.close()
    return n_nuc

################################################################################
def get_n_inst(file_):
    fh, hdf5 = open_file(file_)
    n_nuc, n_inst = get_header(fh, hdf5)
    fh.close()
    return n_inst

################################################################################
def get_header(file_, hdf5):
    if hdf5:
        header = (get_long(file_, path='n_nuclides')[0],
                  get_long(file_, path='n_instances')[0])
    else:
        header = get_long(file_, 2)
    return header

################################################################################
def get_data(file_, n, typeCode, size):
    s = file_.read(n*size)
    return list(struct.unpack('={0}{1}'.format(n,typeCode), s))

################################################################################
def get_long(file_, n=1, path=None):
    if path is not None:
        return file_[path].value
    else:
        return get_data(file_, n, 'q', 8)

################################################################################
def get_double(file_, n=1, path=None):
    if path is not None:
        return file_[path].value
    else:
        return get_data(file_, n, 'd', 8)

################################################################################
def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

################################################################################
if __name__ == '__main__':
    (options, args) = parse_options()
    if args:
        main(args[0],options)
