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
  fh = open(file_,'rb')
  n_nuc, n_inst = get_header(fh)
  dummy = get_long(fh, n_nuc-2) # the rest of the first record is padding
  print('n_nuclides: %s' % n_nuc)
  print('n_instances: %s' % n_inst)
  for i in range(n_inst):
    vals = get_double(fh, n_nuc)
    #print(vals)
    print(vals[0]/sum(vals),vals[1]/sum(vals))

################################################################################
def get_header(file_):
  header = get_long(file_, 2)
  return header

################################################################################
def get_data(file_, n, typeCode, size):
  s = file_.read(n*size)
  return list(struct.unpack('={0}{1}'.format(n,typeCode), s))

################################################################################
def get_long(file_, n=1, path=None):
  return get_data(file_, n, 'q', 8)

################################################################################
def get_double(file_, n=1, path=None):
  return get_data(file_, n, 'd', 8)

################################################################################
if __name__ == '__main__':
  (options, args) = parse_options()
  if args:
    main(args[0],options)
