#!/usr/bin/env python 

from __future__ import division

import struct
import sys

import silomesh # https://github.com/nhorelik/silomesh/

################################################################################
def parse_options():
  """Process command line arguments"""

  from optparse import OptionParser
  usage = r"""%prog [options] <voxel_file>"""
  p = OptionParser(usage=usage)
  p.add_option('-o', '--output', action='store', dest='output',
             default='plot.silo', help='Path to output SILO file.')
  parsed = p.parse_args()
  if not parsed[1]:
    p.print_help()
    return parsed
  return parsed

################################################################################
def main(file_, o):

  print file_
  fh = open(file_,'rb')
  header = get_header(fh)
  meshparms = header['dimension'] + header['lower_left'] + header['upper_right']
  nx,ny,nz = meshparms[0], meshparms[1], meshparms[2]
  
  silomesh.init_silo(o.output)
  silomesh.init_mesh('plot', *meshparms)
  
  silomesh.init_var("cell_or_mat_id")
  
  for x in range(1,nx+1):
    sys.stdout.write(" {0}%\r".format(int(x/nx*100)))
    sys.stdout.flush()
    for y in range(1,ny+1):
      for z in range(1,nz+1):
        id_ = get_int(fh)[0]
        silomesh.set_value(float(id_), x, y, z)
  print
  silomesh.finalize_var()
  silomesh.finalize_mesh()
  silomesh.finalize_silo()  

################################################################################
def get_header(file_):
  nx,ny,nz = get_int(file_, 3)
  wx,wy,wz = get_double(file_, 3)
  lx,ly,lz = get_double(file_, 3)
  header = {'dimension':[nx,ny,nz], 'width':[wx,wy,wz], 'lower_left':[lx,ly,lz],
            'upper_right': [lx+wx*nx,ly+wy*ny,lz+wz*nz]}
  return header

################################################################################
def get_data(file_, n, typeCode, size):
  return list(struct.unpack('={0}{1}'.format(n,typeCode),
                            file_.read(n*size)))

################################################################################
def get_int(file_, n=1, path=None):
  return get_data(file_, n, 'i', 4)

################################################################################
def get_double(file_, n=1, path=None):
  return get_data(file_, n, 'd', 8)

################################################################################
if __name__ == '__main__':
  (options, args) = parse_options()
  if args:
    main(args[0],options)
