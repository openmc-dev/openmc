#!/usr/bin/env python2

from __future__ import division, print_function

import struct
import sys

################################################################################
def parse_options():
  """Process command line arguments"""

  from optparse import OptionParser
  usage = r"""%prog [options] <voxel_file>"""
  p = OptionParser(usage=usage)
  p.add_option('-o', '--output', action='store', dest='output',
             default='plot', help='Path to output SILO or VTK file.')
  p.add_option('-v', '--vtk', action='store_true', dest='vtk',
           default=False, help='Flag to convert to VTK instead of SILO.')
  parsed = p.parse_args()
  if not parsed[1]:
    p.print_help()
    return parsed
  return parsed

################################################################################
def main(file_, o):

  print(file_)
  fh = open(file_,'rb')
  header = get_header(fh)
  meshparms = header['dimension'] + header['lower_left'] + header['upper_right']
  nx,ny,nz = meshparms[0], meshparms[1], meshparms[2]
  ll = header['lower_left']

  if o.vtk:
    try:
      import vtk
    except:
      print('The vtk python bindings do not appear to be installed properly.\n'
            'On Ubuntu: sudo apt-get install python-vtk\n'
            'See: http://www.vtk.org/')
      return

    origin = [(l+w*n/2.) for n,l,w in zip((nx,ny,nz),ll,header['width'])]

    grid = vtk.vtkImageData()
    grid.SetDimensions(nx+1,ny+1,nz+1)
    grid.SetOrigin(*ll)
    grid.SetSpacing(*header['width'])

    data = vtk.vtkDoubleArray()
    data.SetName("id")
    data.SetNumberOfTuples(nx*ny*nz)
    for x in range(nx):
      sys.stdout.write(" {0}%\r".format(int(x/nx*100)))
      sys.stdout.flush()
      for y in range(ny):
        for z in range(nz):
          i = z*nx*ny + y*nx + x
          id_ = get_int(fh)[0]
          data.SetValue(i, id_)
    grid.GetCellData().AddArray(data)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetInput(grid)
    if not o.output[-4:] == ".vti": o.output += ".vti"
    writer.SetFileName(o.output)
    writer.Write()

  else:

    try:
      import silomesh
    except:
      print('The silomesh package does not appear to be installed properly.\n'
            'See: https://github.com/nhorelik/silomesh/')
      return
    if not o.output[-5:] == ".silo": o.output += ".silo"
    silomesh.init_silo(o.output)
    silomesh.init_mesh('plot', *meshparms)
    silomesh.init_var("id")
    for x in range(1,nx+1):
      sys.stdout.write(" {0}%\r".format(int(x/nx*100)))
      sys.stdout.flush()
      for y in range(1,ny+1):
        for z in range(1,nz+1):
          id_ = get_int(fh)[0]
          silomesh.set_value(float(id_), x, y, z)
    print()
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
