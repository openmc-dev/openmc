#!/usr/bin/env python2
"""Convert binary particle track to VTK poly data.

Usage information can be obtained by running 'track.py --help':

    usage: track.py [-h] [-o OUT] IN [IN ...]

    Convert particle track file to a .pvtp file.

    positional arguments:
      IN                    Input particle track data filename(s).

    optional arguments:
      -h, --help            show this help message and exit
      -o OUT, -out OUT, --out OUT
                            Output VTK poly data filename.

"""

import os
import argparse
import struct
import vtk


def _parse_args():
    # Create argument parser.
    parser = argparse.ArgumentParser(
        description='Convert particle track file to a .pvtp file.')
    parser.add_argument('input', metavar='IN', type=str, nargs='+',
                        help='Input particle track data filename(s).')
    parser.add_argument('-o', '-out', '--out', metavar='OUT', type=str,
                        dest='out',
                        help='Output VTK poly data filename.')

    # Parse and return commandline arguments.
    return parser.parse_args()


def main():
    # Parse commandline arguments.
    args = _parse_args()

    # Check input file extensions.
    for fname in args.input:
        if len(fname) > 3 and fname[-3:] == '.h5': continue
        if len(fname) > 7 and fname[-7:] == '.binary': continue
        raise ValueError("Input file names must either end with '.h5' or "
                         "'.binary'.")
    
    # Make sure that the output filename ends with '.pvtp'.
    if not args.out:
        args.out = 'tracks.pvtp'
    elif os.path.splitext(args.out)[1] != '.pvtp':
        args.out = ''.join([args.out, '.pvtp'])

    # Import HDF library if HDF files are present
    for fname in args.input:
        if fname[-3:] != '.h5': continue
        import h5py
        break

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    j = 0
    k = 0
    for fname in args.input:
        if len(fname) > 7 and fname[-7:] == '.binary':
            track = open(fname, 'rb').read()
            coords = [struct.unpack("ddd", track[24*i : 24*(i+1)])
                      for i in range(len(track)/24)]
            n_points = len(coords)
            for triplet in coords: points.InsertNextPoint(triplet)
        elif fname[-3:] == '.h5':
            coords = h5py.File(fname).get('coordinates')
            n_points = coords.shape[0]
            for i in range(n_points): points.InsertNextPoint(coords[i,:])
                
        # Create VTK line and assign points to line.
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(n_points)
        for i in range(n_points):
            line.GetPointIds().SetId(i, j+i)
        
        cells.InsertNextCell(line)
        j += n_points
        k += 1
    global data
    data = vtk.vtkPolyData()
    data.SetPoints(points)
    data.SetLines(cells)
    
    writer = vtk.vtkXMLPPolyDataWriter()
    writer.SetInput(data)
    writer.SetFileName(args.out)
    writer.Write()



if __name__ == '__main__':
    main()


#### The following code is retained for debugging purposes:

##poly_mapper = vtk.vtkPolyDataMapper()
##poly_mapper.SetInputConnection(data.GetProducerPort())
##
##poly_actor = vtk.vtkActor()
##poly_actor.SetMapper(poly_mapper)
##
##ren1 = vtk.vtkRenderer()
##ren1.AddActor(poly_actor)
##ren1.SetBackground(0.0, 0.0, 0.0)
##
##ren1.ResetCamera()
##
##ren_win = vtk.vtkRenderWindow()
##ren_win.AddRenderer(ren1)
##ren_win.SetSize(400, 400)
##
##iren = vtk.vtkRenderWindowInteractor()
##iren.SetRenderWindow(ren_win)
##iren.Initialize()
##iren.Start()
