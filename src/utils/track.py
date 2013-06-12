#!/usr/bin/env python2
"""Convert binary particle track to VTK poly data.

Usage information can be obtained by running 'track.py --help':

    usage: track.py [-h] [-o OUT] IN [IN ...]

    Convert binary particle track file to a .pvtp file.

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
        description='Convert binary particle track file to a .pvtp file.')
    parser.add_argument('input', metavar='IN', type=argparse.FileType('rb'),
                        nargs='+',
                        help='Input particle track data filename(s).')
    parser.add_argument('-o', '-out', '--out', metavar='OUT', type=str,
                        dest='out',
                        help='Output VTK poly data filename.')

    # Parse and return commandline arguments.
    return parser.parse_args()


def main():
    # Parse commandline argumetns.
    args = _parse_args()
    
    # Make sure that the output filename ends with '.pvtp'.
    if not args.out:
        args.out = 'tracks.pvtp'
    elif os.path.splitext(args.out)[1] != '.pvtp':
        args.out = ''.join([args.out, '.pvtp'])

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    j = 0
    k = 0
    for fin in args.input:
        track = fin.read()
        coords = [struct.unpack("ddd", track[24*i : 24*(i+1)])
                  for i in range(len(track)/24)]
        n_points = len(coords)
        
        for triplet in coords:
            points.InsertNextPoint(triplet)
        
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
