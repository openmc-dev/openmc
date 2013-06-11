#!/usr/bin/env python2
"""Convert binary particle track to VTK poly data.

Run 'track.py -help' for usage.

"""

import os
import argparse
import struct
import vtk

if __name__ == '__main__':
    # Create argument parser.
    parser = argparse.ArgumentParser(
        description='Convert binary particle track file to a .pvtp file.')
    parser.add_argument('input', metavar='IN', type=argparse.FileType('rb'),
                        help='Input particle track data filename.')
    parser.add_argument('-o', '-out', '--out', metavar='OUT', type=str, dest='out',
                        help='Output VTK poly data filename.')
    # Parse commandline arguments.
    args = parser.parse_args()

    # Make sure that the output filename ends with '.pvtp'.
    if not args.out:
        args.out = ''.join([os.path.splitext(args.input.name)[0], '.pvtp'])
    elif os.path.splitext(args.out)[1] != '.pvtp':
        args.out = ''.join([args.out, '.pvtp'])

    # Convert binary data into a list of coordinate triplets.
    track = args.input.read()
    coords = [struct.unpack("ddd", track[24*i : 24*(i+1)])
              for i in range(len(track)/24)]

    # Create VTK points.
    points = vtk.vtkPoints()
    for triplet in coords:
        points.InsertNextPoint(triplet)

    # Create VTK line and assign points to line.
    line = vtk.vtkPolyLine()
    line.GetPointIds().SetNumberOfIds(points.GetNumberOfPoints())
    for i in range(points.GetNumberOfPoints()):
        line.GetPointIds().SetId(i,i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(line)

    data = vtk.vtkPolyData()
    data.SetPoints(points)
    data.SetLines(cells)

    writer = vtk.vtkXMLPPolyDataWriter()
    writer.SetInput(data)
    writer.SetFileName(args.out)
    writer.Write()

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
