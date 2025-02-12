#!/usr/bin/env python3

"""Convert HDF5 particle track to VTK poly data.

"""

import argparse

import openmc
import vtk


def _parse_args():
    # Create argument parser.
    parser = argparse.ArgumentParser(
        description='Convert particle track file(s) to a .pvtp file.')
    parser.add_argument('input', metavar='IN', type=str,
                        help='Input particle track data filename.')
    parser.add_argument('-o', '--out', metavar='OUT', type=str, dest='out',
                        help='Output VTK poly data filename.')

    # Parse and return commandline arguments.
    return parser.parse_args()


def main(out=None):
    
    if out is None:
        # Parse commandline arguments.
        args = _parse_args()
        out = args.out
        input = args.input


    # Make sure that the output filename ends with '.pvtp'.
    if not out:
        out = 'tracks.pvtp'
    elif not out.endswith('.pvtp'):
        out += '.pvtp'

    # Write coordinate values to points array.
    track_file = openmc.Tracks(input)
    track_file.write_to_vtk(out)

if __name__ == '__main__':
    main()
