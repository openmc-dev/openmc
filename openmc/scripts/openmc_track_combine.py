#!/usr/bin/env python3

"""Combine multiple HDF5 particle track files."""

import argparse

import openmc


def main(input=None, out=None):
    
    if input is None or out is None:
        # Parse command-line arguments
        parser = argparse.ArgumentParser(
            description='Combine particle track files into a single .h5 file.')
        parser.add_argument('input', metavar='IN', nargs='+',
                            help='Input HDF5 particle track filename(s).')
        parser.add_argument('-o', '--out', metavar='OUT', default='tracks.h5',
                            help='Output HDF5 particle track file.')

        args = parser.parse_args()
        input = args.input
        out = args.out

    openmc.Tracks.combine(input, out)


if __name__ == '__main__':
    main()
