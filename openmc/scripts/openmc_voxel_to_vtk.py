#!/usr/bin/env python3

from argparse import ArgumentParser
import openmc

def main(voxel_file=None, output=None):
    if voxel_file is None or output is None:
        # Process command line arguments
        parser = ArgumentParser('Converts a voxel HDF5 file to a VTK file')
        parser.add_argument("voxel_file", help="Path to voxel h5 file")
        parser.add_argument(
            "-o",
            "--output",
            action="store",
            default="plot",
            help="Path to output VTK file.",
        )
        args = parser.parse_args()
        voxel_file = args.voxel_file
        output = args.output

    openmc.voxel_to_vtk(voxel_file, output)
    print(f"Written VTK file {output}...")

if __name__ == "__main__":
    main()