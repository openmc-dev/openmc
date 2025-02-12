#!/usr/bin/env python3

"""This script can be used to create HDF5 nuclear data libraries used by
OpenMC. There are four different ways you can specify ACE libraries that are to
be converted:

1. List each ACE library as a positional argument. This is very useful in
   conjunction with the usual shell utilities (ls, find, etc.).
2. Use the --xsdir option to specify a MCNP xsdir file.
3. Use the --xsdata option to specify a Serpent xsdata file.

The script does not use any extra information from xsdir/xsdata files to
determine whether the nuclide is metastable. Instead, the --metastable argument
can be used to specify whether the ZAID naming convention follows the NNDC data
convention (1000*Z + A + 300 + 100*m), or the MCNP data convention (essentially
the same as NNDC, except that the first metastable state of Am242 is 95242 and
the ground state is 95642).

"""

import argparse
from functools import partial
import os
from pathlib import Path
import warnings

import openmc.data
from openmc.data.ace import TableType


def main(destination, xsdir, xsdata, libraries, metastable, libver):

    if not destination.is_dir():
        destination.mkdir(parents=True, exist_ok=True)

    ace_libraries = []
    if xsdir is not None:
        ace_libraries.extend(openmc.data.ace.get_libraries_from_xsdir(xsdir))
    elif xsdata is not None:
        ace_libraries.extend(openmc.data.ace.get_libraries_from_xsdata(xsdata))
    else:
        ace_libraries = [Path(lib) for lib in libraries]

    converted = {}
    library = openmc.data.DataLibrary()

    for path in ace_libraries:
        # Check that ACE library exists
        if not os.path.exists(path):
            warnings.warn(f"ACE library '{path}' does not exist.")
            continue

        lib = openmc.data.ace.Library(path)
        for table in lib.tables:
            # Check type of the ACE table and determine appropriate class /
            # conversion function
            if table.data_type == TableType.NEUTRON_CONTINUOUS:
                name = table.zaid
                cls = openmc.data.IncidentNeutron
                converter = partial(cls.from_ace, metastable_scheme=metastable)
            elif table.data_type == TableType.THERMAL_SCATTERING:
                # Adjust name to be the new thermal scattering name
                name = openmc.data.get_thermal_name(table.zaid)
                cls = openmc.data.ThermalScattering
                converter = cls.from_ace
            else:
                print(f"Can't convert ACE table {table.name}")
                continue

            if name not in converted:
                try:
                    data = converter(table)
                except Exception as e:
                    print(f"Failed to convert {table.name}: {e}")
                    continue

                print(f"Converting {table.name} (ACE) to {data.name} (HDF5)")

                # Determine output filename
                outfile = destination / (data.name.replace(".", "_") + ".h5")
                data.export_to_hdf5(outfile, "w", libver=libver)

                # Register with library
                library.register_file(outfile)

                # Add nuclide to list
                converted[name] = outfile
            else:
                # Read existing HDF5 file
                data = cls.from_hdf5(converted[name])

                # Add data for new temperature
                try:
                    print(f"Converting {table.name} (ACE) to {data.name} (HDF5)")
                    if table.data_type == TableType.NEUTRON_CONTINUOUS:
                        data.add_temperature_from_ace(table, metastable)
                    else:
                        data.add_temperature_from_ace(table)
                except Exception as e:
                    print(f"Failed to convert {table.name}: {e}")
                    continue

                # Re-export
                data.export_to_hdf5(converted[name], "w", libver=libver)

    # Write cross_sections.xml
    library.export_to_xml(destination / "cross_sections.xml")


if __name__ == "__main__":

    class CustomFormatter(
        argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
    ):
        pass

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=CustomFormatter
    )
    parser.add_argument("libraries", nargs="*", help="ACE libraries to convert to HDF5")
    parser.add_argument(
        "-d",
        "--destination",
        type=Path,
        default=Path.cwd(),
        help="Directory to create new library in",
    )
    parser.add_argument(
        "-m",
        "--metastable",
        choices=["mcnp", "nndc"],
        default="nndc",
        help="How to interpret ZAIDs for metastable nuclides",
    )
    parser.add_argument("--xsdir", help="MCNP xsdir file that lists " "ACE libraries")
    parser.add_argument(
        "--xsdata", help="Serpent xsdata file that lists " "ACE libraries"
    )
    parser.add_argument(
        "--libver",
        choices=["earliest", "latest"],
        default="earliest",
        help="Output HDF5 versioning. Use "
        "'earliest' for backwards compatibility or 'latest' "
        "for performance",
    )
    args = parser.parse_args()

    main(
        destination=args.destination,
        xsdir=args.xsdir,
        xsdata=args.xsdata,
        libraries=args.libraries,
        metastable=args.metastable,
        libver=args.libver,
    )
