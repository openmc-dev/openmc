#!/usr/bin/env python
"""Update lattices in geometry .xml files to the latest format.

Usage information can be obtained by running 'update_lattices.py --help':

usage: update_lattices.py [-h] IN [IN ...]

Update lattices in geometry.xml files to the latest format. This will remove
'outside' attributes/elements and replace them with 'outer' attributes. Note
that this script will not delete the given files; it will append '.original'
to the given files and write new ones.

positional arguments:
  IN          Input geometry.xml file(s).

optional arguments:
  -h, --help  show this help message and exit

"""

from __future__ import print_function

import argparse
from random import randint
from shutil import move
import xml.etree.ElementTree as ET


def parse_args():
    """Read the input files from the commandline."""
    # Create argument parser.
    parser = argparse.ArgumentParser(description="Update lattices in "
         "geometry.xml files to the latest format.  This will remove 'outside' "
         "attributes/elements and replace them with 'outer' attributes.  Note "
         "that this script will not delete the given files; it will append "
         "'.original' to the given files and write new ones.")
    parser.add_argument('input', metavar='IN', type=str, nargs='+',
                        help='Input geometry.xml file(s).')

    # Parse and return commandline arguments.
    return parser.parse_args()


def get_universe_ids(geometry_root):
    """Return a set of universe id numbers."""
    root = geometry_root
    out = {0}

    # Get the ids of universes defined by cells.
    for cell in root.iter('cell'):
        # Get universe attributes.
        if 'universe' in cell.attrib:
            uid = cell.attrib['universe']
            out.add(int(uid))

        # Get universe elements.
        elif cell.find('universe') is not None:
            elem = cell.find('universe')
            uid = elem.text
            out.add(int(uid))

    # Get the ids of universes defined by lattices.
    for lat in root.iter('lattice'):
        # Get id attributes.
        if 'id' in lat.attrib:
            uid = lat.attrib['id']
            out.add(int(uid))

        # Get id elements.
        elif lat.find('id') is not None:
            elem = lat.find('id')
            uid = elem.text
            out.add(int(uid))

    return out


def get_cell_ids(geometry_root):
    """Return a set of cell id numbers."""
    root = geometry_root
    out = set()

    # Get the ids of universes defined by cells.
    for cell in root.iter('cell'):
        # Get id attributes.
        if 'id' in cell.attrib:
            cid = cell.attrib['id']
            out.add(int(cid))

        # Get id elements.
        elif cell.find('id') is not None:
            elem = cell.find('id')
            cid = elem.text
            out.add(int(cid))

    return out


def find_new_id(current_ids, preferred=None):
    """Return a new id that is not already present in current_ids."""
    distance_from_preferred = 21
    max_random_attempts = 10000

    # First, try to find an id near the preferred number.
    if preferred is not None:
        assert isinstance(preferred, int)
        for i in range(1, distance_from_preferred):
            if (preferred - i not in current_ids) and (preferred - i > 0):
                return preferred - i
            if (preferred + i not in current_ids) and (preferred + i > 0):
                return preferred + i

    # If that was unsuccessful, attempt to randomly guess a new id number.
    for i in range(max_random_attempts):
        num = randint(1, 2147483647)
        if num not in current_inds:
            return num

    # Raise an error if an id was not found.
    raise RuntimeError('Could not find a unique id number for a new universe.')


def get_lat_id(lattice_element):
    """Return the id integer of the lattice_element."""
    assert isinstance(lattice_element, ET.Element)
    if 'id' in lat.attrib:
        return int(lat.attrib['id'].strip())
    elif any([child.tag == 'id' for child in lat]):
        elem = lat.find('id')
        return int(elem.text.strip())
    else:
        raise RuntimeError('Could not find the id for a lattice.')


def pop_lat_outside(lattice_element):
    """Return lattice's outside material and remove from attributes/elements."""
    assert isinstance(lattice_element, ET.Element)

    # Check attributes.
    if 'outside' in lat.attrib:
        material = lat.attrib['outside'].strip()
        del lat.attrib['outside']

    # Check subelements.
    elif any([child.tag == 'outside' for child in lat]):
        elem = lat.find('outside')
        material = elem.text.strip()
        lat.remove(elem)

    # No 'outside' specified.  This means the outside is a void.
    else:
        material = 'void'

    return material



if __name__ == '__main__':
    args = parse_args()
    for fname in args.input:
        # Parse the XML data.
        tree = ET.parse(fname)
        root = tree.getroot()

        # Ignore files that do not contain lattices.
        if all([child.tag != 'lattice' for child in root]): continue

        # Get a set of already-used universe and cell ids.
        uids = get_universe_ids(root)
        cids = get_cell_ids(root)
        taken_ids = uids.union(cids)

        # Update the definitions of each lattice
        for lat in root.iter('lattice'):
            # Get the lattice's id.
            lat_id = get_lat_id(lat)

            # Pop the 'outside' material.
            material = pop_lat_outside(lat)

            # Get an id number for a new outer universe.  Ideally, the id should
            # be close to the lattice's id.
            new_uid = find_new_id(taken_ids, preferred=lat_id)
            assert new_uid not in taken_ids

            # Add the new universe filled with the old 'outside' material to the
            # geometry.
            new_cell = ET.Element('cell')
            new_cell.attrib['id'] = str(new_uid)
            new_cell.attrib['universe'] = str(new_uid)
            new_cell.attrib['material'] = material
            root.append(new_cell)
            taken_ids.add(new_uid)

            # Add the new universe to the lattice's 'outer' attribute.
            lat.attrib['outer'] = str(new_uid)

        # Move the original geometry file to preserve it.
        move(fname, ''.join((fname, '.original')))

        # Write a new geometry file.
        tree.write(fname)
