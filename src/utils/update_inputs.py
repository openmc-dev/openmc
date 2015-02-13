#!/usr/bin/env python
"""Update OpenMC's input XML files to the latest format.

Usage information can be obtained by running 'update_inputs.py --help':

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


description = "Update OpenMC's input XML files to the latest format."
epilog = """\
If any of the given files do not match the most up-to-date formatting, then they
will be automatically rewritten.  The old out-of-date files will not be deleted;
they will be moved to a new file with '.original' appended to their name.

Formatting changes that will be made:

geometry.xml:  Lattices containing 'outside' attributes/tags will be replaced
  with lattices containing 'outer' attributes, and the appropriate
  cells/universes will be added.
"""


def parse_args():
    """Read the input files from the commandline."""
    # Create argument parser.
    parser = argparse.ArgumentParser(
         description=description,
         epilog=epilog,
         formatter_class=argparse.RawTextHelpFormatter)
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
    if 'id' in lattice_element.attrib:
        return int(lattice_element.attrib['id'].strip())
    elif any([child.tag == 'id' for child in lattice_element]):
        elem = lattice_element.find('id')
        return int(elem.text.strip())
    else:
        raise RuntimeError('Could not find the id for a lattice.')


def pop_lat_outside(lattice_element):
    """Return lattice's outside material and remove from attributes/elements."""
    assert isinstance(lattice_element, ET.Element)

    # Check attributes.
    if 'outside' in lattice_element.attrib:
        material = lattice_element.attrib['outside'].strip()
        del lattice_element.attrib['outside']

    # Check subelements.
    elif any([child.tag == 'outside' for child in lattice_element]):
        elem = lattice_element.find('outside')
        material = elem.text.strip()
        lattice_element.remove(elem)

    # No 'outside' specified.  This means the outside is a void.
    else:
        material = 'void'

    return material


def update_geometry(geometry_root):
    """Update the given XML geometry tree. Return True if changes were made."""
    root = geometry_root
    was_updated = False

    # Ignore files that do not contain lattices.
    if all([child.tag != 'lattice' for child in root]): return False

    # Get a set of already-used universe and cell ids.
    uids = get_universe_ids(root)
    cids = get_cell_ids(root)
    taken_ids = uids.union(cids)

    # Update the definitions of each lattice
    for lat in root.iter('lattice'):
        # Get the lattice's id.
        lat_id = get_lat_id(lat)

        # Ignore lattices that have 'outer' specified.
        if any([child.tag == 'outer' for child in lat]): continue

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

        was_updated = True

    return was_updated


if __name__ == '__main__':
    args = parse_args()
    for fname in args.input:
        # Parse the XML data.
        tree = ET.parse(fname)
        root = tree.getroot()
        was_updated = False

        if root.tag == 'geometry':
            was_updated = update_geometry(root)

        if was_updated:
            # Move the original geometry file to preserve it.
            move(fname, fname + '.original')

            # Write a new geometry file.
            tree.write(fname)
