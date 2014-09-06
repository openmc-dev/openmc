#!/usr/bin/env python
"""Update lattices in geometry .xml files to the post-hex_lattice format.

Usage information can be obtained by running 'track.py --help':

    usage: update_lattices.py [-h] IN [IN ...]

    Update lattices in geometry.xml files to the latest format.

    positional arguments:
      IN          Input geometry.xml file(s).

    optional arguments:
      -h, --help  show this help message and exit

"""

from __future__ import print_function

import argparse


def parse_args():
    """Read the input files from the commandline."""
    # Create argument parser.
    parser = argparse.ArgumentParser(description=\
                  'Update lattices in geometry.xml files to the latest format.')
    parser.add_argument('input', metavar='IN', type=str, nargs='+',
                        help='Input geometry.xml file(s).')

    # Parse and return commandline arguments.
    return parser.parse_args()


def read_file(fname):
    """Open a file and return a list of lines."""
    fin = open(fname)
    lines = fin.readlines()
    fin.close()
    return lines


def process_lattice(lines):
    """Delete 'type' and replace 'width' with 'pitch' in lattice string."""
    assert 'type' in lines
    assert 'width' in lines

    # Delete the 'type' attribute or tag.
    if '<type>' in lines:
        # 'type' is a subelement.
        start = lines.index('<type>')
        end = lines.index('</type>') + 6
        lines = ''.join([lines[:start], lines[end+1:]])
    else:
        # 'type' is an attribute.
        start = lines.index('type')
        end = lines.index('rectangular') + 11  # adjusted to include end quote
        lines = ''.join([lines[:start], lines[end+1:]])

    # Change 'width' to 'pitch'.
    lines = lines.replace('width', 'pitch')

    return lines


def process_xml(file_lines):
    """Update the lattices in the given geometry.xml lines."""
    # Join the list of lines into one big string to simplify the indexing of
    # the xml tabs.
    lines = ''.join(file_lines)
    # Seperate the lattice elements.
    lats = lines.split('<lattice')
    # Add the lattice tags back to the lattices.
    for i in range(1, len(lats)):
        lats[i] = ''.join(['<lattice', lats[i]])

    # Update the lattices.
    for i in range(1, len(lats)):
         lats[i] = process_lattice(lats[i])

    # Put the xml back together into one big string.
    lines = ''.join(lats)
    # Seperate the xml back into the original lines.
    lines = lines.split('\n')
    # Add the newline characters back on to the end of each line.
    for i in range(len(lines) - 1):
        lines[i] = ''.join([lines[i], '\n'])

    return lines


if __name__ == '__main__':
    args = parse_args()
    for fname in args.input:
        input_lines = read_file(fname)

        # Process files with lattices.
        if any(['<lattice' in x for x in input_lines]):
            output_lines = process_xml(input_lines)
        # Leave files without lattices unchanged.
        else:
            output_lines = input_lines

        fout = open(fname, mode='w')
        fout.write(''.join(output_lines))
        fout.close()

