"""This module contains a few utility functions for reading ENDF_ data. It is by
no means enough to read an entire ENDF file.  For a more complete ENDF reader,
see Pyne_.

.. _ENDF: http://www.nndc.bnl.gov/endf
.. _Pyne: http://www.pyne.io

"""

import re

def read_float(float_string):
    """Parse ENDF 6E11.0 formatted string into a float."""
    assert len(float_string) == 11
    pattern = r'([\s\-]\d+\.\d+)([\+\-]\d+)'
    return float(re.sub(pattern, r'\1e\2', float_string))


def read_CONT_line(line):
    """Parse 80-column line from ENDF CONT record into floats and ints."""
    return (read_float(line[0:11]), read_float(line[11:22]), int(line[22:33]),
            int(line[33:44]), int(line[44:55]), int(line[55:66]),
            int(line[66:70]), int(line[70:72]), int(line[72:75]),
            int(line[75:80]))


def identify_nuclide(fname):
    """Read the header of an ENDF file and extract identifying information."""
    with open(fname, 'r') as fh:
        # Skip the tape id (TPID).
        line = fh.readline()

        # Read the first HEAD and CONT info.
        line = fh.readline()
        ZA, AW, LRP, LFI, NLIB, NMOD, MAT, MF, MT, NS = read_CONT_line(line)
        line = fh.readline()
        ELIS, STA, LIS, LISO, junk, NFOR, MAT, MF, MT, NS = read_CONT_line(line)

    # Return dictionary of the most important identifying information.
    return {'Z': int(ZA) // 1000,
            'A': int(ZA) % 1000,
            'LFI': bool(LFI),
            'LIS': LIS,
            'LISO': LISO}
