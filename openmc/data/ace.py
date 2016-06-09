"""This module is for reading ACE-format cross sections. ACE stands for "A
Compact ENDF" format and originated from work on MCNP_. It is used in a number
of other Monte Carlo particle transport codes.

ACE-format cross sections are typically generated from ENDF_ files through a
cross section processing program like NJOY_. The ENDF data consists of tabulated
thermal data, ENDF/B resonance parameters, distribution parameters in the
unresolved resonance region, and tabulated data in the fast region. After the
ENDF data has been reconstructed and Doppler-broadened, the ACER module
generates ACE-format cross sections.

.. _MCNP: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/
.. _NJOY: http://t2.lanl.gov/codes.shtml
.. _ENDF: http://www.nndc.bnl.gov/endf

"""

from __future__ import division, unicode_literals
import io
from os import SEEK_CUR
import struct
import sys
from warnings import warn

import numpy as np


if sys.version_info[0] >= 3:
    basestring = str


def ascii_to_binary(ascii_file, binary_file):
    """Convert an ACE file in ASCII format (type 1) to binary format (type 2).

    Parameters
    ----------
    ascii_file : str
        Filename of ASCII ACE file
    binary_file : str
        Filename of binary ACE file to be written

    """

    # Open ASCII file
    ascii = open(ascii_file, 'r')

    # Set default record length
    record_length = 4096

    # Read data from ASCII file
    lines = ascii.readlines()
    ascii.close()

    # Open binary file
    binary = open(binary_file, 'wb')

    idx = 0

    while idx < len(lines):
        # check if it's a > 2.0.0 version header
        if lines[idx].split()[0][1] == '.':
            if lines[idx + 1].split()[3] == '3':
                idx = idx + 3
            else:
                raise NotImplementedError('Only backwards compatible ACE'
                                          'headers currently supported')
        # Read/write header block
        hz = lines[idx][:10].encode('UTF-8')
        aw0 = float(lines[idx][10:22])
        tz = float(lines[idx][22:34])
        hd = lines[idx][35:45].encode('UTF-8')
        hk = lines[idx + 1][:70].encode('UTF-8')
        hm = lines[idx + 1][70:80].encode('UTF-8')
        binary.write(struct.pack(str('=10sdd10s70s10s'), hz, aw0, tz, hd, hk, hm))

        # Read/write IZ/AW pairs
        data = ' '.join(lines[idx + 2:idx + 6]).split()
        iz = list(map(int, data[::2]))
        aw = list(map(float, data[1::2]))
        izaw = [item for sublist in zip(iz, aw) for item in sublist]
        binary.write(struct.pack(str('=' + 16*'id'), *izaw))

        # Read/write NXS and JXS arrays. Null bytes are added at the end so
        # that XSS will start at the second record
        nxs = list(map(int, ' '.join(lines[idx + 6:idx + 8]).split()))
        jxs = list(map(int, ' '.join(lines[idx + 8:idx + 12]).split()))
        binary.write(struct.pack(str('=16i32i{0}x'.format(record_length - 500)),
                                 *(nxs + jxs)))

        # Read/write XSS array. Null bytes are added to form a complete record
        # at the end of the file
        n_lines = (nxs[0] + 3)//4
        xss = list(map(float, ' '.join(lines[
            idx + 12:idx + 12 + n_lines]).split()))
        extra_bytes = record_length - ((len(xss)*8 - 1) % record_length + 1)
        binary.write(struct.pack(str('={0}d{1}x'.format(nxs[0], extra_bytes)),
                                 *xss))

        # Advance to next table in file
        idx += 12 + n_lines

    # Close binary file
    binary.close()


def get_table(filename, name=None):
    """Read a single table from an ACE file

    Parameters
    ----------
    filename : str
        Path of the ACE library to load table from
    name : str, optional
        Name of table to load, e.g. '92235.71c'

    Returns
    -------
    openmc.data.ace.Table
        ACE table with specified name. If no name is specified, the first table
        in the file is returned.

    """

    lib = Library(filename)
    if name is None:
        return lib.tables[0]
    else:
        for table in lib.tables:
            if table.name == name:
                return table
        else:
            raise ValueError('Could not find ACE table with name: {}'
                             .format(name))


class Library(object):
    """A Library objects represents an ACE-formatted file which may contain
    multiple tables with data.

    Parameters
    ----------
    filename : str
        Path of the ACE library file to load.
    table_names : None, str, or iterable, optional
        Tables from the file to read in.  If None, reads in all of the
        tables. If str, reads in only the single table of a matching name.
    verbose : bool, optional
        Determines whether output is printed to the stdout when reading a
        Library

    Attributes
    ----------
    tables : list
        List of :class:`Table` instances

    """

    def __init__(self, filename, table_names=None, verbose=False):
        if isinstance(table_names, basestring):
            table_names = [table_names]
        if table_names is not None:
            table_names = set(table_names)

        self.tables = []

        # Determine whether file is ASCII or binary
        try:
            fh = io.open(filename, 'rb')
            # Grab 10 lines of the library
            sb = b''.join([fh.readline() for i in range(10)])

            # Try to decode it with ascii
            sd = sb.decode('ascii')

            # No exception so proceed with ASCII - reopen in non-binary
            fh.close()
            fh = io.open(filename, 'r')
            fh.seek(0)
            self._read_ascii(fh, table_names, verbose)
        except UnicodeDecodeError:
            fh.close()
            fh = open(filename, 'rb')
            self._read_binary(fh, table_names, verbose)

    def _read_binary(self, fh, table_names, verbose=False,
                     recl_length=4096, entries=512):
        """Read a binary (Type 2) ACE table.

        Parameters
        ----------
        fh : file
            Open ACE file
        table_names : None, str, or iterable
            Tables from the file to read in.  If None, reads in all of the
            tables. If str, reads in only the single table of a matching name.
        verbose : str, optional
            Whether to display what tables are being read. Defaults to False.
        recl_length : int, optional
            Fortran record length in binary file. Default value is 4096 bytes.
        entries : int, optional
            Number of entries per record. The default is 512 corresponding to a
            record length of 4096 bytes with double precision data.

        """

        while True:
            start_position = fh.tell()

            # Check for end-of-file
            if len(fh.read(1)) == 0:
                return
            fh.seek(start_position)

            # Read name, atomic mass ratio, temperature, date, comment, and
            # material
            name, atomic_weight_ratio, temperature, date, comment, mat = \
                struct.unpack(str('=10sdd10s70s10s'), fh.read(116))
            name = name.decode().strip()

            # Read ZAID/awr combinations
            data = struct.unpack(str('=' + 16*'id'), fh.read(192))
            pairs = list(zip(data[::2], data[1::2]))

            # Read NXS
            nxs = list(struct.unpack(str('=16i'), fh.read(64)))

            # Determine length of XSS and number of records
            length = nxs[0]
            n_records = (length + entries - 1)//entries

            # verify that we are supposed to read this table in
            if (table_names is not None) and (name not in table_names):
                fh.seek(start_position + recl_length*(n_records + 1))
                continue

            if verbose:
                temperature_in_K = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temperature_in_K))

            # Read JXS
            jxs = list(struct.unpack(str('=32i'), fh.read(128)))

            # Read XSS
            fh.seek(start_position + recl_length)
            xss = list(struct.unpack(str('={0}d'.format(length)),
                                     fh.read(length*8)))

            # Insert zeros at beginning of NXS, JXS, and XSS arrays so that the
            # indexing will be the same as Fortran. This makes it easier to
            # follow the ACE format specification.
            nxs.insert(0, 0)
            nxs = np.array(nxs, dtype=int)

            jxs.insert(0, 0)
            jxs = np.array(jxs, dtype=int)

            xss.insert(0, 0.0)
            xss = np.array(xss)

            # Create ACE table with data read in
            table = Table(name, atomic_weight_ratio, temperature, pairs,
                          nxs, jxs, xss)
            self.tables.append(table)

            # Advance to next record
            fh.seek(start_position + recl_length*(n_records + 1))

    def _read_ascii(self, fh, table_names, verbose=False):
        """Read an ASCII (Type 1) ACE table.

        Parameters
        ----------
        fh : file
            Open ACE file
        table_names : None, str, or iterable
            Tables from the file to read in.  If None, reads in all of the
            tables. If str, reads in only the single table of a matching name.
        verbose : str, optional
            Whether to display what tables are being read. Defaults to False.

        """

        tables_seen = set()

        lines = [fh.readline() for i in range(13)]

        while (0 != len(lines)) and (lines[0] != ''):
            # Read name of table, atomic mass ratio, and temperature. If first
            # line is empty, we are at end of file

            # check if it's a 2.0 style header
            if lines[0].split()[0][1] == '.':
                words = lines[0].split()
                version = words[0]
                name = words[1]
                if len(words) == 3:
                    source = words[2]
                words = lines[1].split()
                atomic_weight_ratio = float(words[0])
                temperature = float(words[1])
                commentlines = int(words[3])
                for i in range(commentlines):
                    lines.pop(0)
                    lines.append(fh.readline())
            else:
                words = lines[0].split()
                name = words[0]
                atomic_weight_ratio = float(words[1])
                temperature = float(words[2])

            datastr = ' '.join(lines[2:6]).split()
            pairs = list(zip(map(int, datastr[::2]),
                             map(float, datastr[1::2])))

            datastr = '0 ' + ' '.join(lines[6:8])
            nxs = np.fromstring(datastr, sep=' ', dtype=int)

            n_lines = (nxs[1] + 3)//4
            n_bytes = len(lines[-1]) * (n_lines - 2) + 1

            # Ensure that we have more tables to read in
            if (table_names is not None) and (table_names < tables_seen):
                break
            tables_seen.add(name)

            # verify that we are suppossed to read this table in
            if (table_names is not None) and (name not in table_names):
                fh.seek(n_bytes, SEEK_CUR)
                fh.readline()
                lines = [fh.readline() for i in range(13)]
                continue

            # read and fix over-shoot
            lines += fh.readlines(n_bytes)
            if 12 + n_lines < len(lines):
                goback = sum([len(line) for line in lines[12+n_lines:]])
                lines = lines[:12+n_lines]
                fh.seek(-goback, SEEK_CUR)

            if verbose:
                temperature_in_K = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temperature_in_K))

            # Read comment
            comment = lines[1].strip()

            # Insert zeros at beginning of NXS, JXS, and XSS arrays so that the
            # indexing will be the same as Fortran. This makes it easier to
            # follow the ACE format specification.
            datastr = '0 ' + ' '.join(lines[8:12])
            jxs = np.fromstring(datastr, dtype=int, sep=' ')

            datastr = '0.0 ' + ''.join(lines[12:12+n_lines])
            xss = np.fromstring(datastr, sep=' ')

            table = Table(name, atomic_weight_ratio, temperature, pairs,
                          nxs, jxs, xss)
            self.tables.append(table)

            # Read all data blocks
            lines = [fh.readline() for i in range(13)]


class Table(object):
    """ACE cross section table

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. '92235.70c'.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in MeV.
    pairs : list of tuple
        16 pairs of ZAIDs and atomic weight ratios. Used for thermal scattering
        tables to indicate what isotopes scattering is applied to.
    nxs : numpy.ndarray
        Array that defines various lengths with in the table
    jxs : numpy.ndarray
        Array that gives locations in the ``xss`` array for various blocks of
        data
    xss : numpy.ndarray
        Raw data for the ACE table

    """
    def __init__(self, name, atomic_weight_ratio, temperature, pairs,
                 nxs, jxs, xss):
        self.name = name
        self.atomic_weight_ratio = atomic_weight_ratio
        self.temperature = temperature
        self.pairs = pairs
        self.nxs = nxs
        self.jxs = jxs
        self.xss = xss

    def __repr__(self):
        return "<ACE Table: {}>".format(self.name)
