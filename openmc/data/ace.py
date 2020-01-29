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

from collections import OrderedDict
import enum
from pathlib import PurePath, Path
import struct
import sys

import numpy as np

from openmc.mixin import EqualityMixin
import openmc.checkvalue as cv
from .data import ATOMIC_SYMBOL, gnd_name
from .endf import ENDF_FLOAT_RE


def get_metadata(zaid, metastable_scheme='nndc'):
    """Return basic identifying data for a nuclide with a given ZAID.

    Parameters
    ----------
    zaid : int
        ZAID (1000*Z + A) obtained from a library
    metastable_scheme : {'nndc', 'mcnp'}
        Determine how ZAID identifiers are to be interpreted in the case of
        a metastable nuclide. Because the normal ZAID (=1000*Z + A) does not
        encode metastable information, different conventions are used among
        different libraries. In MCNP libraries, the convention is to add 400
        for a metastable nuclide except for Am242m, for which 95242 is
        metastable and 95642 (or 1095242 in newer libraries) is the ground
        state. For NNDC libraries, ZAID is given as 1000*Z + A + 100*m.

    Returns
    -------
    name : str
        Name of the table
    element : str
        The atomic symbol of the isotope in the table; e.g., Zr.
    Z : int
        Number of protons in the nucleus
    mass_number : int
        Number of nucleons in the nucleus
    metastable : int
        Metastable state of the nucleus. A value of zero indicates ground state.

    """

    cv.check_type('zaid', zaid, int)
    cv.check_value('metastable_scheme', metastable_scheme, ['nndc', 'mcnp'])

    Z = zaid // 1000
    mass_number = zaid % 1000

    if metastable_scheme == 'mcnp':
        if zaid > 1000000:
            # New SZA format
            Z = Z % 1000
            if zaid == 1095242:
                metastable = 0
            else:
                metastable = zaid // 1000000
        else:
            if zaid == 95242:
                metastable = 1
            elif zaid == 95642:
                metastable = 0
            else:
                metastable = 1 if mass_number > 300 else 0
    elif metastable_scheme == 'nndc':
        metastable = 1 if mass_number > 300 else 0

    while mass_number > 3 * Z:
        mass_number -= 100

    # Determine name
    element = ATOMIC_SYMBOL[Z]
    name = gnd_name(Z, mass_number, metastable)

    return (name, element, Z, mass_number, metastable)


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
    ascii = open(str(ascii_file), 'r')

    # Set default record length
    record_length = 4096

    # Read data from ASCII file
    lines = ascii.readlines()
    ascii.close()

    # Open binary file
    binary = open(str(binary_file), 'wb')

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

    if name is None:
        return Library(filename).tables[0]
    else:
        lib = Library(filename, name)
        if lib.tables:
            return lib.tables[0]
        else:
            raise ValueError('Could not find ACE table with name: {}'
                             .format(name))


class Library(EqualityMixin):
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
        if isinstance(table_names, str):
            table_names = [table_names]
        if table_names is not None:
            table_names = set(table_names)

        self.tables = []

        # Determine whether file is ASCII or binary
        filename = str(filename)
        try:
            fh = open(filename, 'rb')
            # Grab 10 lines of the library
            sb = b''.join([fh.readline() for i in range(10)])

            # Try to decode it with ascii
            sb.decode('ascii')

            # No exception so proceed with ASCII - reopen in non-binary
            fh.close()
            with open(filename, 'r') as fh:
                self._read_ascii(fh, table_names, verbose)
        except UnicodeDecodeError:
            fh.close()
            with open(filename, 'rb') as fh:
                self._read_binary(fh, table_names, verbose)

    def _read_binary(self, ace_file, table_names, verbose=False,
                     recl_length=4096, entries=512):
        """Read a binary (Type 2) ACE table.

        Parameters
        ----------
        ace_file : file
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
            start_position = ace_file.tell()

            # Check for end-of-file
            if len(ace_file.read(1)) == 0:
                return
            ace_file.seek(start_position)

            # Read name, atomic mass ratio, temperature, date, comment, and
            # material
            name, atomic_weight_ratio, temperature, date, comment, mat = \
                struct.unpack(str('=10sdd10s70s10s'), ace_file.read(116))
            name = name.decode().strip()

            # Read ZAID/awr combinations
            data = struct.unpack(str('=' + 16*'id'), ace_file.read(192))
            pairs = list(zip(data[::2], data[1::2]))

            # Read NXS
            nxs = list(struct.unpack(str('=16i'), ace_file.read(64)))

            # Determine length of XSS and number of records
            length = nxs[0]
            n_records = (length + entries - 1)//entries

            # verify that we are supposed to read this table in
            if (table_names is not None) and (name not in table_names):
                ace_file.seek(start_position + recl_length*(n_records + 1))
                continue

            if verbose:
                kelvin = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, kelvin))

            # Read JXS
            jxs = list(struct.unpack(str('=32i'), ace_file.read(128)))

            # Read XSS
            ace_file.seek(start_position + recl_length)
            xss = list(struct.unpack(str('={0}d'.format(length)),
                                     ace_file.read(length*8)))

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
            ace_file.seek(start_position + recl_length*(n_records + 1))

    def _read_ascii(self, ace_file, table_names, verbose=False):
        """Read an ASCII (Type 1) ACE table.

        Parameters
        ----------
        ace_file : file
            Open ACE file
        table_names : None, str, or iterable
            Tables from the file to read in.  If None, reads in all of the
            tables. If str, reads in only the single table of a matching name.
        verbose : str, optional
            Whether to display what tables are being read. Defaults to False.

        """

        tables_seen = set()

        lines = [ace_file.readline() for i in range(13)]

        while len(lines) != 0 and lines[0].strip() != '':
            # Read name of table, atomic mass ratio, and temperature. If first
            # line is empty, we are at end of file

            # check if it's a 2.0 style header
            if lines[0].split()[0][1] == '.':
                words = lines[0].split()
                name = words[1]
                words = lines[1].split()
                atomic_weight_ratio = float(words[0])
                temperature = float(words[1])
                commentlines = int(words[3])
                for _ in range(commentlines):
                    lines.pop(0)
                    lines.append(ace_file.readline())
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

            # Ensure that we have more tables to read in
            if (table_names is not None) and (table_names <= tables_seen):
                break
            tables_seen.add(name)

            # verify that we are supposed to read this table in
            if (table_names is not None) and (name not in table_names):
                for _ in range(n_lines - 1):
                    ace_file.readline()
                lines = [ace_file.readline() for i in range(13)]
                continue

            # Read lines corresponding to this table
            lines += [ace_file.readline() for i in range(n_lines - 1)]

            if verbose:
                kelvin = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, kelvin))

            # Insert zeros at beginning of NXS, JXS, and XSS arrays so that the
            # indexing will be the same as Fortran. This makes it easier to
            # follow the ACE format specification.
            datastr = '0 ' + ' '.join(lines[8:12])
            jxs = np.fromstring(datastr, dtype=int, sep=' ')

            datastr = '0.0 ' + ''.join(lines[12:12+n_lines])
            xss = np.fromstring(datastr, sep=' ')

            # When NJOY writes an ACE file, any values less than 1e-100 actually
            # get written without the 'e'. Thus, what we do here is check
            # whether the xss array is of the right size (if a number like
            # 1.0-120 is encountered, np.fromstring won't capture any numbers
            # after it). If it's too short, then we apply the ENDF float regular
            # expression. We don't do this by default because it's expensive!
            if xss.size != nxs[1] + 1:
                datastr = ENDF_FLOAT_RE.sub(r'\1e\2\3', datastr)
                xss = np.fromstring(datastr, sep=' ')
                assert xss.size == nxs[1] + 1

            table = Table(name, atomic_weight_ratio, temperature, pairs,
                          nxs, jxs, xss)
            self.tables.append(table)

            # Read all data blocks
            lines = [ace_file.readline() for i in range(13)]


class TableType(enum.Enum):
    """Type of ACE data table."""
    NEUTRON_CONTINUOUS = 'c'
    NEUTRON_DISCRETE = 'd'
    THERMAL_SCATTERING = 't'
    DOSIMETRY = 'y'
    PHOTOATOMIC = 'p'
    PHOTONUCLEAR = 'u'
    PROTON = 'h'
    DEUTERON = 'o'
    TRITON = 'r'
    HELIUM3 = 's'
    ALPHA = 'a'

    @classmethod
    def from_suffix(cls, suffix):
        """Determine ACE table type from a suffix.

        Parameters
        ----------
        suffix : str
            Single letter ACE table designator, e.g., 'c'

        Returns
        -------
        TableType
            ACE table type

        """
        for member in cls:
            if suffix.endswith(member.value):
                return member
        raise ValueError("Suffix '{}' has no corresponding ACE table type."
                         .format(suffix))


class Table(EqualityMixin):
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

    Attributes
    ----------
    data_type : TableType
        Type of the ACE data

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

    @property
    def zaid(self):
        return self.name.split('.')[0]

    @property
    def data_type(self):
        xs = self.name.split('.')[1]
        return TableType.from_suffix(xs[-1])

    def __repr__(self):
        return "<ACE Table: {}>".format(self.name)


def get_libraries_from_xsdir(path):
    """Determine paths to ACE files from an MCNP xsdir file.

    Parameters
    ----------
    path : str or path-like
        Path to xsdir file

    Returns
    -------
    list
        List of paths to ACE libraries
    """
    xsdir = Path(path)

    # Find 'directory' section
    with open(path, 'r') as fh:
        lines = fh.readlines()
    for index, line in enumerate(lines):
        if line.strip().lower() == 'directory':
            break
    else:
        raise RuntimeError("Could not find 'directory' section in MCNP xsdir file")

    # Handle continuation lines indicated by '+' at end of line
    lines = lines[index + 1:]
    continue_lines = [i for i, line in enumerate(lines)
                      if line.strip().endswith('+')]
    for i in reversed(continue_lines):
        lines[i] += lines[i].strip()[:-1] + lines.pop(i + 1)

    # Create list of ACE libraries -- we use an ordered dictionary while
    # building to get O(1) membership checks while retaining insertion order
    libraries = OrderedDict()
    for line in lines:
        words = line.split()
        if len(words) < 3:
            continue

        lib = (xsdir.parent / words[2]).resolve()
        if lib not in libraries:
            # Value in dictionary is not used, so we just assign None. Below a
            # list is created from the keys alone
            libraries[lib] = None

    return list(libraries.keys())


def get_libraries_from_xsdata(path):
    """Determine paths to ACE files from a Serpent xsdata file.

    Parameters
    ----------
    path : str or path-like
        Path to xsdata file

    Returns
    -------
    list
        List of paths to ACE libraries
    """
    xsdata = Path(path)
    with open(xsdata, 'r') as xsdata:
        # As in get_libraries_from_xsdir, we use a dict for O(1) membership
        # check while retaining insertion order
        libraries = OrderedDict()
        for line in xsdata:
            words = line.split()
            if len(words) >= 9:
                lib = (xsdata.parent / words[8]).resolve()
                if lib not in libraries:
                    libraries[lib] = None
    return list(libraries.keys())
