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
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from numpy.polynomial import Polynomial
import h5py

from . import atomic_number, atomic_symbol, reaction_name
from .container import Tabulated1D, interpolation_scheme
from .angle_distribution import AngleDistribution
from .energy_distribution import *
from .product import Product
from .angle_energy import AngleEnergy
from .kalbach_mann import KalbachMann
from .uncorrelated import UncorrelatedAngleEnergy
from .correlated import CorrelatedAngleEnergy
from .nbody import NBodyPhaseSpace
from .thermal import CoherentElastic
from .urr import ProbabilityTables
from openmc.stats import Tabular, Discrete, Uniform, Mixture

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


def _get_tabulated_1d(array, idx=0):
    """Create a Tabulated1D object from array.

    Parameters
    ----------
    array : numpy.ndarray
        Array is formed as a 1 dimensional array as follows: [number of regions,
        final pair for each region, interpolation parameters, number of pairs,
        x-values, y-values]
    idx : int, optional
        Offset to read from in array (default of zero)

    Returns
    -------
    openmc.data.Tabulated1D
        Tabulated data object

    """

    # Get number of regions and pairs
    n_regions = int(array[idx])
    n_pairs = int(array[idx + 1 + 2*n_regions])

    # Get interpolation information
    idx += 1
    if n_regions > 0:
        nbt = np.asarray(array[idx:idx + n_regions], dtype=int)
        interp = np.asarray(array[idx + n_regions:idx + 2*n_regions], dtype=int)
    else:
        # NR=0 regions implies linear-linear interpolation by default
        nbt = np.array([n_pairs])
        interp = np.array([2])

    # Get (x,y) pairs
    idx += 2*n_regions + 1
    x = array[idx:idx + n_pairs]
    y = array[idx + n_pairs:idx + 2*n_pairs]

    return Tabulated1D(x, y, nbt, interp)


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
        return list(lib.tables.values())[0]
    else:
        return lib.tables[name]


def get_all_tables(filename):
    """Read all tables from an ACE file

    Parameters
    ----------
    filename : str
        Path of the ACE library to load table from
    name : str, optional
        Name of table to load, e.g. '92235.71c'

    Returns
    -------
    list of openmc.data.ace.Table
        ACE tables read from the file

    """

    lib = Library(filename)
    return list(lib.tables.values())


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
    tables : dict
        Dictionary whose keys are the names of the ACE tables and whose values
        are the instances of subclasses of :class:`Table`
        (e.g. :class:`NeutronTable`)

    """

    def __init__(self, filename, table_names=None, verbose=False):
        if isinstance(table_names, basestring):
            table_names = [table_names]
        if table_names is not None:
            table_names = set(table_names)

        self.tables = {}

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
            name = name.strip()

            # Read ZAID/awr combinations
            izaw_pairs = struct.unpack(str('=' + 16*'id'), fh.read(192))

            # Read NXS
            nxs = list(struct.unpack(str('=16i'), fh.read(64)))

            # Determine length of XSS and number of records
            length = nxs[0]
            n_records = (length + entries - 1)//entries

            # name is bytes, make it a string
            name = name.decode()
            # verify that we are supposed to read this table in
            if (table_names is not None) and (name not in table_names):
                fh.seek(start_position + recl_length*(n_records + 1))
                continue

            # ensure we have a valid table type
            if len(name) == 0 or name[-1] not in table_types:
                warn("Unsupported table: " + name, RuntimeWarning)
                fh.seek(start_position + recl_length*(n_records + 1))
                continue

            # get the table
            table = table_types[name[-1]](name, atomic_weight_ratio, temperature)

            if verbose:
                temperature_in_K = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temperature_in_K))
            self.tables[name] = table

            # If table is S(a,b), add zaids
            zaids = np.array(izaw_pairs[::2])
            table.zaids = zaids[np.nonzero(zaids)]

            # Read JXS
            jxs = list(struct.unpack(str('=32i'), fh.read(128)))

            # Read XSS
            fh.seek(start_position + recl_length)
            xss = list(struct.unpack(str('={0}d'.format(length)),
                                     fh.read(length*8)))

            # Insert empty object at beginning of NXS, JXS, and XSS arrays so
            # that the indexing will be the same as Fortran. This makes it
            # easier to follow the ACE format specification.
            nxs.insert(0, 0)
            table._nxs = np.array(nxs, dtype=int)

            jxs.insert(0, 0)
            table._jxs = np.array(jxs, dtype=int)

            xss.insert(0, 0.0)
            table._xss = np.array(xss)

            # Read all data blocks
            table._read_all()

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

            izaw_pairs = (' '.join(lines[2:6])).split()

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

            # ensure we have a valid table type
            if len(name) == 0 or name[-1] not in table_types:
                warn("Unsupported table: " + name, RuntimeWarning)
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

            # get the table
            table = table_types[name[-1]](name, atomic_weight_ratio, temperature)

            if verbose:
                temperature_in_K = round(temperature * 1e6 / 8.617342e-5)
                print("Loading nuclide {0} at {1} K".format(name, temperature_in_K))
            self.tables[name] = table

            # Read comment
            table.comment = lines[1].strip()

            # If table is S(a,b), add zaids
            if isinstance(table, SabTable):
                zaids = np.fromiter(map(int, izaw_pairs[::2]), int)
                table.zaids = zaids[np.nonzero(zaids)]

            # Add NXS, JXS, and XSS arrays to table Insert empty object at
            # beginning of NXS, JXS, and XSS arrays so that the indexing will be
            # the same as Fortran. This makes it easier to follow the ACE format
            # specification.
            table._nxs = nxs

            datastr = '0 ' + ' '.join(lines[8:12])
            table._jxs = np.fromstring(datastr, dtype=int, sep=' ')

            datastr = '0.0 ' + ''.join(lines[12:12+n_lines])
            table._xss = np.fromstring(datastr, sep=' ')

            # Read all data blocks
            table._read_all()
            lines = [fh.readline() for i in range(13)]


class Table(object):
    """Abstract superclass of all other classes for cross section tables.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. '92235.70c'.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    name : str
        ZAID identifier of the table, e.g. '92235.70c'.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in eV.

    """

    def __init__(self, name, atomic_weight_ratio, temperature):
        self.name = name
        self.atomic_weight_ratio = atomic_weight_ratio
        self.temperature = temperature

    def _read_all(self):
        raise NotImplementedError

    def _get_continuous_tabular(self, idx, ldis):
        """Get continuous tabular energy distribution (ACE law 4) starting at specified
        index in the XSS array.

        Parameters
        ----------
        idx : int
             Index in XSS array of the start of the energy distribution data
             (LDIS + LOCC - 1)
        ldis : int
             Index in XSS array of the start of the energy distribution block
             (e.g. JXS[11])

        Returns
        -------
        openmc.data.energy_distribution.ContinuousTabular
            Continuous tabular energy distribution

        """

        # Read number of interpolation regions and incoming energies
        n_regions = int(self._xss[idx])
        n_energy_in = int(self._xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = np.asarray(self._xss[idx:idx + n_regions], dtype=int)
            interpolation = np.asarray(self._xss[idx + n_regions:
                                                 idx + 2*n_regions], dtype=int)
        else:
            breakpoints = np.array([n_energy_in])
            interpolation = np.array([2])

        # Incoming energies at which distributions exist
        idx += 2 * n_regions + 1
        energy = self._xss[idx:idx + n_energy_in]

        # Location of distributions
        idx += n_energy_in
        loc_dist = np.asarray(self._xss[idx:idx + n_energy_in], dtype=int)

        # Initialize variables
        energy_out = []

        # Read each outgoing energy distribution
        for i in range(n_energy_in):
            idx = ldis + loc_dist[i] - 1

            # intt = interpolation scheme (1=hist, 2=lin-lin)
            INTTp = int(self._xss[idx])
            intt = INTTp % 10
            n_discrete_lines = (INTTp - intt)//10
            if intt not in (1, 2):
                warn("Interpolation scheme for continuous tabular distribution "
                     "is not histogram or linear-linear.")
                intt = 2

            n_energy_out = int(self._xss[idx + 1])
            data = self._xss[idx + 2:idx + 2 + 3*n_energy_out]
            data.shape = (3, n_energy_out)

            # Create continuous distribution
            eout_continuous = Tabular(data[0][n_discrete_lines:],
                             data[1][n_discrete_lines:],
                             interpolation_scheme[intt])
            eout_continuous.c = data[2][n_discrete_lines:]

            # If discrete lines are present, create a mixture distribution
            if n_discrete_lines > 0:
                eout_discrete = Discrete(data[0][:n_discrete_lines],
                                         data[1][:n_discrete_lines])
                eout_discrete.c = data[2][:n_discrete_lines]
                if n_discrete_lines == n_energy_out:
                    eout_i = eout_discrete
                else:
                    p_discrete = min(sum(eout_discrete.p), 1.0)
                    eout_i = Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            else:
                eout_i = eout_continuous

            energy_out.append(eout_i)

        return ContinuousTabular(breakpoints, interpolation, energy,
                                 energy_out)

    def _get_general_evaporation(self, idx):
        # Read nuclear temperature as Tabulated1D
        theta = _get_tabulated_1d(array, idx)

        # X-function
        nr = int(array[idx])
        ne = int(array[idx + 1 + 2*nr])
        idx += 2 + 2*nr + 2*ne
        net = int(array[idx])
        x = array[idx + 1:idx + 1 + net]

        raise NotImplementedError("Where'd you get this ACE file from?")

    def _get_maxwell_energy(self, idx):
        # Read nuclear temperature as Tabulated1D
        theta = _get_tabulated_1d(self._xss, idx)

        # Restriction energy
        nr = int(self._xss[idx])
        ne = int(self._xss[idx + 1 + 2*nr])
        u = self._xss[idx + 2 + 2*nr + 2*ne]

        return MaxwellEnergy(theta, u)

    def _get_evaporation(self, idx):
        # Read nuclear temperature as Tabulated1D
        theta = _get_tabulated_1d(self._xss, idx)

        # Restriction energy
        nr = int(self._xss[idx])
        ne = int(self._xss[idx + 1 + 2*nr])
        u = self._xss[idx + 2 + 2*nr + 2*ne]

        return Evaporation(theta, u)

    def _get_watt_energy(self, idx):
        # Energy-dependent a parameter
        a = _get_tabulated_1d(self._xss, idx)

        # Advance index
        nr = int(self._xss[idx])
        ne = int(self._xss[idx + 1 + 2*nr])
        idx += 2 + 2*nr + 2*ne

        # Energy-dependent b parameter
        b = _get_tabulated_1d(self._xss, idx)

        # Advance index
        nr = int(self._xss[idx])
        ne = int(self._xss[idx + 1 + 2*nr])
        idx += 2 + 2*nr + 2*ne

        # Restriction energy
        u = self._xss[idx]

        return WattEnergy(a, b, u)

    def _get_kalbach_mann(self, idx, ldis):
        # Read number of interpolation regions and incoming energies
        n_regions = int(self._xss[idx])
        n_energy_in = int(self._xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = np.asarray(self._xss[idx:idx + n_regions], dtype=int)
            interpolation = np.asarray(self._xss[idx + n_regions:
                                                 idx + 2*n_regions], dtype=int)
        else:
            breakpoints = np.array([n_energy_in])
            interpolation = np.array([2])

        # Incoming energies at which distributions exist
        idx += 2 * n_regions + 1
        energy = self._xss[idx:idx + n_energy_in]

        # Location of distributions
        idx += n_energy_in
        loc_dist = np.asarray(self._xss[idx:idx + n_energy_in], dtype=int)

        # Initialize variables
        energy_out = []
        km_r = []
        km_a = []

        # Read each outgoing energy distribution
        for i in range(n_energy_in):
            idx = ldis + loc_dist[i] - 1

            # intt = interpolation scheme (1=hist, 2=lin-lin)
            INTTp = int(self._xss[idx])
            intt = INTTp % 10
            n_discrete_lines = (INTTp - intt)//10
            if intt not in (1, 2):
                warn("Interpolation scheme for continuous tabular distribution "
                     "is not histogram or linear-linear.")
                intt = 2

            n_energy_out = int(self._xss[idx + 1])
            data = self._xss[idx + 2:idx + 2 + 5*n_energy_out]
            data.shape = (5, n_energy_out)

            # Create continuous distribution
            eout_continuous = Tabular(data[0][n_discrete_lines:],
                             data[1][n_discrete_lines:],
                             interpolation_scheme[intt])
            eout_continuous.c = data[2][n_discrete_lines:]

            # If discrete lines are present, create a mixture distribution
            if n_discrete_lines > 0:
                eout_discrete = Discrete(data[0][:n_discrete_lines],
                                         data[1][:n_discrete_lines])
                eout_discrete.c = data[2][:n_discrete_lines]
                if n_discrete_lines == n_energy_out:
                    eout_i = eout_discrete
                else:
                    p_discrete = min(sum(eout_discrete.p), 1.0)
                    eout_i = Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            else:
                eout_i = eout_continuous

            energy_out.append(eout_i)
            km_r.append(Tabulated1D(data[0], data[3]))
            km_a.append(Tabulated1D(data[0], data[4]))

        return KalbachMann(breakpoints, interpolation, energy, energy_out,
                           km_r, km_a)

    def _get_correlated(self, idx, ldis):
        # Read number of interpolation regions and incoming energies
        n_regions = int(self._xss[idx])
        n_energy_in = int(self._xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = np.asarray(self._xss[idx:idx + n_regions], dtype=int)
            interpolation = np.asarray(self._xss[idx + n_regions:
                                                 idx + 2*n_regions], dtype=int)
        else:
            breakpoints = np.array([n_energy_in])
            interpolation = np.array([2])

        # Incoming energies at which distributions exist
        idx += 2 * n_regions + 1
        energy = self._xss[idx:idx + n_energy_in]

        # Location of distributions
        idx += n_energy_in
        loc_dist = np.asarray(self._xss[idx:idx + n_energy_in], dtype=int)

        # Initialize list of distributions
        energy_out = []
        mu = []

        # Read each outgoing energy distribution
        for i in range(n_energy_in):
            idx = ldis + loc_dist[i] - 1

            # intt = interpolation scheme (1=hist, 2=lin-lin)
            INTTp = int(self._xss[idx])
            intt = INTTp % 10
            n_discrete_lines = (INTTp - intt)//10
            if intt not in (1, 2):
                warn("Interpolation scheme for continuous tabular distribution "
                     "is not histogram or linear-linear.")
                intt = 2

            # Secondary energy distribution
            n_energy_out = int(self._xss[idx + 1])
            data = self._xss[idx + 2:idx + 2 + 4*n_energy_out]
            data.shape = (4, n_energy_out)

            # Create continuous distribution
            eout_continuous = Tabular(data[0][n_discrete_lines:],
                                      data[1][n_discrete_lines:],
                                      interpolation_scheme[intt],
                                      ignore_negative=True)
            eout_continuous.c = data[2][n_discrete_lines:]

            # If discrete lines are present, create a mixture distribution
            if n_discrete_lines > 0:
                eout_discrete = Discrete(data[0][:n_discrete_lines],
                                         data[1][:n_discrete_lines])
                eout_discrete.c = data[2][:n_discrete_lines]
                if n_discrete_lines == n_energy_out:
                    eout_i = eout_discrete
                else:
                    p_discrete = min(sum(eout_discrete.p), 1.0)
                    eout_i = Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            else:
                eout_i = eout_continuous

            energy_out.append(eout_i)

            lc = np.asarray(data[3], dtype=int)

            # Secondary angular distributions
            mu_i = []
            for j in range(n_energy_out):
                if lc[j] > 0:
                    idx = ldis + abs(lc[j]) - 1

                    intt = int(self._xss[idx])
                    n_cosine = int(self._xss[idx + 1])
                    data = self._xss[idx + 2:idx + 2 + 3*n_cosine]
                    data.shape = (3, n_cosine)

                    mu_ij = Tabular(data[0], data[1], interpolation_scheme[intt])
                    mu_ij.c = data[2]
                else:
                    # Isotropic distribution
                    mu_ij = Uniform(-1., 1.)

                mu_i.append(mu_ij)

            # Add cosine distributions for this incoming energy to list
            mu.append(mu_i)

        return CorrelatedAngleEnergy(breakpoints, interpolation, energy,
                                     energy_out, mu)

    def _get_energy_distribution(self, location_dist, location_start, rx=None):
        """Returns an EnergyDistribution object from data read in starting at
        location_start.

        Parameters
        ----------
        location_dist : int
            Index in the XSS array corresponding to the start of a block,
            e.g. JXS(11) for the the DLW block.
        location_start : int
            Index in the XSS array corresponding to the start of an energy
            distribution array
        rx : Reaction
            Reaction this energy distribution will be associated with

        Returns
        -------
        distribution : openmc.data.AngleEnergy
            Secondary angle-energy distribution

        """

        # Set starting index for energy distribution
        idx = location_dist + location_start - 1

        law = int(self._xss[idx + 1])
        location_data = int(self._xss[idx + 2])

        # Position index for reading law data
        idx = location_dist + location_data - 1

        # Parse energy distribution data
        if law == 2:
            primary_flag = int(self._xss[idx])
            energy = self._xss[idx + 1]
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = DiscretePhoton(primary_flag, energy,
                                                 self.atomic_weight_ratio)
        elif law in (3, 33):
            threshold, mass_ratio = self._xss[idx:idx + 2]
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = LevelInelastic(threshold, mass_ratio)
        elif law == 4:
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = self._get_continuous_tabular(idx, location_dist)
        elif law == 5:
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = self._get_general_evaporation(idx)
        elif law == 7:
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = self._get_maxwell_energy(idx)
        elif law == 9:
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = self._get_evaporation(idx)
        elif law == 11:
            distribution = UncorrelatedAngleEnergy()
            distribution.energy = self._get_watt_energy(idx)
        elif law == 44:
            distribution = self._get_kalbach_mann(idx, location_dist)
        elif law == 61:
            distribution = self._get_correlated(idx, location_dist)
        elif law == 66:
            n_particles = int(self._xss[idx])
            total_mass = self._xss[idx + 1]
            distribution = NBodyPhaseSpace(total_mass, n_particles,
                                           self.atomic_weight_ratio, rx.Q_value)
        else:
            raise IOError("Unsupported ACE secondary energy "
                          "distribution law {0}".format(law))

        return distribution


class NeutronTable(Table):
    """A NeutronTable object contains continuous-energy neutron interaction data
    read from an ACE-formatted table. These objects are not normally
    instantiated by the user but rather created when reading data using a
    Library object and stored within the :attr:`Library.tables` attribute.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. '92235.70c'.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    absorption_xs : numpy.ndarray
        The microscopic absorption cross section for each value on the energy
        grid.
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide.
    energy : numpy.ndarray
        The energy values (MeV) at which reaction cross-sections are tabulated.
    heating_number : numpy.ndarray
        The total heating number for each value on the energy grid in MeV-b.
    name : str
        ZAID identifier of the table, e.g. 92235.70c.
    reactions : collections.OrderedDict
        Contains the cross sections, secondary angle and energy distributions,
        and other associated data for each reaction. The keys are the MT values
        and the values are Reaction objects.
    temperature : float
        Temperature of the target nuclide in eV.
    total_xs : numpy.ndarray
        The microscopic total cross section for each value on the energy grid in b.
    urr : None or openmc.data.ProbabilityTables
        Unresolved resonance region probability tables

    """

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(NeutronTable, self).__init__(name, atomic_weight_ratio, temperature)
        self.absorption_xs = None
        self.energy = None
        self.heating_number = None
        self.total_xs = None
        self.reactions = OrderedDict()
        self.urr = None

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Continuous-E Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Continuous-E Neutron Table>"

    def __iter__(self):
        return iter(self.reactions.values())

    def _read_all(self):
        self._read_cross_sections()
        self._read_nu()
        self._read_secondaries()
        self._read_photon_production_data()
        self._read_unr()

    def _read_cross_sections(self):
        """Read reaction cross sections and other data.

        Reads and parses the ESZ, MTR, LQR, TRY, LSIG, and SIG blocks. These
        blocks contain the energy grid, all reaction cross sections, the total
        cross section, average heating numbers, and a list of reactions with
        their Q-values and multiplicites.
        """

        # Determine number of energies on nuclide grid and number of reactions
        # excluding elastic scattering
        n_energies = self._nxs[3]
        n_reactions = self._nxs[4]

        # Read energy grid and total, absorption, elastic scattering, and
        # heating cross sections -- note that this appear separate from the rest
        # of the reaction cross sections
        arr = self._xss[self._jxs[1]:self._jxs[1] + 5*n_energies]
        arr.shape = (5, n_energies)
        self.energy, self.total_xs, self.absorption_xs, \
            elastic_xs, self.heating_number = arr

        # Create elastic scattering reaction
        elastic_scatter = Reaction(2, self)
        elastic_scatter.products.append(Product('neutron'))
        elastic_scatter.xs = Tabulated1D(self.energy, elastic_xs)
        self.reactions[2] = elastic_scatter

        # Create all other reactions with MT values
        mts = np.asarray(self._xss[self._jxs[3]:self._jxs[3] + n_reactions], dtype=int)
        qvalues = np.asarray(self._xss[self._jxs[4]:self._jxs[4] +
                                       n_reactions], dtype=float)
        tys = np.asarray(self._xss[self._jxs[5]:self._jxs[5] + n_reactions], dtype=int)

        # Create all reactions other than elastic scatter
        reactions = [(mt, Reaction(mt, self)) for mt in mts]
        self.reactions.update(reactions)

        # Loop over all reactions other than elastic scattering
        for i, rx in enumerate(list(self.reactions.values())[1:]):
            # Copy Q values determine if scattering should be treated in the
            # center-of-mass or lab system
            rx.Q_value = qvalues[i]
            rx.center_of_mass = (tys[i] < 0)

            # For neutron-producing reactions, get yield
            if rx.MT < 100:
                if tys[i] != 19:
                    if abs(tys[i]) > 100:
                        # Energy-dependent neutron yield
                        idx = self._jxs[11] + abs(tys[i]) - 101
                        yield_ = _get_tabulated_1d(self._xss, idx)
                    else:
                        yield_ = abs(tys[i])

                    neutron = Product('neutron')
                    neutron.yield_ = yield_
                    rx.products.append(neutron)

            # Get locator for cross-section data
            loc = int(self._xss[self._jxs[6] + i])

            # Determine starting index on energy grid
            rx.threshold_idx = int(self._xss[self._jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energies = int(self._xss[self._jxs[7] + loc])
            energy = self.energy[rx.threshold_idx:rx.threshold_idx + n_energies]

            # Read reaction cross section
            xs = self._xss[self._jxs[7] + loc + 1:self._jxs[7] + loc + 1 + n_energies]
            rx.xs = Tabulated1D(energy, xs)

    def _read_nu(self):
        """Read the NU block -- this contains information on the prompt and delayed
        neutron precursor yields, decay constants, etc

        """
        # No NU block
        if self._jxs[2] == 0:
            return

        products = []
        derived_products = []

        # Either prompt nu or total nu is given
        if self._xss[self._jxs[2]] > 0:
            whichnu = 'prompt' if self._jxs[24] > 0 else 'total'

            neutron = Product('neutron')
            neutron.emission_mode = whichnu

            idx = self._jxs[2]
            LNU = int(self._xss[idx])
            if LNU == 1:
                # Polynomial function form of nu
                NC = int(self._xss[idx+1])
                coefficients = self._xss[idx+2 : idx+2+NC]
                neutron.yield_ = Polynomial(coefficients)
            elif LNU == 2:
                # Tabular data form of nu
                neutron.yield_ = _get_tabulated_1d(self._xss, idx + 1)

            products.append(neutron)

        # Both prompt nu and total nu
        elif self._xss[self._jxs[2]] < 0:
            # Read prompt neutron yield
            prompt_neutron = Product('neutron')
            prompt_neutron.emission_mode = 'prompt'

            idx = self._jxs[2] + 1
            LNU = int(self._xss[idx])
            if LNU == 1:
                # Polynomial function form of nu
                NC = int(self._xss[idx+1])
                coefficients = self._xss[idx+2 : idx+2+NC]
                prompt_neutron.yield_ = Polynomial(coefficients)
            elif LNU == 2:
                # Tabular data form of nu
                prompt_neutron.yield_ = _get_tabulated_1d(self._xss, idx + 1)

            # Read total neutron yield
            total_neutron = Product('neutron')
            total_neutron.emission_mode = 'total'

            idx = self._jxs[2] + int(abs(self._xss[self._jxs[2]])) + 1
            LNU = int(self._xss[idx])

            if LNU == 1:
                # Polynomial function form of nu
                NC = int(self._xss[idx+1])
                coefficients = self._xss[idx+2 : idx+2+NC]
                total_neutron.yield_ = Polynomial(coefficients)
            elif LNU == 2:
                # Tabular data form of nu
                total_neutron.yield_ = _get_tabulated_1d(self._xss, idx + 1)

            products.append(prompt_neutron)
            derived_products.append(total_neutron)

        # Check for delayed nu data
        if self._jxs[24] > 0:
            yield_delayed = _get_tabulated_1d(self._xss, self._jxs[24] + 1)

            # Delayed neutron precursor distribution
            idx = self._jxs[25]
            n_group = self._nxs[8]
            total_group_probability = 0.
            for i, group in enumerate(range(n_group)):
                delayed_neutron = Product('neutron')
                delayed_neutron.emission_mode = 'delayed'
                delayed_neutron.decay_rate = self._xss[idx]

                group_probability = _get_tabulated_1d(self._xss, idx + 1)
                if np.all(group_probability.y == group_probability.y[0]):
                    delayed_neutron.yield_ = deepcopy(yield_delayed)
                    delayed_neutron.yield_.y *= group_probability.y[0]
                    total_group_probability += group_probability.y[0]
                else:
                    raise NotImplementedError(
                        'Delayed neutron with energy-dependent '
                        'group probability')

                # Advance position
                nr = int(self._xss[idx + 1])
                ne = int(self._xss[idx + 2 + 2*nr])
                idx += 3 + 2*nr + 2*ne

                # Energy distribution for delayed fission neutrons
                location_start = int(self._xss[self._jxs[26] + group])
                delayed_neutron.distribution.append(
                    self._get_energy_distribution(self._jxs[27], location_start))

                products.append(delayed_neutron)

            # Renormalize delayed neutron yields to reflect fact that in ACE
            # file, the sum of the group probabilities is not exactly one
            for product in products[1:]:
                product.yield_.y /= total_group_probability

        # Copy fission neutrons to reactions
        for MT, rx in self.reactions.items():
            if MT in (18, 19, 20, 21, 38):
                rx.products = deepcopy(products)
                if derived_products:
                    rx.derived_products = deepcopy(derived_products)

    def _get_angle_distribution(self, location_dist, location_start):
        # Set starting index for angle distribution
        idx = location_dist + location_start - 1

        # Number of energies at which angular distributions are tabulated
        n_energies = int(self._xss[idx])
        idx += 1

        # Incoming energy grid
        energy = self._xss[idx:idx + n_energies]
        idx += n_energies

        # Read locations for angular distributions
        lc = np.asarray(self._xss[idx:idx + n_energies], dtype=int)
        idx += n_energies

        mu = []
        for i in range(n_energies):
            if lc[i] > 0:
                # Equiprobable 32 bin distribution
                idx = location_dist + abs(lc[i]) - 1
                cos = self._xss[idx:idx + 33]
                pdf = np.zeros(33)
                pdf[:32] = 1.0/(32.0*np.diff(cos))
                cdf = np.linspace(0.0, 1.0, 33)

                mu_i = Tabular(cos, pdf, 'histogram', ignore_negative=True)
                mu_i.c = cdf
            elif lc[i] < 0:
                # Tabular angular distribution
                idx = location_dist + abs(lc[i]) - 1
                intt = int(self._xss[idx])
                n_points = int(self._xss[idx + 1])
                data = self._xss[idx + 2:idx + 2 + 3*n_points]
                data.shape = (3, n_points)

                mu_i = Tabular(data[0], data[1], interpolation_scheme[intt])
                mu_i.c = data[2]
            else:
                # Isotropic angular distribution
                mu_i = Uniform(-1., 1.)

            mu.append(mu_i)

        return AngleDistribution(energy, mu)

    def _read_secondaries(self):
        """Read angle/energy distributions for each reaction MT
        """

        # Number of reactions with secondary neutrons (including elastic
        # scattering)
        n_reactions = self._nxs[5] + 1

        for i, rx in enumerate(list(self.reactions.values())[:n_reactions]):
            if rx.MT == 18:
                for p in rx.products:
                    if p.emission_mode == 'prompt':
                        neutron = p
                        break
            else:
                neutron = rx.products[0]

            if i > 0:
                # Determine locator for ith energy distribution
                lnw = int(self._xss[self._jxs[10] + i - 1])

                while lnw > 0:
                    # Applicability of this distribution
                    neutron.applicability.append(_get_tabulated_1d(
                        self._xss, self._jxs[11] + lnw + 2))

                    # Read energy distribution data
                    neutron.distribution.append(self._get_energy_distribution(
                        self._jxs[11], lnw, rx))

                    lnw = int(self._xss[self._jxs[11] + lnw - 1])
            else:
                # No energy distribution for elastic scattering
                neutron.distribution.append(UncorrelatedAngleEnergy())

            # Check if angular distribution data exist
            loc = int(self._xss[self._jxs[8] + i])
            if loc == -1:
                # Angular distribution data are given as part of product
                # angle-energy distribution
                continue
            elif loc == 0:
                # No angular distribution data are given for this
                # reaction, isotropic scattering is asssumed
                angle_dist = None
            else:
                angle_dist = self._get_angle_distribution(self._jxs[9], loc)

            # Apply angular distribution to each uncorrelated angle-energy
            # distribution
            if angle_dist is not None:
                for d in neutron.distribution:
                    d.angle = angle_dist

    def _read_photon_production_data(self):
        """Read cross sections for each photon-production reaction"""

        n_photon_reactions = self._nxs[6]
        photon_mts = np.asarray(self._xss[self._jxs[13]:self._jxs[13] +
                                          n_photon_reactions], dtype=int)

        for i, rx in enumerate(photon_mts):
            # Determine corresponding reaction
            mt = photon_mts[i] // 1000
            reactions = []
            if mt not in self.reactions:
                # If the photon is assigned to MT=18 but the file splits fission
                # into MT=19,20,21,38, assign the photon product to each of the
                # individual reactions
                if mt == 18:
                    for mt_fiss in (19, 20, 21, 38):
                        if mt_fiss in self.reactions:
                            reactions.append(self.reactions[mt_fiss])
                if not reactions:
                    reactions.append(Reaction(mt, self))
            else:
                reactions.append(self.reactions[mt])

            # Create photon product and assign to reactions
            photon = Product('photon')
            for rx in reactions:
                rx.products.append(photon)

            # ==================================================================
            # Read photon yield / production cross section

            loca = int(self._xss[self._jxs[14] + i])
            idx = self._jxs[15] + loca - 1
            mftype = int(self._xss[idx])
            idx += 1

            if mftype in (12, 16):
                # Yield data taken from ENDF File 12 or 6
                mtmult = int(self._xss[idx])
                assert mtmult == mt

                # Read photon yield as function of energy
                photon.yield_ = _get_tabulated_1d(self._xss, idx + 1)

            elif mftype == 13:
                # Cross section data from ENDF File 13

                # Energy grid index at which data starts
                threshold_idx = int(self._xss[idx]) - 1

                # Get photon production cross section
                n_energy = int(self._xss[idx + 1])
                photon._xs = self._xss[idx + 2:idx + 2 + n_energy]

                # Determine yield based on ratio of cross sections
                energy = self.energy[threshold_idx:threshold_idx + n_energy]
                photon.yield_ = Tabulated1D(energy, photon._xs)

            else:
                raise ValueError("MFTYPE must be 12, 13, 16. Got {0}".format(
                        mftype))

            # ==================================================================
            # Read photon energy distribution

            location_start = int(self._xss[self._jxs[18] + i])

            # Read energy distribution data
            distribution = self._get_energy_distribution(
                self._jxs[19], location_start)
            assert isinstance(distribution, UncorrelatedAngleEnergy)

            # ==================================================================
            # Read photon angular distribution
            loc = int(self._xss[self._jxs[16] + i])

            if loc == 0:
                # No angular distribution data are given for this reaction,
                # isotropic scattering is asssumed in LAB
                energy = np.array([photon.yield_.x[0], photon.yield_.x[-1]])
                mu_isotropic = Uniform(-1., 1.)
                distribution.angle = AngleDistribution(
                    energy, [mu_isotropic, mu_isotropic])
            else:
                distribution.angle = self._get_angle_distribution(self._jxs[17], loc)

            # Add to list of distributions
            photon.distribution.append(distribution)

    def _read_unr(self):
        """Read the unresolved resonance range probability tables if present.
        """

        # Check if URR probability tables are present
        idx = self._jxs[23]
        if idx == 0:
            return

        N = int(self._xss[idx])      # Number of incident energies
        M = int(self._xss[idx+1])    # Length of probability table
        interpolation = int(self._xss[idx+2])
        inelastic_flag = int(self._xss[idx+3])
        absorption_flag = int(self._xss[idx+4])
        multiply_smooth = (int(self._xss[idx+5]) == 1)
        idx += 6

        # Get energies at which tables exist
        energy = self._xss[idx : idx+N]
        idx += N

        # Get probability tables
        table = self._xss[idx:idx+N*6*M]
        table.shape = (N, 6, M)

        # Create object
        self.urr = ProbabilityTables(energy, table, interpolation, inelastic_flag,
                                     absorption_flag, multiply_smooth)

    def export_to_hdf5(self, path, element=None, mass_number=None, metastable=0):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        element : str or None
            Elemental symbol, e.g. Zr. If not specified, the atomic
            number/symbol are inferred from the name of the table.
        mass_number : int or None
            Mass number of the nuclide. For natural elements, a value of zero
            should be given. If not specified, the mass number is inferred from
            the name of the table.
        metastable : int
            Metastable level of the nuclide. Defaults to 0.

        """

        f = h5py.File(path, 'a')

        # If element and/or mass number haven't been specified, make an educated
        # guess
        zaid, xs = self.name.split('.')
        if element is None:
            Z = int(zaid) // 1000
            element = atomic_symbol[Z]
        else:
            Z = atomic_number[element]
        if mass_number is None:
            mass_number = int(zaid) % 1000

        # Write basic data
        if metastable > 0:
            name = '{}{}_m{}.{}'.format(element, mass_number, metastable, xs)
        else:
            name = '{}{}.{}'.format(element, mass_number, xs)
        g = f.create_group(name)
        g.attrs['Z'] = Z
        g.attrs['A'] = mass_number
        g.attrs['metastable'] = metastable
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        g.attrs['temperature'] = self.temperature
        g.attrs['n_reaction'] = len(self.reactions)

        # Write energy grid
        g.create_dataset('energy', data=self.energy)

        # Write reaction data
        for i, rx in enumerate(self.reactions.values()):
            rx_group = g.create_group('reaction_{}'.format(i))
            rx.to_hdf5(rx_group)

            # Write total nu data if available
            if hasattr(rx, 'derived_products') and 'total_nu' not in g:
                tgroup = g.create_group('total_nu')
                rx.derived_products[0].to_hdf5(tgroup)

        # Write unresolved resonance probability tables
        if self.urr is not None:
            urr_group = g.create_group('urr')
            self.urr.to_hdf5(urr_group)

        f.close()

    @classmethod
    def from_hdf5(self, group):
        """Generate continuous-energy neutron interaction data from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group containing interaction data

        Returns
        -------
        openmc.data.ace.NeutronTable
            Continuous-energy neutron interaction data

        """
        name = group.name[1:]
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        temperature = group.attrs['temperature']
        table = NeutronTable(name, atomic_weight_ratio, temperature)

        # Read energy grid
        table.energy = group['energy'].value

        # Read reaction data
        n_reaction = group.attrs['n_reaction']

        # Write reaction data
        for i in range(n_reaction):
            rx_group = group['reaction_{}'.format(i)]
            rx = Reaction.from_hdf5(rx_group, table)
            table.reactions[rx.MT] = rx

            # Read total nu data if available
            if 'total_nu' in rx_group:
                tgroup = rx_group['total_nu']
                rx.derived_products = [Product.from_hdf5(tgroup)]

        # Read unresolved resonance probability tables
        if 'urr' in group:
            urr_group = group['urr']
            table.urr = ProbabilityTables.from_hdf5(urr_group)

        return table


class SabTable(Table):
    """A SabTable object contains thermal scattering data as represented by
    an S(alpha, beta) table.

    Parameters
    ----------
    name : str
        ZAID identifier of the table, e.g. lwtr.10t.
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    temperature : float
        Temperature of the target nuclide in eV.

    Attributes
    ----------
    atomic_weight_ratio : float
        Atomic mass ratio of the target nuclide.
    elastic_xs : openmc.data.Tabulated1D or openmc.data.CoherentElastic
        Elastic scattering cross section derived in the coherent or incoherent
        approximation
    inelastic_xs : openmc.data.Tabulated1D
        Inelastic scattering cross section derived in the incoherent
        approximation
    name : str
        ZAID identifier of the table, e.g. 92235.70c.
    temperature : float
        Temperature of the target nuclide in eV.

    """

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(SabTable, self).__init__(name, atomic_weight_ratio, temperature)
        self.elastic_xs = None
        self.elastic_mu_out = None

        self.inelastic_xs = None
        self.inelastic_e_out = None
        self.inelastic_mu_out = None
        self.secondary_mode = None

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Thermal S(a,b) Table: {0}>".format(self.name)
        else:
            return "<ACE Thermal S(a,b) Table>"

    def _read_all(self):
        self._read_itie()
        self._read_itce()
        self._read_itxe()
        self._read_itca()

    def _read_itie(self):
        """Read energy-dependent inelastic scattering cross sections.
        """
        idx = self._jxs[1]
        n_energies = int(self._xss[idx])
        energy = self._xss[idx+1 : idx+1+n_energies]
        xs = self._xss[idx+1+n_energies : idx+1+2*n_energies]
        self.inelastic_xs = Tabulated1D(energy, xs)

    def _read_itce(self):
        """Read energy-dependent elastic scattering cross sections.
        """
        # Determine if ITCE block exists
        idx = self._jxs[4]
        if idx == 0:
            return

        # Read values
        n_energies = int(self._xss[idx])
        energy = self._xss[idx+1 : idx+1+n_energies]
        P = self._xss[idx+1+n_energies : idx+1+2*n_energies]

        if self._nxs[5] == 4:
            self.elastic_xs = CoherentElastic(energy, P)
        else:
            self.elastic_xs = Tabulated1D(energy, P)

    def _read_itxe(self):
        """Read coupled energy/angle distributions for inelastic scattering.
        """
        # Determine number of energies and angles
        NE_in = len(self.inelastic_xs)
        NE_out = self._nxs[4]

        if self._nxs[7] == 0:
            self.secondary_mode = 'equal'
        elif self._nxs[7] == 1:
            self.secondary_mode = 'skewed'
        elif self._nxs[7] == 2:
            self.secondary_mode = 'continuous'

        if self.secondary_mode in ('equal', 'skewed'):
            NMU = self._nxs[3]
            idx = self._jxs[3]
            self.inelastic_e_out = self._xss[idx:idx+NE_in*NE_out*(NMU+2):NMU+2]
            self.inelastic_e_out.shape = (NE_in, NE_out)

            self.inelastic_mu_out = self._xss[idx:idx+NE_in*NE_out*(NMU+2)]
            self.inelastic_mu_out.shape = (NE_in, NE_out, NMU+2)
            self.inelastic_mu_out = self.inelastic_mu_out[:, :, 1:]
        else:
            NMU = self._nxs[3] - 1
            idx = self._jxs[3]
            locc = self._xss[idx:idx + NE_in].astype(int)
            NE_out = self._xss[idx + NE_in:idx + 2*NE_in].astype(int)
            energy_out = []
            mu_out = []
            for i in range(NE_in):
                idx = locc[i]

                # Outgoing energy distribution for incoming energy i
                e = self._xss[idx + 1:idx + 1 + NE_out[i]*(NMU + 3):NMU + 3]
                p = self._xss[idx + 2:idx + 2 + NE_out[i]*(NMU + 3):NMU + 3]
                c = self._xss[idx + 3:idx + 3 + NE_out[i]*(NMU + 3):NMU + 3]
                eout_i = Tabular(e, p, 'linear-linear', ignore_negative=True)
                eout_i.c = c

                # Outgoing angle distribution for each (incoming, outgoing) energy pair
                mu_i = []
                for j in range(NE_out[i]):
                    mu = self._xss[idx + 4:idx + 4 + NMU]
                    p_mu = 1./NMU*np.ones(NMU)
                    mu_ij = Discrete(mu, p_mu)
                    mu_ij.c = np.cumsum(p_mu)
                    mu_i.append(mu_ij)
                    idx += 3 + NMU

                energy_out.append(eout_i)
                mu_out.append(mu_i)

            # Create correlated angle-energy distribution
            breakpoints = [NE_in]
            interpolation = [2]
            energy = self.inelastic_xs.x
            self.inelastic_dist = CorrelatedAngleEnergy(
                breakpoints, interpolation, energy, energy_out, mu_out)

    def _read_itca(self):
        """Read angular distributions for elastic scattering.
        """
        NMU = self._nxs[6]
        if self._jxs[4] == 0 or NMU == -1:
            return
        idx = self._jxs[6]

        NE = len(self.elastic_xs)
        self.elastic_mu_out = self._xss[idx:idx+NE*NMU]
        self.elastic_mu_out.shape = (NE, NMU)

    def export_to_hdf5(self, path, name):
        """Export table to an HDF5 file.

        Parameters
        ----------
        path : str
            Path to write HDF5 file to
        name : str
            Name of compound (used as first group in HDF5 file)

        """

        f = h5py.File(path, 'a')

        # Write basic data
        g = f.create_group('{}.{}'.format(name, self.name.split('.')[1]))
        g.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio
        g.attrs['temperature'] = self.temperature
        g.attrs['zaids'] = self.zaids

        # Write thermal elastic scattering
        if self.elastic_xs is not None:
            elastic_group = g.create_group('elastic')
            self.elastic_xs.to_hdf5(elastic_group, 'xs')
            if self.elastic_mu_out is not None:
                elastic_group.create_dataset('mu_out', data=self.elastic_mu_out)

        # Write thermal inelastic scattering
        if self.inelastic_xs is not None:
            inelastic_group = g.create_group('inelastic')
            self.inelastic_xs.to_hdf5(inelastic_group, 'xs')
            inelastic_group.attrs['secondary_mode'] = np.string_(self.secondary_mode)
            if self.secondary_mode in ('equal', 'skewed'):
                inelastic_group.create_dataset('energy_out', data=self.inelastic_e_out)
                inelastic_group.create_dataset('mu_out', data=self.inelastic_mu_out)
            elif self.secondary_mode == 'continuous':
                self.inelastic_dist.to_hdf5(inelastic_group)

    @classmethod
    def from_hdf5(self, group):
        """Generate thermal scattering data from HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.SabTable
            Neutron thermal scattering data

        """
        name = group.name[1:]
        atomic_weight_ratio = group.attrs['atomic_weight_ratio']
        temperature = group.attrs['temperature']
        table = SabTable(name, atomic_weight_ratio, temperature)
        table.zaids = group.attrs['zaids']

        # Read thermal elastic scattering
        if 'elastic' in group:
            elastic_group = group['elastic']

            # Cross section
            elastic_xs_type = elastic_group['xs'].attrs['type'].decode()
            if elastic_xs_type == 'tab1':
                table.elastic_xs = Tabulated1D.from_hdf5(elastic_group['xs'])
            elif elastic_xs_type == 'bragg':
                table.elastic_xs = CoherentElastic.from_hdf5(elastic_group['xs'])

            # Angular distribution
            if 'mu_out' in elastic_group:
                table.elastic_mu_out = elastic_group['mu_out'].value

        # Read thermal inelastic scattering
        if 'inelastic' in group:
            inelastic_group = group['inelastic']
            table.secondary_mode = inelastic_group.attrs['secondary_mode'].decode()
            table.inelastic_xs = Tabulated1D.from_hdf5(inelastic_group['xs'])
            if table.secondary_mode in ('equal', 'skewed'):
                table.inelastic_e_out = inelastic_group['energy_out']
                table.inelastic_mu_out = inelastic_group['mu_out']
            elif table.secondary_mode == 'continuous':
                table.inelastic_dist = AngleEnergy.from_hdf5(inelastic_group)

        return table


class Reaction(object):
    """Reaction(MT, table=None)

    A Reaction object represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution. These objects are stored within the ``reactions`` attribute on
    subclasses of Table, e.g. NeutronTable.

    Parameters
    ----------
    MT : int
        The ENDF MT number for this reaction. On occasion, MCNP uses MT numbers
        that don't correspond exactly to the ENDF specification.
    table : openmc.data.ace.Table
        The ACE table which contains this reaction. This is useful if data on
        the parent nuclide is needed (for instance, the energy grid at which
        cross sections are tabulated)

    Attributes
    ----------
    center_of_mass : bool
        Indicates whether scattering kinematics should be performed in the
        center-of-mass or laboratory reference frame.
        grid above the threshold value in barns.
    MT : int
        The ENDF MT number for this reaction.
    Q_value : float
        The Q-value of this reaction in MeV.
    table : openmc.data.ace.Table
        The ACE table which contains this reaction.
    threshold : float
        Threshold of the reaction in MeV
    threshold_idx : int
        The index on the energy grid corresponding to the threshold of this
        reaction.
    xs : openmc.data.Tabulated1D
        Microscopic cross section for this reaction as a function of incident
        energy
    products : Iterable of openmc.data.Product
        Reaction products

    """

    def __init__(self, MT, table=None):
        self.center_of_mass = True
        self.table = table
        self.MT = MT
        self.Q_value = 0.
        self.threshold_idx = 0
        self._xs = None
        self.products = []

    def __repr__(self):
        if self.MT in reaction_name:
            return "<ACE Reaction: MT={} {}>".format(self.MT, reaction_name[self.MT])
        else:
            return "<ACE Reaction: MT={}>".format(self.MT)

    @property
    def center_of_mass(self):
        return self._center_of_mass

    @property
    def products(self):
        return self._products

    @property
    def threshold(self):
        return self.xs.x[0]

    @property
    def xs(self):
        return self._xs

    @center_of_mass.setter
    def center_of_mass(self, center_of_mass):
        cv.check_type('center of mass', center_of_mass, (bool, np.bool_))
        self._center_of_mass = center_of_mass

    @products.setter
    def products(self, products):
        cv.check_type('reaction products', products, Iterable, Product)
        self._products = products

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section', xs, Tabulated1D)
        for y in xs.y:
            cv.check_greater_than('reaction cross section', y, 0.0, True)
        self._xs = xs

    def to_hdf5(self, group):
        """Write reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['MT'] = self.MT
        if self.MT in reaction_name:
            group.attrs['label'] = np.string_(reaction_name[self.MT])
        else:
            group.attrs['label'] = np.string_(self.MT)
        group.attrs['Q_value'] = self.Q_value
        group.attrs['threshold_idx'] = self.threshold_idx + 1
        group.attrs['center_of_mass'] = 1 if self.center_of_mass else 0
        group.attrs['n_product'] = len(self.products)
        if self.xs is not None:
            group.create_dataset('xs', data=self.xs.y)
        for i, p in enumerate(self.products):
            pgroup = group.create_group('product_{}'.format(i))
            p.to_hdf5(pgroup)

    @classmethod
    def from_hdf5(cls, group, table):
        """Generate reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        Returns
        -------
        openmc.data.ace.Reaction
            Reaction data

        """
        MT = group.attrs['MT']
        rxn = cls(MT)
        rxn.table = table
        rxn.Q_value = group.attrs['Q_value']
        rxn.threshold_idx = group.attrs['threshold_idx'] - 1
        rxn.center_of_mass = bool(group.attrs['center_of_mass'])

        # Read cross section
        if 'xs' in group:
            xs = group['xs'].value
            rxn.xs = Tabulated1D(table.energy, xs)

        # Read reaction products
        n_product = group.attrs['n_product']
        products = []
        for i in range(n_product):
            pgroup = group['product_{}'.format(i)]
            products.append(Product.from_hdf5(pgroup))
        rxn.products = products

        return rxn


class DosimetryTable(Table):
    def __init__(self, name, atomic_weight_ratio, temperature):
        super(DosimetryTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Dosimetry Table: {0}>".format(self.name)
        else:
            return "<ACE Dosimetry Table>"


class NeutronDiscreteTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(NeutronDiscreteTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Discrete-E Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Discrete-E Neutron Table>"


class NeutronMGTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(NeutronMGTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Multigroup Neutron Table: {0}>".format(self.name)
        else:
            return "<ACE Multigroup Neutron Table>"


class PhotoatomicTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(PhotoatomicTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Continuous-E Photoatomic Table: {0}>".format(self.name)
        else:
            return "<ACE Continuous-E Photoatomic Table>"

    def _read_all(self):
        self._read_eszg()
        self._read_jinc()
        self._read_jcoh()
        self._read_heating()
        self._read_compton_data()

    def _read_eszg(self):
        # Determine number of energies on common energy grid
        n_energies = self._nxs[3]

        # Read cross sections
        idx = self._jxs[1]
        data = np.asarray(self._xss[idx:idx + 5*n_energies])
        data.shape = (5, n_energies)
        self.energy = data[0]
        self.incoherent = data[1]
        self.coherent = data[2]
        self.photoelectric = data[3]
        self.pairproduction = data[4]

    def _read_jinc(self):
        # Read incoherent scattering function
        idx = self._jxs[2]
        self.incoherent_scattering = self._xss[idx:idx + 21]

    def _read_jcoh(self):
        # Read coherent form factors and integrated coherent form factors
        idx = self._jxs[3]
        self.int_coherent_form_factors = self._xss[idx:idx + 55]
        self.coherent_form_factors = self._xss[idx + 55:idx + 2*55]

    def _read_jflo(self):
        raise NotImplementedError

    def _read_heating(self):
        idx = self._jxs[5]
        self.avg_heating = self._xss[idx:idx + self._nxs[3]]

    def _read_compton_data(self):
        # Determine number of Compton profiles
        n_shells = self._nxs[5]

        if n_shells > 0:
            # Number of electrons per shell
            idx = self._jxs[6]
            self.electrons_per_shell = np.asarray(
                self._xss[idx:idx + n_shells], dtype=int)

            # Binding energy per shell
            idx = self._jxs[7]
            self.binding_energy_per_shell = self._xss[idx:idx + n_shells]

            # Probability of interaction per shell
            idx = self._jxs[8]
            self.probability_per_shell = self._xss[idx:idx + n_shells]

            # Initialize arrays for Compton profile data
            self.compton_profile_interp = np.zeros(n_shells)
            self.compton_profile_momentum = []
            self.compton_profile_pdf = []
            self.compton_profile_cdf = []

            for i in range(n_shells):
                # Get locator for SWD block
                loca = int(self._xss[self._jxs[9] + i])
                idx = self._jxs[10] + loca - 1

                # Get interpolation parameter and number of momentum entries
                self.compton_profile_interp[i] = int(self._xss[idx])
                n_momentum = int(self._xss[idx + 1])
                idx += 2

                # Get momentum entries, PDF, and CDF
                data = self._xss[idx:idx + 3*n_momentum]
                data.shape = (3, n_momentum)
                self.compton_profile_momentum.append(data[0])
                self.compton_profile_pdf.append(data[1])
                self.compton_profile_cdf.append(data[2])


class PhotoatomicMGTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(PhotoatomicMGTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Multigroup Photoatomic Table: {0}>".format(self.name)
        else:
            return "<ACE Multigroup Photoatomic Table>"


class ElectronTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(ElectronTable, self).__init__(
            name, atomic_weight_ratio, temperature)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Electron Table: {0}>".format(self.name)
        else:
            return "<ACE Electron Table>"


class PhotonuclearTable(Table):

    def __init__(self, name, atomic_weight_ratio, temperature):
        super(PhotonuclearTable, self).__init__(
            name, atomic_weight_ratio, temperature)
        self.reactions = OrderedDict()

    def __repr__(self):
        if hasattr(self, 'name'):
            return "<ACE Photonuclear Table: {0}>".format(self.name)
        else:
            return "<ACE Photonuclear Table>"

    def _read_all(self):
        self._read_basic()
        self._read_cross_sections()
        self._read_secondaries()
        self._read_angular_distributions()
        self._read_energy_distributions()

    def _read_basic(self):
        n_energies = self._nxs[3]

        # Read energy mesh
        idx = self._jxs[1]
        self.energy = self._xss[idx:idx + n_energies]

        # Read total cross section
        idx = self._jxs[2]
        self.total_xs = self._xss[idx:idx + n_energies]

        # Read non-elastic and elastic cross section
        if self._jxs[4] > 0:
            idx = self._jxs[3]
            self.non_elastic_xs = self._xss[idx:idx + n_energies]
            idx = self._jxs[4]
            self.elastic_xs = self._xss[idx:idx + n_energies]
        else:
            self.non_elastic_xs = self.total_xs.copy()
            self.elastic_xs = np.zeros(n_energies)

        # Read heating numbers
        idx = self._jxs[5]
        if idx > 0:
            self.heating_number = self._xss[idx:idx + n_energies]
        else:
            self.heating_number = np.zeros(n_energies)

    def _read_cross_sections(self):
        # Determine number of reactions
        n_reactions = self._nxs[4]

        # Read list of MT numbers and Q values
        mts = np.asarray(self._xss[self._jxs[6]:self._jxs[6] +
                                  n_reactions], dtype=int)
        qvalues = np.asarray(self._xss[self._jxs[7]:self._jxs[7] +
                                      n_reactions])

        # Create reactions in dictionary
        reactions = [(mt, Reaction(mt, self)) for mt in mts]
        self.reactions.update(reactions)

        for i, rx in enumerate(self.reactions.values()):
            # Copy Q values
            rx.Q_value = qvalues[i]

            # Determine starting index on energy grid and number of energies
            idx = self._jxs[9] + int(self._xss[self._jxs[8] + i]) - 1
            rx.threshold_idx = int(self._xss[idx])
            n_energies = int(self._xss[idx + 1])
            energy = self.energy[rx.threshold_idx:rx.threshold_idx + n_energies]
            idx += 2

            # Read reaction cross setion
            xs = self._xss[idx:idx + n_energies]
            rx.xs = Tabulated1D(energy, xs, [], [])

    def _read_secondaries(self):
        names = {1: 'neutron', 2: 'photon', 3: 'electron',
                 9: 'proton', 31: 'deuteron', 32: 'triton',
                 33: 'helium3', 34: 'alpha'}

        n_particles = self._nxs[5]
        n_entries = self._nxs[7]

        idx = self._jxs[10]
        ixs = np.asarray(self._xss[idx:idx + n_particles*n_entries], dtype=int)
        ixs.shape = (n_particles, n_entries)
        self.ixs = ixs.transpose()

        self.particles = []

        for j in range(n_particles):
            # Create dictionary for particle
            particle = {}
            self.particles.append(particle)

            # Get secondary particle type/name
            particle['ipt'] = self.ixs[0, j]
            particle['name'] = names[particle['ipt']]

            # Number of reactions that produce secondary particle
            n_producing = self.ixs[1, j]

            # Particle-production cross section
            idx = self.ixs[2, j]
            particle['ie_production'] = int(self._xss[idx])
            ne = int(self._xss[idx + 1])
            idx += 2
            particle['production'] = self._xss[idx:idx + ne]

            # Average heating numbers
            idx = self.ixs[3, j]
            particle['ie_heating'] = int(self._xss[idx])
            ne = int(self._xss[idx + 1])
            idx += 2
            particle['heating_number'] = self._xss[idx:idx + ne]

            # MTs of particle production reactions
            idx = self.ixs[4, j]
            particle['mt_producing'] = np.asarray(
                self._xss[idx:idx + n_producing], dtype=int)

            # Coordinate system of reaction producing secondary particle
            idx = self.ixs[5, j]
            particle['center_of_mass'] = [i < 0 for i in
                                          self._xss[idx:idx + n_producing]]

            # Reaction yields
            particle['yield'] = {}
            for k in range(n_producing):
                # Create dictionary for yield data
                yieldData = {}

                # Read reaction yield data for a single MT
                idx = self.ixs[7, j] + int(self._xss[self.ixs[6, j] + k]) - 1

                yieldData['mftype'] = int(self._xss[idx])
                idx += 1

                if yieldData['mftype'] in (6, 12, 16):
                    # Yield data from ENDF File 6 or 12
                    mtmult = int(self._xss[idx])
                    assert mtmult == particle['mt_producing'][k]

                    # Read yield as function of energy
                    yieldData['multiplicity'] = _get_tabulated_1d(
                        self._xss, idx + 1)

                elif yieldData['mftype'] == 13:
                    # Production cross section for corresponding MT
                    yieldData['ie'] = int(self._xss[idx])
                    ne = int(self._xss[idx + 1])
                    idx += 2
                    yieldData['cross_section'] = self._xss[idx:idx + ne]

                # Add reaction yield data to dictionary
                mt = particle['mt_producing'][k]
                particle['yield'][mt] = yieldData

    def _read_angular_distributions(self):
        for j, particle in enumerate(self.particles):
            # Create dictionary for angular distributions
            angular_dists = {}
            particle['angular_distribution'] = angular_dists

            for k, mt in enumerate(particle['mt_producing']):
                landp = int(self._xss[self.ixs[8, j] + k])

                # check if angular distribution data exists
                if landp == -1:
                    # Angular distribution data are specified through the
                    # DLWP block
                    continue
                elif landp == 0:
                    # No angular distribution data are given for this
                    # reaction, isotropic scattering is assumed
                    ie = self.reactions[mt].threshold_idx
                    ne = len(self.reactions[mt].sigma)
                    angular_dists[mt] = AngularDistribution.isotropic(
                        np.array([self.energy[ie], self.energy[ie + ne - 1]]))
                    continue

            idx = self.ixs[9, j] + landp - 1

            angular_dists[mt] = AngularDistribution()
            angular_dists[mt].read(self._xss, idx, self.ixs[9, j])

    def _read_energy_distributions(self):
        for j, particle in enumerate(self.particles):
            # Create dictionary for energy distributions
            energy_dists = {}
            particle['energy_distribution'] = energy_dists

            for k, mt in enumerate(particle['mt_producing']):
                # Determine locator for kth energy distribution
                ldlwp = int(self._xss[self.ixs[10, j] + k])

                # Read energy distribution data
                energy_dists[mt] = self._get_energy_distribution(
                    self.ixs[11, j], ldlwp)
table_types = {
    "c": NeutronTable,
    "t": SabTable,
    "y": DosimetryTable,
    "d": NeutronDiscreteTable,
    "p": PhotoatomicTable,
    "m": NeutronMGTable,
    "g": PhotoatomicMGTable,
    "e": ElectronTable,
    "u": PhotonuclearTable}
