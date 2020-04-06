"""Module for parsing and manipulating data from ENDF evaluations.

All the classes and functions in this module are based on document
ENDF-102 titled "Data Formats and Procedures for the Evaluated Nuclear
Data File ENDF-6". The latest version from June 2009 can be found at
http://www-nds.iaea.org/ndspub/documents/endf/endf102/endf102.pdf

"""
import io
import re
import os
from math import pi
from pathlib import PurePath
from collections import OrderedDict
from collections.abc import Iterable

import numpy as np
from numpy.polynomial.polynomial import Polynomial

from .data import ATOMIC_SYMBOL, gnd_name
from .function import Tabulated1D, INTERPOLATION_SCHEME
from openmc.stats.univariate import Uniform, Tabular, Legendre
try:
    from ._endf import float_endf
    _CYTHON = True
except ImportError:
    _CYTHON = False


_LIBRARY = {0: 'ENDF/B', 1: 'ENDF/A', 2: 'JEFF', 3: 'EFF',
            4: 'ENDF/B High Energy', 5: 'CENDL', 6: 'JENDL',
            17: 'TENDL', 18: 'ROSFOND', 21: 'SG-21', 31: 'INDL/V',
            32: 'INDL/A', 33: 'FENDL', 34: 'IRDF', 35: 'BROND',
            36: 'INGDB-90', 37: 'FENDL/A', 41: 'BROND'}

_SUBLIBRARY = {
    0: 'Photo-nuclear data',
    1: 'Photo-induced fission product yields',
    3: 'Photo-atomic data',
    4: 'Radioactive decay data',
    5: 'Spontaneous fission product yields',
    6: 'Atomic relaxation data',
    10: 'Incident-neutron data',
    11: 'Neutron-induced fission product yields',
    12: 'Thermal neutron scattering data',
    19: 'Neutron standards',
    113: 'Electro-atomic data',
    10010: 'Incident-proton data',
    10011: 'Proton-induced fission product yields',
    10020: 'Incident-deuteron data',
    10030: 'Incident-triton data',
    20030: 'Incident-helion (3He) data',
    20040: 'Incident-alpha data'
}

SUM_RULES = {1: [2, 3],
             3: [4, 5, 11, 16, 17, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35,
                 36, 37, 41, 42, 44, 45, 152, 153, 154, 156, 157, 158, 159, 160,
                 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                 173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185,
                 186, 187, 188, 189, 190, 194, 195, 196, 198, 199, 200],
             4: list(range(50, 92)),
             16: list(range(875, 892)),
             18: [19, 20, 21, 38],
             27: [18, 101],
             101: [102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114,
                   115, 116, 117, 155, 182, 191, 192, 193, 197],
             103: list(range(600, 650)),
             104: list(range(650, 700)),
             105: list(range(700, 750)),
             106: list(range(750, 800)),
             107: list(range(800, 850))}

ENDF_FLOAT_RE = re.compile(r'([\s\-\+]?\d*\.\d+)([\+\-]) ?(\d+)')


def py_float_endf(s):
    """Convert string of floating point number in ENDF to float.

    The ENDF-6 format uses an 'e-less' floating point number format,
    e.g. -1.23481+10. Trying to convert using the float built-in won't work
    because of the lack of an 'e'. This function allows such strings to be
    converted while still allowing numbers that are not in exponential notation
    to be converted as well.

    Parameters
    ----------
    s : str
        Floating-point number from an ENDF file

    Returns
    -------
    float
        The number

    """
    return float(ENDF_FLOAT_RE.sub(r'\1e\2\3', s))


if not _CYTHON:
    float_endf = py_float_endf


def int_endf(s):
    """Convert string of integer number in ENDF to int.

    The ENDF-6 format technically allows integers to be represented by a field
    of all blanks. This function acts like int(s) except when s is a string of
    all whitespace, in which case zero is returned.

    Parameters
    ----------
    s : str
        Integer or spaces

    Returns
    -------
    integer
        The number or 0
    """
    return 0 if s.isspace() else int(s)


def get_text_record(file_obj):
    """Return data from a TEXT record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    str
        Text within the TEXT record

    """
    return file_obj.readline()[:66]


def get_cont_record(file_obj, skip_c=False):
    """Return data from a CONT record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from
    skip_c : bool
        Determine whether to skip the first two quantities (C1, C2) of the CONT
        record.

    Returns
    -------
    tuple
        The six items within the CONT record

    """
    line = file_obj.readline()
    if skip_c:
        C1 = None
        C2 = None
    else:
        C1 = float_endf(line[:11])
        C2 = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    N1 = int_endf(line[44:55])
    N2 = int_endf(line[55:66])
    return (C1, C2, L1, L2, N1, N2)


def get_head_record(file_obj):
    """Return data from a HEAD record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    tuple
        The six items within the HEAD record

    """
    line = file_obj.readline()
    ZA = int(float_endf(line[:11]))
    AWR = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    N1 = int_endf(line[44:55])
    N2 = int_endf(line[55:66])
    return (ZA, AWR, L1, L2, N1, N2)


def get_list_record(file_obj):
    """Return data from a LIST record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    list
        The six items within the header
    list
        The values within the list

    """
    # determine how many items are in list
    items = get_cont_record(file_obj)
    NPL = items[4]

    # read items
    b = []
    for i in range((NPL - 1)//6 + 1):
        line = file_obj.readline()
        n = min(6, NPL - 6*i)
        for j in range(n):
            b.append(float_endf(line[11*j:11*(j + 1)]))

    return (items, b)


def get_tab1_record(file_obj):
    """Return data from a TAB1 record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    list
        The six items within the header
    openmc.data.Tabulated1D
        The tabulated function

    """
    # Determine how many interpolation regions and total points there are
    line = file_obj.readline()
    C1 = float_endf(line[:11])
    C2 = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    n_regions = int_endf(line[44:55])
    n_pairs = int_endf(line[55:66])
    params = [C1, C2, L1, L2]

    # Read the interpolation region data, namely NBT and INT
    breakpoints = np.zeros(n_regions, dtype=int)
    interpolation = np.zeros(n_regions, dtype=int)
    m = 0
    for i in range((n_regions - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_regions - m)
        for j in range(to_read):
            breakpoints[m] = int_endf(line[0:11])
            interpolation[m] = int_endf(line[11:22])
            line = line[22:]
            m += 1

    # Read tabulated pairs x(n) and y(n)
    x = np.zeros(n_pairs)
    y = np.zeros(n_pairs)
    m = 0
    for i in range((n_pairs - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_pairs - m)
        for j in range(to_read):
            x[m] = float_endf(line[:11])
            y[m] = float_endf(line[11:22])
            line = line[22:]
            m += 1

    return params, Tabulated1D(x, y, breakpoints, interpolation)


def get_tab2_record(file_obj):
    # Determine how many interpolation regions and total points there are
    params = get_cont_record(file_obj)
    n_regions = params[4]

    # Read the interpolation region data, namely NBT and INT
    breakpoints = np.zeros(n_regions, dtype=int)
    interpolation = np.zeros(n_regions, dtype=int)
    m = 0
    for i in range((n_regions - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_regions - m)
        for j in range(to_read):
            breakpoints[m] = int(line[0:11])
            interpolation[m] = int(line[11:22])
            line = line[22:]
            m += 1

    return params, Tabulated2D(breakpoints, interpolation)


def get_intg_record(file_obj):
    """
    Return data from an INTG record in an ENDF-6 file. Used to store the
    covariance matrix in a compact format.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    numpy.ndarray
        The correlation matrix described in the INTG record
    """
    # determine how many items are in list and NDIGIT
    items = get_cont_record(file_obj)
    ndigit = items[2]
    npar = items[3]    # Number of parameters
    nlines = items[4]  # Lines to read
    NROW_RULES = {2: 18, 3: 12, 4: 11, 5: 9, 6: 8}
    nrow = NROW_RULES[ndigit]

    # read lines and build correlation matrix
    corr = np.identity(npar)
    for i in range(nlines):
        line = file_obj.readline()
        ii = int_endf(line[:5]) - 1  # -1 to account for 0 indexing
        jj = int_endf(line[5:10]) - 1
        factor = 10**ndigit
        for j in range(nrow):
            if jj+j >= ii:
                break
            element = int_endf(line[11+(ndigit+1)*j:11+(ndigit+1)*(j+1)])
            if element > 0:
                corr[ii, jj] = (element+0.5)/factor
            elif element < 0:
                corr[ii, jj] = (element-0.5)/factor

    # Symmetrize the correlation matrix
    corr = corr + corr.T - np.diag(corr.diagonal())
    return corr


def get_evaluations(filename):
    """Return a list of all evaluations within an ENDF file.

    Parameters
    ----------
    filename : str
        Path to ENDF-6 formatted file

    Returns
    -------
    list
        A list of :class:`openmc.data.endf.Evaluation` instances.

    """
    evaluations = []
    with open(str(filename), 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if line[66:70] == '  -1':
                break
            fh.seek(pos)
            evaluations.append(Evaluation(fh))
    return evaluations


class Evaluation:
    """ENDF material evaluation with multiple files/sections

    Parameters
    ----------
    filename_or_obj : str or file-like
        Path to ENDF file to read or an open file positioned at the start of an
        ENDF material

    Attributes
    ----------
    info : dict
        Miscellaneous information about the evaluation.
    target : dict
        Information about the target material, such as its mass, isomeric state,
        whether it's stable, and whether it's fissionable.
    projectile : dict
        Information about the projectile such as its mass.
    reaction_list : list of 4-tuples
        List of sections in the evaluation. The entries of the tuples are the
        file (MF), section (MT), number of records (NC), and modification
        indicator (MOD).

    """
    def __init__(self, filename_or_obj):
        if isinstance(filename_or_obj, (str, PurePath)):
            fh = open(str(filename_or_obj), 'r')
        else:
            fh = filename_or_obj
        self.section = {}
        self.info = {}
        self.target = {}
        self.projectile = {}
        self.reaction_list = []

        # Skip TPID record. Evaluators sometimes put in TPID records that are
        # ill-formated because they lack MF/MT values or put them in the wrong
        # columns.
        if fh.tell() == 0:
            fh.readline()
        MF = 0

        # Determine MAT number for this evaluation
        while MF == 0:
            position = fh.tell()
            line = fh.readline()
            MF = int(line[70:72])
        self.material = int(line[66:70])
        fh.seek(position)

        while True:
            # Find next section
            while True:
                position = fh.tell()
                line = fh.readline()
                MAT = int(line[66:70])
                MF = int(line[70:72])
                MT = int(line[72:75])
                if MT > 0 or MAT == 0:
                    fh.seek(position)
                    break

            # If end of material reached, exit loop
            if MAT == 0:
                fh.readline()
                break

            section_data = ''
            while True:
                line = fh.readline()
                if line[72:75] == '  0':
                    break
                else:
                    section_data += line
            self.section[MF, MT] = section_data

        self._read_header()

    def __repr__(self):
        if 'zsymam' in self.target:
            name = self.target['zsymam'].replace(' ', '')
        else:
            name = 'Unknown'
        return '<{} for {} {}>'.format(self.info['sublibrary'], name,
                                       self.info['library'])

    def _read_header(self):
        file_obj = io.StringIO(self.section[1, 451])

        # Information about target/projectile
        items = get_head_record(file_obj)
        Z, A = divmod(items[0], 1000)
        self.target['atomic_number'] = Z
        self.target['mass_number'] = A
        self.target['mass'] = items[1]
        self._LRP = items[2]
        self.target['fissionable'] = (items[3] == 1)
        try:
            library = _LIBRARY[items[4]]
        except KeyError:
            library = 'Unknown'
        self.info['modification'] = items[5]

        # Control record 1
        items = get_cont_record(file_obj)
        self.target['excitation_energy'] = items[0]
        self.target['stable'] = (int(items[1]) == 0)
        self.target['state'] = items[2]
        self.target['isomeric_state'] = m = items[3]
        self.info['format'] = items[5]
        assert self.info['format'] == 6

        # Set correct excited state for Am242_m1, which is wrong in ENDF/B-VII.1
        if Z == 95 and A == 242 and m == 1:
            self.target['state'] = 2

        # Control record 2
        items = get_cont_record(file_obj)
        self.projectile['mass'] = items[0]
        self.info['energy_max'] = items[1]
        library_release = items[2]
        self.info['sublibrary'] = _SUBLIBRARY[items[4]]
        library_version = items[5]
        self.info['library'] = (library, library_version, library_release)

        # Control record 3
        items = get_cont_record(file_obj)
        self.target['temperature'] = items[0]
        self.info['derived'] = (items[2] > 0)
        NWD = items[4]
        NXC = items[5]

        # Text records
        text = [get_text_record(file_obj) for i in range(NWD)]
        if len(text) >= 5:
            self.target['zsymam'] = text[0][0:11]
            self.info['laboratory'] = text[0][11:22]
            self.info['date'] = text[0][22:32]
            self.info['author'] = text[0][32:66]
            self.info['reference'] = text[1][1:22]
            self.info['date_distribution'] = text[1][22:32]
            self.info['date_release'] = text[1][33:43]
            self.info['date_entry'] = text[1][55:63]
            self.info['identifier'] = text[2:5]
            self.info['description'] = text[5:]

        # File numbers, reaction designations, and number of records
        for i in range(NXC):
            _, _, mf, mt, nc, mod = get_cont_record(file_obj, skip_c=True)
            self.reaction_list.append((mf, mt, nc, mod))

    @property
    def gnd_name(self):
        return gnd_name(self.target['atomic_number'],
                        self.target['mass_number'],
                        self.target['isomeric_state'])


class Tabulated2D:
    """Metadata for a two-dimensional function.

    This is a dummy class that is not really used other than to store the
    interpolation information for a two-dimensional function. Once we refactor
    to adopt GND-like data containers, this will probably be removed or
    extended.

    Parameters
    ----------
    breakpoints : Iterable of int
        Breakpoints for interpolation regions
    interpolation : Iterable of int
        Interpolation scheme identification number, e.g., 3 means y is linear in
        ln(x).

    """
    def __init__(self, breakpoints, interpolation):
        self.breakpoints = breakpoints
        self.interpolation = interpolation
