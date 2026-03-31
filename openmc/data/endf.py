"""Module for parsing and manipulating data from ENDF evaluations.

All the classes and functions in this module are based on document ENDF-102
titled "Data Formats and Procedures for the Evaluated Nuclear Data File ENDF-6".
The version from September 2023 can be found at
https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf

"""
import io
from pathlib import PurePath
import re

from .data import gnds_name
from .function import Tabulated1D
from endf.material import _LIBRARY, _SUBLIBRARY, get_materials as get_evaluations
from endf.incident_neutron import SUM_RULES
from endf.records import (
    float_endf,
    py_float_endf,
    int_endf,
    get_text_record,
    get_cont_record,
    get_head_record,
    get_list_record,
    get_tab1_record as _get_tab1_record,
    get_tab2_record,
    get_intg_record,
)


def get_tab1_record(file_obj):
    """Return data from a TAB1 record in an ENDF-6 file.

    This wraps the endf package's get_tab1_record to return an
    openmc.data.Tabulated1D (which has HDF5 support and is a Function1D)
    instead of endf.Tabulated1D.
    """
    params, tab = _get_tab1_record(file_obj)
    return params, Tabulated1D(tab.x, tab.y, tab.breakpoints, tab.interpolation)


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
            need_to_close = True
        else:
            fh = filename_or_obj
            need_to_close = False
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

        if need_to_close:
            fh.close()

        self._read_header()

    def __repr__(self):
        name = self.target['zsymam'].replace(' ', '')
        return f"<{self.info['sublibrary']} for {name} {self.info['library']}>"

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
        else:
            self.target['zsymam'] = 'Unknown'

        # File numbers, reaction designations, and number of records
        for i in range(NXC):
            _, _, mf, mt, nc, mod = get_cont_record(file_obj, skip_c=True)
            self.reaction_list.append((mf, mt, nc, mod))

    @property
    def gnds_name(self):
        return gnds_name(self.target['atomic_number'],
                         self.target['mass_number'],
                         self.target['isomeric_state'])

