from collections.abc import Iterable, Callable, MutableMapping
from copy import deepcopy
from dataclasses import dataclass, field
from io import StringIO
import math
from numbers import Real
from warnings import warn

import h5py
import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
from openmc.stats import Uniform, Tabular, Legendre
from .angle_distribution import AngleDistribution
from .angle_energy import AngleEnergy
from .correlated import CorrelatedAngleEnergy
from .data import ATOMIC_SYMBOL, K_BOLTZMANN, EV_PER_MEV
from .endf import as_evaluation, get_head_record, get_tab1_record, \
    get_list_record, get_tab2_record, get_cont_record
from .energy_distribution import EnergyDistribution, LevelInelastic, \
    DiscretePhoton
from .function import Tabulated1D, Polynomial
from .kalbach_mann import KalbachMann
from .laboratory import LaboratoryAngleEnergy
from .nbody import NBodyPhaseSpace
from .photon import _SUBSHELLS
from .product import Product
from .uncorrelated import UncorrelatedAngleEnergy


REACTION_NAME = {1: '(n,total)', 2: '(n,elastic)', 3: "(n,nonelastic)",
                 4: '(n,level)', 5: '(n,misc)', 11: '(n,2nd)', 16: '(n,2n)',
                 17: '(n,3n)', 18: '(n,fission)', 19: '(n,f)', 20: '(n,nf)',
                 21: '(n,2nf)', 22: '(n,na)', 23: '(n,n3a)', 24: '(n,2na)',
                 25: '(n,3na)', 27: '(n,absorption)', 28: '(n,np)', 29: '(n,n2a)',
                 30: '(n,2n2a)', 32: '(n,nd)', 33: '(n,nt)', 34: '(n,n3He)',
                 35: '(n,nd2a)', 36: '(n,nt2a)', 37: '(n,4n)', 38: '(n,3nf)',
                 41: '(n,2np)', 42: '(n,3np)', 44: '(n,n2p)', 45: '(n,npa)',
                 91: '(n,nc)', 101: '(n,disappear)', 102: '(n,gamma)',
                 103: '(n,p)', 104: '(n,d)', 105: '(n,t)', 106: '(n,3He)',
                 107: '(n,a)', 108: '(n,2a)', 109: '(n,3a)', 111: '(n,2p)',
                 112: '(n,pa)', 113: '(n,t2a)', 114: '(n,d2a)', 115: '(n,pd)',
                 116: '(n,pt)', 117: '(n,da)', 152: '(n,5n)', 153: '(n,6n)',
                 154: '(n,2nt)', 155: '(n,ta)', 156: '(n,4np)', 157: '(n,3nd)',
                 158: '(n,nda)', 159: '(n,2npa)', 160: '(n,7n)', 161: '(n,8n)',
                 162: '(n,5np)', 163: '(n,6np)', 164: '(n,7np)', 165: '(n,4na)',
                 166: '(n,5na)', 167: '(n,6na)', 168: '(n,7na)', 169: '(n,4nd)',
                 170: '(n,5nd)', 171: '(n,6nd)', 172: '(n,3nt)', 173: '(n,4nt)',
                 174: '(n,5nt)', 175: '(n,6nt)', 176: '(n,2n3He)',
                 177: '(n,3n3He)', 178: '(n,4n3He)', 179: '(n,3n2p)',
                 180: '(n,3n2a)', 181: '(n,3npa)', 182: '(n,dt)',
                 183: '(n,npd)', 184: '(n,npt)', 185: '(n,ndt)',
                 186: '(n,np3He)', 187: '(n,nd3He)', 188: '(n,nt3He)',
                 189: '(n,nta)', 190: '(n,2n2p)', 191: '(n,p3He)',
                 192: '(n,d3He)', 193: '(n,3Hea)', 194: '(n,4n2p)',
                 195: '(n,4n2a)', 196: '(n,4npa)', 197: '(n,3p)',
                 198: '(n,n3p)', 199: '(n,3n2pa)', 200: '(n,5n2p)', 203: '(n,Xp)',
                 204: '(n,Xd)', 205: '(n,Xt)', 206: '(n,X3He)', 207: '(n,Xa)',
                 301: 'heating', 444: 'damage-energy',
                 501: 'photon-total', 502: 'coherent-scatter',
                 504: 'incoherent-scatter', 515: 'pair-production-electron',
                 516: 'pair-production', 517: 'pair-production-nuclear',
                 522: 'photoelectric',
                 649: '(n,pc)', 699: '(n,dc)', 749: '(n,tc)', 799: '(n,3Hec)',
                 849: '(n,ac)', 891: '(n,2nc)', 901: 'heating-local'}
REACTION_NAME.update({i: f'(n,n{i - 50})' for i in range(51, 91)})
REACTION_NAME.update({i: f'(n,p{i - 600})' for i in range(600, 649)})
REACTION_NAME.update({i: f'(n,d{i - 650})' for i in range(650, 699)})
REACTION_NAME.update({i: f'(n,t{i - 700})' for i in range(700, 749)})
REACTION_NAME.update({i: f'(n,3He{i - 750})' for i in range(750, 799)})
REACTION_NAME.update({i: f'(n,a{i - 800})' for i in range(800, 849)})
REACTION_NAME.update({i: f'(n,2n{i - 875})' for i in range(875, 891)})
REACTION_NAME.update(
    {534 + i: f'photoelectric-{shell}' for i, shell in enumerate(_SUBSHELLS[1:])}
)

REACTION_MT = {name: mt for mt, name in REACTION_NAME.items()}
REACTION_MT['total'] = 1
REACTION_MT['elastic'] = 2
REACTION_MT['fission'] = 18
REACTION_MT['absorption'] = 27
REACTION_MT['capture'] = 102

FISSION_MTS = (18, 19, 20, 21, 38)


_METADATA_SEMANTIC_STATUS = (
    'complete', 'unsupported_lang', 'unsupported_ltp', 'unsupported_frame',
    'invalid_combination')
_METADATA_STRUCTURAL_FAILURE = (
    'unknown_law', 'unproven_cursor', 'truncated_record', 'count_mismatch',
    'premature_end', 'trailing_record', 'invalid_header',
    'invalid_law6_control')
_METADATA_JP = (0, 1, 2, 10, 11, 12, 20, 21, 22)
_METADATA_SCHEMA = 'openmc-evaluated-product-metadata/v2'
_METADATA_RESULT = 'published_structurally_complete'
_METADATA_UTF8 = h5py.string_dtype(encoding='utf-8')


class _MetadataParseError(ValueError):
    """Internal typed terminal outcome from the MF=6 cursor reader."""

    def __init__(self, reason):
        self.reason = reason
        super().__init__(reason)


def _metadata_declared_record_lines(reader, line):
    """Return the number of physical lines declared by a record header."""
    try:
        n1 = int(line[44:55].strip() or 0)
        n2 = int(line[55:66].strip() or 0)
    except ValueError:
        return None
    if n1 < 0 or n2 < 0:
        return (-1 if reader in (
            get_list_record, get_tab1_record, get_tab2_record) else 1)
    if reader is get_list_record:
        return 1 + (n1 + 5)//6
    if reader is get_tab1_record:
        return 1 + (n1 + 2)//3 + (n2 + 2)//3
    if reader is get_tab2_record:
        return 1 + (n1 + 2)//3
    return 1


def _metadata_record(reader, file_obj, failure='truncated_record',
                     eof_failure=None, parse_failure=None):
    """Read one declared record envelope or raise a typed terminal outcome."""
    if file_obj.tell() >= len(file_obj.getvalue()):
        raise _MetadataParseError(eof_failure or failure)
    remaining = file_obj.getvalue()[file_obj.tell():]
    lines = remaining.splitlines(keepends=True)
    expected_lines = _metadata_declared_record_lines(reader, lines[0])
    if expected_lines == -1:
        raise _MetadataParseError('count_mismatch')
    if (expected_lines is not None and
            (len(lines) < expected_lines or any(
                len(line.rstrip('\r\n')) < 66
                for line in lines[:expected_lines]))):
        raise _MetadataParseError(failure)
    try:
        return reader(file_obj)
    except (EOFError, IndexError, OSError, ValueError, TypeError) as exc:
        raise _MetadataParseError(parse_failure or failure) from exc


@dataclass(frozen=True)
class EvaluatedProductMetadataEntry:
    """Immutable evaluated MF=6 product-provenance entry.

    The attributes are evaluated-data metadata only. They are not an OpenMC
    transport product and do not describe a sampled final state.
    """

    subsection_index: int
    zap: int
    awp: float
    awp_interpretation: str
    lip: int
    law: int
    lang_values: tuple
    lidp: int | None
    ltp_values: tuple
    semantic_status: str
    identity_status: str
    derived_particle: str | None
    apsx: float | None = None
    npsx: int | None = None

    def __post_init__(self):
        """Require LAW=6 controls only for LAW=6 entries."""
        if self.law == 6:
            if self.apsx is None or self.npsx is None:
                raise ValueError('LAW=6 metadata requires APSX and NPSX')
            if (isinstance(self.apsx, (bool, np.bool_)) or
                    not isinstance(self.apsx, Real)):
                raise ValueError('invalid LAW=6 metadata APSX')
            apsx = float(self.apsx)
            if not math.isfinite(apsx) or apsx <= 0.0:
                raise ValueError('invalid LAW=6 metadata APSX')
            if (isinstance(self.npsx, (bool, np.bool_)) or
                    not isinstance(self.npsx, (int, np.integer)) or
                    not 3 <= self.npsx <= 2**31 - 1):
                raise ValueError('invalid LAW=6 metadata NPSX')
        elif self.apsx is not None or self.npsx is not None:
            raise ValueError('non-LAW=6 metadata cannot contain APSX or NPSX')


@dataclass(frozen=True)
class EvaluatedProductMetadata:
    """Immutable, transport-isolated published MF=6 section metadata."""

    mt: int
    jp: int
    lct: int
    section_semantic_status: str
    entries: tuple


@dataclass(frozen=True)
class EvaluatedProductMetadataResult:
    """Result of direct MF=6 metadata parsing before any attachment.

    A structural failure and an absent MF=6 section intentionally contain no
    published metadata. This prevents a partially consumed cursor from being
    attached to an ACE-derived reaction.
    """

    mt: int
    result_status: str
    structural_failure: str | None = None
    jp: int | None = None
    lct: int | None = None
    section_semantic_status: str | None = None
    hdf5_group_version: int | None = None
    entries: tuple = ()
    schema: str = field(default=_METADATA_SCHEMA, init=False)
    source_format: str = field(default='ENDF-6', init=False)
    mf: int = field(default=6, init=False)

    def __post_init__(self):
        if (isinstance(self.mt, (bool, np.bool_)) or
                not isinstance(self.mt, (int, np.integer)) or
                not 1 <= self.mt <= 999):
            raise ValueError('invalid evaluated product metadata result MT')
        empty = (self.jp is None and self.lct is None and
                 self.section_semantic_status is None and
                 self.hdf5_group_version is None and not self.entries)
        if self.result_status == 'absent':
            if self.structural_failure is not None or not empty:
                raise ValueError('invalid absent metadata result')
        elif self.result_status == 'structural_failure':
            if self.structural_failure not in _METADATA_STRUCTURAL_FAILURE:
                raise ValueError('invalid metadata structural failure')
            if not empty:
                raise ValueError('invalid structural-failure metadata result')
        elif self.result_status == _METADATA_RESULT:
            if (self.structural_failure is not None or
                    self.hdf5_group_version != 1):
                raise ValueError('invalid published metadata result')
            _metadata_validate_published(self.metadata)
        else:
            raise ValueError('invalid evaluated product metadata result status')

    @property
    def metadata(self):
        """Return the attachable object only for a published result."""
        if self.result_status != _METADATA_RESULT:
            return None
        return EvaluatedProductMetadata(
            self.mt, self.jp, self.lct, self.section_semantic_status,
            self.entries)


def _metadata_int(value, failure='count_mismatch'):
    """Return an ENDF control integer or raise for a non-integral value."""
    try:
        valid = not isinstance(value, (bool, np.bool_)) and int(value) == value
    except (OverflowError, TypeError, ValueError):
        valid = False
    if not valid:
        raise _MetadataParseError(failure)
    return int(value)


def _metadata_status(statuses):
    """Return the canonical maximum semantic status."""
    return max(statuses, key=_METADATA_SEMANTIC_STATUS.index, default='complete')


def _metadata_identity(zap, lip):
    """Derive only unambiguous convenience identities from raw ZAP/LIP."""
    if zap == 0:
        return 'special_particle', 'photon'
    if zap == 1:
        return 'special_particle', 'neutron'
    return 'raw_only', None


def _metadata_identity_is_valid(zap, lip, status, derived):
    """Return whether a raw or optionally derived identity is conservative."""
    if zap in (0, 1):
        return (status, derived) == _metadata_identity(zap, lip)
    if status == 'raw_only':
        return derived is None
    if status != 'ground_state_nuclide' or derived is None or lip != 0:
        return False
    z, a = divmod(zap, 1000)
    return (zap != 1000 and 1 <= z <= a < 1000 and
            z < len(ATOMIC_SYMBOL) and derived == f'{ATOMIC_SYMBOL[z]}{a}')


def _metadata_entry(index, params):
    """Create a raw entry and its mutable semantic-reason collection."""
    try:
        awp = float(params[1])
    except (OverflowError, TypeError, ValueError) as exc:
        raise _MetadataParseError('invalid_header') from exc
    zap, lip, law = (_metadata_int(params[0], 'invalid_header'),
                     _metadata_int(params[2], 'invalid_header'),
                     _metadata_int(params[3], 'invalid_header'))
    if not math.isfinite(awp) or awp < 0.0:
        raise _MetadataParseError('invalid_header')
    identity_status, derived_particle = _metadata_identity(zap, lip)
    return {
        'subsection_index': index, 'zap': zap, 'awp': awp, 'lip': lip,
        'law': law, 'lang_values': [], 'lidp': None, 'ltp_values': [],
        'semantic_statuses': [], 'identity_status': identity_status,
        'derived_particle': derived_particle, 'apsx': None, 'npsx': None,
    }


def _metadata_law1(file_obj, entry, lct):
    params, _ = _metadata_record(get_tab2_record, file_obj)
    lang, ne = _metadata_int(params[2]), _metadata_int(params[5])
    if ne < 0:
        raise _MetadataParseError('count_mismatch')
    entry['lang_values'].append(lang)
    for _ in range(ne):
        items, values = _metadata_record(get_list_record, file_obj)
        nd, na, nw, nep = (_metadata_int(items[2]), _metadata_int(items[3]),
                           _metadata_int(items[4]), _metadata_int(items[5]))
        if (min(nd, na, nw, nep) < 0 or nd > nep or
                nw != nep*(na + 2) or len(values) != nw):
            raise _MetadataParseError('count_mismatch')
        if lang == 2 and na not in (0, 1, 2):
            raise _MetadataParseError('count_mismatch')
        if lang in (11, 12, 13, 14, 15) and (na < 2 or na % 2):
            raise _MetadataParseError('count_mismatch')
    if lang not in (1, 2, 11, 12, 13, 14, 15):
        entry['semantic_statuses'].append('unsupported_lang')
    if lang == 2 and lct != 2:
        entry['semantic_statuses'].append('invalid_combination')


def _metadata_law2(file_obj, entry, lct):
    params, _ = _metadata_record(get_tab2_record, file_obj)
    ne = _metadata_int(params[5])
    if ne < 0:
        raise _MetadataParseError('count_mismatch')
    for _ in range(ne):
        items, values = _metadata_record(get_list_record, file_obj)
        lang, nw, nl = (_metadata_int(items[2]), _metadata_int(items[4]),
                        _metadata_int(items[5]))
        if min(nw, nl) < 0 or len(values) != nw:
            raise _MetadataParseError('count_mismatch')
        entry['lang_values'].append(lang)
        if lang == 0 and (nl < 1 or nw != nl):
            raise _MetadataParseError('count_mismatch')
        if lang in (12, 14) and (nl < 2 or nw != 2*nl):
            raise _MetadataParseError('count_mismatch')
        if lang not in (0, 12, 14):
            entry['semantic_statuses'].append('unsupported_lang')
    if lct != 2:
        entry['semantic_statuses'].append('unsupported_frame')


def _metadata_law5(file_obj, entry):
    params, _ = _metadata_record(get_tab2_record, file_obj)
    lidp, ne = _metadata_int(params[2]), _metadata_int(params[5])
    if ne < 0:
        raise _MetadataParseError('count_mismatch')
    entry['lidp'] = lidp
    if lidp not in (0, 1):
        entry['semantic_statuses'].append('invalid_combination')
    for _ in range(ne):
        items, values = _metadata_record(get_list_record, file_obj)
        ltp, nw, nl = (_metadata_int(items[2]), _metadata_int(items[4]),
                       _metadata_int(items[5]))
        if min(nw, nl) < 0 or len(values) != nw:
            raise _MetadataParseError('count_mismatch')
        entry['ltp_values'].append(ltp)
        if ltp == 1 and lidp in (0, 1):
            expected = 4*nl + 3 if lidp == 0 else 3*nl + 3
        elif ltp == 2:
            expected = nl + 1
        elif ltp in (12, 14, 15):
            expected = 2*nl
        else:
            expected = None
        if expected is not None and nw != expected:
            raise _MetadataParseError('count_mismatch')
        if ltp not in (1, 2, 12, 14, 15):
            entry['semantic_statuses'].append('unsupported_ltp')


def _metadata_law6(file_obj, entry):
    apsx, c2, l1, l2, n1, npsx = _metadata_record(
        get_cont_record, file_obj, parse_failure='invalid_law6_control')
    try:
        apsx = float(apsx)
        c2 = float(c2)
        controls = (_metadata_int(l1, 'invalid_law6_control'),
                    _metadata_int(l2, 'invalid_law6_control'),
                    _metadata_int(n1, 'invalid_law6_control'))
        npsx = _metadata_int(npsx, 'invalid_law6_control')
    except _MetadataParseError as exc:
        raise _MetadataParseError('invalid_law6_control') from exc
    except (OverflowError, TypeError, ValueError) as exc:
        raise _MetadataParseError('invalid_law6_control') from exc
    if (not math.isfinite(apsx) or apsx <= 0.0 or c2 != 0.0 or
            controls != (0, 0, 0) or not 3 <= npsx <= 2**31 - 1):
        raise _MetadataParseError('invalid_law6_control')
    entry['apsx'] = apsx
    entry['npsx'] = npsx


def _metadata_law7(file_obj, entry, lct):
    params, _ = _metadata_record(get_tab2_record, file_obj)
    ne = _metadata_int(params[5])
    if ne < 0:
        raise _MetadataParseError('count_mismatch')
    for _ in range(ne):
        inner, _ = _metadata_record(get_tab2_record, file_obj)
        nmu = _metadata_int(inner[5])
        if nmu < 0:
            raise _MetadataParseError('count_mismatch')
        for _ in range(nmu):
            _metadata_record(get_tab1_record, file_obj)
    if lct != 1:
        entry['semantic_statuses'].append('unsupported_frame')


def _metadata_fission_semantics(mt, jp, entries):
    """Apply the ENDF-102 MT=18 multiplicity-product sequence constraints."""
    if jp == 0:
        if any(entry['law'] < 0 for entry in entries):
            for entry in entries:
                entry['semantic_statuses'].append('invalid_combination')
        return
    if mt != 18:
        for entry in entries:
            entry['semantic_statuses'].append('invalid_combination')
        return
    jpp, jpn = divmod(jp, 10)
    kinds = []
    for entry in entries:
        if entry['zap'] == 1:
            kinds.append('neutron')
        elif entry['zap'] == 0:
            kinds.append('photon')
        else:
            entry['semantic_statuses'].append('invalid_combination')
            kinds.append('other')
    if kinds != sorted(kinds, key=lambda value: value == 'photon') or 'other' in kinds:
        for entry in entries:
            entry['semantic_statuses'].append('invalid_combination')
    for kind, digit, negative in (('neutron', jpn, -5), ('photon', jpp, -15)):
        subset = [
            entry for entry, entry_kind in zip(entries, kinds)
            if entry_kind == kind]
        if digit == 0:
            if any(entry['law'] < 0 for entry in subset):
                for entry in subset:
                    entry['semantic_statuses'].append('invalid_combination')
            continue
        if not subset:
            for entry in entries:
                entry['semantic_statuses'].append('invalid_combination')
            continue
        if subset[0]['law'] != 0:
            for entry in subset:
                entry['semantic_statuses'].append('invalid_combination')
            continue
        later = subset[1:]
        valid = bool(later) and (
            all(entry['law'] == negative for entry in later) if digit == 1
            else all(entry['law'] > 0 for entry in later))
        if not valid:
            for entry in subset:
                entry['semantic_statuses'].append('invalid_combination')


def _metadata_validate_published(metadata, expected_mt=None):
    """Validate the closed published object independently of its HDF5 wire."""
    mt = _metadata_int(metadata.mt, 'invalid_header')
    jp = _metadata_int(metadata.jp, 'invalid_header')
    lct = _metadata_int(metadata.lct, 'invalid_header')
    if not 1 <= mt <= 999 or jp not in _METADATA_JP or lct not in (1, 2, 3, 4):
        raise ValueError('invalid evaluated product metadata controls')
    if expected_mt is not None and mt != expected_mt:
        raise ValueError('evaluated product metadata MT does not match reaction')
    if jp != 0 and not metadata.entries:
        raise ValueError('active JP requires an MF=6 product subsection')
    raw_entries = []
    for index, entry in enumerate(metadata.entries):
        controls = (entry.subsection_index, entry.zap, entry.lip, entry.law)
        if any(_metadata_int(value, 'invalid_header') != value
               for value in controls):
            raise ValueError('invalid evaluated product metadata entry control')
        if entry.subsection_index != index:
            raise ValueError('non-contiguous evaluated product metadata entries')
        awp = float(entry.awp)
        if not math.isfinite(awp) or awp < 0.0:
            raise ValueError('invalid evaluated product metadata AWP')
        law = entry.law
        if law not in (-15, -5, 0, 1, 2, 3, 4, 5, 6, 7):
            raise ValueError('invalid evaluated product metadata LAW')
        photon_energy = (mt, entry.zap, law) == (102, 0, 2)
        expected_awp = ('primary_photon_energy_eV' if photon_energy
                        else 'mass_ratio')
        if entry.awp_interpretation != expected_awp:
            raise ValueError('invalid evaluated product metadata AWP interpretation')
        if not _metadata_identity_is_valid(
                entry.zap, entry.lip, entry.identity_status,
                entry.derived_particle):
            raise ValueError('invalid evaluated product metadata identity')
        if entry.semantic_status not in _METADATA_SEMANTIC_STATUS:
            raise ValueError('invalid evaluated product metadata status')
        lang = tuple(_metadata_int(value) for value in entry.lang_values)
        ltp = tuple(_metadata_int(value) for value in entry.ltp_values)
        statuses = []
        if lct == 3 and law != 1:
            statuses.append('invalid_combination')
        if law == 1:
            if len(lang) != 1 or ltp or entry.lidp is not None:
                raise ValueError('invalid LAW=1 metadata selectors')
            if lang[0] not in (1, 2, 11, 12, 13, 14, 15):
                statuses.append('unsupported_lang')
            if lang[0] == 2 and lct != 2:
                statuses.append('invalid_combination')
        elif law == 2:
            if ltp or entry.lidp is not None:
                raise ValueError('invalid LAW=2 metadata selectors')
            if any(value not in (0, 12, 14) for value in lang):
                statuses.append('unsupported_lang')
            if lct != 2:
                statuses.append('unsupported_frame')
        elif law == 5:
            if lang or entry.lidp is None:
                raise ValueError('invalid LAW=5 metadata selectors')
            _metadata_int(entry.lidp)
            if entry.lidp not in (0, 1):
                statuses.append('invalid_combination')
            if any(value not in (1, 2, 12, 14, 15) for value in ltp):
                statuses.append('unsupported_ltp')
        elif lang or ltp or entry.lidp is not None:
            raise ValueError('invalid evaluated product metadata selectors')
        if law == 7 and lct != 1:
            statuses.append('unsupported_frame')
        if law == 6:
            if entry.apsx is None or entry.npsx is None:
                raise ValueError('LAW=6 metadata requires APSX and NPSX')
            if (isinstance(entry.apsx, (bool, np.bool_)) or
                    not isinstance(entry.apsx, Real)):
                raise ValueError('invalid LAW=6 metadata APSX')
            apsx = float(entry.apsx)
            if (not math.isfinite(apsx) or apsx <= 0.0 or
                    isinstance(entry.npsx, (bool, np.bool_)) or
                    not isinstance(entry.npsx, (int, np.integer)) or
                    not 3 <= entry.npsx <= 2**31 - 1):
                raise ValueError('invalid LAW=6 metadata controls')
        elif entry.apsx is not None or entry.npsx is not None:
            raise ValueError('non-LAW=6 metadata cannot contain APSX or NPSX')
        raw_entries.append({
            'zap': entry.zap, 'law': law, 'semantic_statuses': statuses})
    _metadata_fission_semantics(mt, jp, raw_entries)
    expected_statuses = tuple(
        _metadata_status(entry['semantic_statuses']) for entry in raw_entries)
    if expected_statuses != tuple(
            entry.semantic_status for entry in metadata.entries):
        raise ValueError('inconsistent evaluated product metadata semantics')
    if _metadata_status(expected_statuses) != metadata.section_semantic_status:
        raise ValueError('invalid evaluated product metadata aggregate status')


def get_evaluated_product_metadata(ev, mt):
    """Read a transport-isolated MF=6 metadata result for one MT.

    The Evaluation section boundary is trusted because Evaluation has removed
    SEND. All known LAW envelopes are consumed by this reader itself; it never
    calls a sampled-distribution constructor and never touches ``Product``.
    Structural completeness is limited to the supported envelopes and declared
    counts; it is not a full ENDF normalization or conservation audit.

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        Evaluation containing the selected MF=6 section text.
    mt : int
        Requested ENDF reaction number.

    Returns
    -------
    EvaluatedProductMetadataResult
        Closed direct result. Absent and structurally failed results retain the
        requested MT but contain null controls and no entries. Only a
        ``published_structurally_complete`` result exposes attachable
        :attr:`EvaluatedProductMetadataResult.metadata`.

    """
    if (6, mt) not in ev.section:
        return EvaluatedProductMetadataResult(mt=mt, result_status='absent')
    try:
        file_obj = StringIO(ev.section[6, mt])
        head = _metadata_record(get_head_record, file_obj, 'invalid_header')
        jp, lct, nk = (_metadata_int(head[2], 'invalid_header'),
                       _metadata_int(head[3], 'invalid_header'),
                       _metadata_int(head[4], 'invalid_header'))
        head_n2 = _metadata_int(head[5], 'invalid_header')
        if (jp not in _METADATA_JP or lct not in (1, 2, 3, 4) or
                nk < 0 or head_n2 != 0):
            return EvaluatedProductMetadataResult(
                mt=mt, result_status='structural_failure',
                structural_failure='invalid_header')
        if jp != 0 and nk == 0:
            return EvaluatedProductMetadataResult(
                mt=mt, result_status='structural_failure',
                structural_failure='invalid_header')
        raw_entries = []
        for index in range(nk):
            params, _ = _metadata_record(
                get_tab1_record, file_obj, eof_failure='premature_end')
            entry = _metadata_entry(index, params)
            law = entry['law']
            if law not in (-15, -5, 0, 1, 2, 3, 4, 5, 6, 7):
                return EvaluatedProductMetadataResult(
                    mt=mt, result_status='structural_failure',
                    structural_failure='unknown_law')
            if lct == 3 and law != 1:
                entry['semantic_statuses'].append('invalid_combination')
            if law == 1:
                _metadata_law1(file_obj, entry, lct)
            elif law == 2:
                _metadata_law2(file_obj, entry, lct)
            elif law == 5:
                _metadata_law5(file_obj, entry)
            elif law == 6:
                _metadata_law6(file_obj, entry)
            elif law == 7:
                _metadata_law7(file_obj, entry, lct)
            raw_entries.append(entry)
        if file_obj.read().strip():
            return EvaluatedProductMetadataResult(
                mt=mt, result_status='structural_failure',
                structural_failure='trailing_record')
        _metadata_fission_semantics(mt, jp, raw_entries)
        entries = []
        for entry in raw_entries:
            law = entry['law']
            interpretation = ('primary_photon_energy_eV'
                              if (mt, entry['zap'], law) == (102, 0, 2)
                              else 'mass_ratio')
            entries.append(EvaluatedProductMetadataEntry(
                entry['subsection_index'], entry['zap'], entry['awp'],
                interpretation, entry['lip'], law, tuple(entry['lang_values']),
                entry['lidp'], tuple(entry['ltp_values']),
                _metadata_status(entry['semantic_statuses']),
                entry['identity_status'], entry['derived_particle'],
                entry['apsx'], entry['npsx']))
        metadata = EvaluatedProductMetadata(
            mt, jp, lct,
            _metadata_status([entry.semantic_status for entry in entries]),
            tuple(entries))
        _metadata_validate_published(metadata)
        return EvaluatedProductMetadataResult(
            mt=mt, result_status=_METADATA_RESULT, jp=jp, lct=lct,
            section_semantic_status=metadata.section_semantic_status,
            hdf5_group_version=1, entries=metadata.entries)
    except _MetadataParseError as exc:
        return EvaluatedProductMetadataResult(
            mt=mt, result_status='structural_failure',
            structural_failure=exc.reason)
    except (EOFError, IndexError, OSError, ValueError, TypeError):
        return EvaluatedProductMetadataResult(
            mt=mt, result_status='structural_failure',
            structural_failure='unproven_cursor')


def _metadata_text(value):
    """Return a UTF-8 HDF5 attribute as text."""
    if isinstance(value, bytes):
        return value.decode('utf-8')
    if not isinstance(value, str):
        raise ValueError('invalid evaluated product metadata UTF-8 scalar')
    return value


def _metadata_hdf5_scalar(group, name, dtype):
    """Return a scalar with the required HDF5 kind and width."""
    value = group.attrs[name]
    actual = group.attrs.get_id(name).dtype
    expected = np.dtype(dtype)
    if (np.asarray(value).shape != () or actual.kind != expected.kind or
            actual.itemsize != expected.itemsize):
        raise ValueError('invalid evaluated product metadata scalar dtype')
    return value


def _metadata_hdf5_text(group, name):
    """Return a scalar variable-length UTF-8 HDF5 attribute."""
    value = group.attrs[name]
    dtype = group.attrs.get_id(name).dtype
    info = h5py.check_string_dtype(dtype)
    if (np.asarray(value).shape != () or info is None or
            info.encoding != 'utf-8' or info.length is not None):
        raise ValueError('invalid evaluated product metadata string dtype')
    return _metadata_text(value)


def _metadata_to_hdf5(metadata, group, expected_mt):
    """Write the closed metadata wire representation below a reaction group."""
    _metadata_validate_published(metadata, expected_mt)
    metadata_group = group.create_group('evaluated_product_metadata')
    metadata_group.attrs.create('version', np.int32(1), dtype=np.int32)
    metadata_group.attrs.create('schema', _METADATA_SCHEMA,
                                dtype=_METADATA_UTF8)
    metadata_group.attrs.create('source_format', 'ENDF-6',
                                dtype=_METADATA_UTF8)
    metadata_group.attrs.create('result_status', _METADATA_RESULT,
                                dtype=_METADATA_UTF8)
    metadata_group.attrs.create('mf', np.int32(6), dtype=np.int32)
    metadata_group.attrs.create('mt', np.int32(metadata.mt), dtype=np.int32)
    metadata_group.attrs.create('jp', np.int32(metadata.jp), dtype=np.int32)
    metadata_group.attrs.create('lct', np.int32(metadata.lct), dtype=np.int32)
    metadata_group.attrs.create(
        'section_semantic_status', metadata.section_semantic_status,
        dtype=_METADATA_UTF8)
    for index, entry in enumerate(metadata.entries):
        entry_group = metadata_group.create_group(f'entry_{index}')
        entry_group.attrs.create(
            'subsection_index', np.int32(entry.subsection_index),
            dtype=np.int32)
        entry_group.attrs.create('zap', np.int32(entry.zap), dtype=np.int32)
        entry_group.attrs.create('awp', np.float64(entry.awp),
                                 dtype=np.float64)
        entry_group.attrs.create(
            'awp_interpretation', entry.awp_interpretation,
            dtype=_METADATA_UTF8)
        entry_group.attrs.create('lip', np.int32(entry.lip), dtype=np.int32)
        entry_group.attrs.create('law', np.int32(entry.law), dtype=np.int32)
        entry_group.attrs.create(
            'lidp', np.int32(-1 if entry.lidp is None else entry.lidp),
            dtype=np.int32)
        entry_group.attrs.create(
            'semantic_status', entry.semantic_status, dtype=_METADATA_UTF8)
        entry_group.attrs.create(
            'identity_status', entry.identity_status, dtype=_METADATA_UTF8)
        entry_group.attrs.create(
            'derived_particle', entry.derived_particle or '',
            dtype=_METADATA_UTF8)
        if entry.law == 6:
            entry_group.attrs.create('apsx', np.float64(entry.apsx),
                                     dtype=np.float64)
            entry_group.attrs.create('npsx', np.int32(entry.npsx),
                                     dtype=np.int32)
        entry_group.create_dataset(
            'lang_values', data=np.asarray(entry.lang_values, dtype=np.int32),
            dtype=np.int32)
        entry_group.create_dataset(
            'ltp_values', data=np.asarray(entry.ltp_values, dtype=np.int32),
            dtype=np.int32)


def _metadata_from_hdf5(group, mt):
    """Read and validate the closed metadata wire representation."""
    required = {'version', 'schema', 'source_format', 'result_status', 'mf',
                'mt', 'jp', 'lct', 'section_semantic_status'}
    if set(group.attrs) != required:
        raise ValueError('invalid evaluated product metadata group attributes')
    if _metadata_int(_metadata_hdf5_scalar(
            group, 'version', np.int32)) != 1:
        raise ValueError('invalid evaluated product metadata group version')
    if (_metadata_hdf5_text(group, 'schema') != _METADATA_SCHEMA or
            _metadata_hdf5_text(group, 'source_format') != 'ENDF-6' or
            _metadata_hdf5_text(group, 'result_status') != _METADATA_RESULT or
            _metadata_int(_metadata_hdf5_scalar(group, 'mf', np.int32)) != 6 or
            _metadata_int(_metadata_hdf5_scalar(
                group, 'mt', np.int32)) != mt):
        raise ValueError('invalid evaluated product metadata group identity')
    jp = _metadata_int(_metadata_hdf5_scalar(group, 'jp', np.int32))
    lct = _metadata_int(_metadata_hdf5_scalar(group, 'lct', np.int32))
    if jp not in _METADATA_JP or lct not in (1, 2, 3, 4):
        raise ValueError('invalid evaluated product metadata controls')
    section_status = _metadata_hdf5_text(
        group, 'section_semantic_status')
    if section_status not in _METADATA_SEMANTIC_STATUS:
        raise ValueError('invalid evaluated product metadata status')
    names = sorted(
        group,
        key=lambda name: int(name[6:])
        if name.startswith('entry_') and name[6:].isdigit() else -1)
    if names != [f'entry_{i}' for i in range(len(names))]:
        raise ValueError('non-contiguous evaluated product metadata entries')
    entries = []
    entry_attrs = {'subsection_index', 'zap', 'awp', 'awp_interpretation', 'lip',
                   'law', 'lidp', 'semantic_status', 'identity_status',
                   'derived_particle'}
    for index, name in enumerate(names):
        entry_group = group[name]
        if (not isinstance(entry_group, h5py.Group) or
                not entry_attrs <= set(entry_group.attrs) or
                set(entry_group) != {'lang_values', 'ltp_values'}):
            raise ValueError('invalid evaluated product metadata entry')
        subsection_index = _metadata_int(_metadata_hdf5_scalar(
            entry_group, 'subsection_index', np.int32))
        zap = _metadata_int(_metadata_hdf5_scalar(
            entry_group, 'zap', np.int32))
        lip = _metadata_int(_metadata_hdf5_scalar(
            entry_group, 'lip', np.int32))
        law = _metadata_int(_metadata_hdf5_scalar(
            entry_group, 'law', np.int32))
        expected_attrs = entry_attrs | ({'apsx', 'npsx'} if law == 6 else set())
        if set(entry_group.attrs) != expected_attrs:
            raise ValueError('invalid evaluated product metadata entry')
        lidp = _metadata_int(_metadata_hdf5_scalar(
            entry_group, 'lidp', np.int32))
        awp = float(_metadata_hdf5_scalar(
            entry_group, 'awp', np.float64))
        if subsection_index != index or not math.isfinite(awp) or awp < 0.0:
            raise ValueError('invalid evaluated product metadata scalar')
        if law not in (-15, -5, 0, 1, 2, 3, 4, 5, 6, 7):
            raise ValueError('invalid evaluated product metadata LAW')
        if law == 6:
            apsx = float(_metadata_hdf5_scalar(
                entry_group, 'apsx', np.float64))
            npsx = _metadata_int(_metadata_hdf5_scalar(
                entry_group, 'npsx', np.int32))
            if (not math.isfinite(apsx) or apsx <= 0.0 or
                    not 3 <= npsx <= 2**31 - 1):
                raise ValueError('invalid LAW=6 metadata controls')
        else:
            apsx = None
            npsx = None
        interpretation = _metadata_hdf5_text(
            entry_group, 'awp_interpretation')
        photon_energy = (mt, zap, law) == (102, 0, 2)
        expected_interpretation = (
            'primary_photon_energy_eV' if photon_energy else 'mass_ratio')
        if interpretation != expected_interpretation:
            raise ValueError('invalid evaluated product metadata AWP interpretation')
        lang_dset = entry_group['lang_values']
        ltp_dset = entry_group['ltp_values']
        if (not isinstance(lang_dset, h5py.Dataset) or
                not isinstance(ltp_dset, h5py.Dataset) or
                lang_dset.ndim != 1 or ltp_dset.ndim != 1 or
                lang_dset.dtype.kind != np.dtype(np.int32).kind or
                ltp_dset.dtype.kind != np.dtype(np.int32).kind or
                lang_dset.dtype.itemsize != np.dtype(np.int32).itemsize or
                ltp_dset.dtype.itemsize != np.dtype(np.int32).itemsize):
            raise ValueError('invalid evaluated product metadata selector dtype')
        lang = tuple(int(value) for value in lang_dset[()])
        ltp = tuple(int(value) for value in ltp_dset[()])
        status = _metadata_hdf5_text(entry_group, 'semantic_status')
        identity = _metadata_hdf5_text(entry_group, 'identity_status')
        derived = _metadata_hdf5_text(
            entry_group, 'derived_particle') or None
        identities = ('raw_only', 'special_particle',
                      'ground_state_nuclide')
        if status not in _METADATA_SEMANTIC_STATUS or identity not in identities:
            raise ValueError('invalid evaluated product metadata entry status')
        if law == 1 and (len(lang) != 1 or ltp or lidp != -1):
            raise ValueError('invalid LAW=1 metadata selectors')
        if law == 2 and (ltp or lidp != -1):
            raise ValueError('invalid LAW=2 metadata selectors')
        if law == 5 and (lang or lidp < 0):
            raise ValueError('invalid LAW=5 metadata selectors')
        if law not in (1, 2, 5) and (lang or ltp or lidp != -1):
            raise ValueError('invalid metadata selectors')
        if not _metadata_identity_is_valid(zap, lip, identity, derived):
            raise ValueError('invalid evaluated product metadata identity')
        entries.append(EvaluatedProductMetadataEntry(
            subsection_index, zap, awp, interpretation, lip, law, lang,
            None if lidp == -1 else lidp, ltp, status, identity, derived,
            apsx, npsx))
    metadata = EvaluatedProductMetadata(
        mt, jp, lct, section_status, tuple(entries))
    _metadata_validate_published(metadata)
    return metadata


def _get_products(ev, mt):
    """Generate products from MF=6 in an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        ENDF evaluation to read from
    mt : int
        The MT value of the reaction to get products for

    Raises
    ------
    IOError
        When the Kalbach-Mann systematics is used, but the product
        is not defined in the 'center-of-mass' system. The breakup logic
        is not implemented which can lead to this error being raised while
        the definition of the product is correct.

    Returns
    -------
    products : list of openmc.data.Product
        Products of the reaction

    """
    file_obj = StringIO(ev.section[6, mt])

    # Read HEAD record
    items = get_head_record(file_obj)
    reference_frame = {1: 'laboratory', 2: 'center-of-mass',
                       3: 'light-heavy', 4: 'breakup'}[items[3]]
    n_products = items[4]

    products = []
    for i in range(n_products):
        # Get yield for this product
        params, yield_ = get_tab1_record(file_obj)

        za = int(params[0])
        awr = params[1]
        law = params[3]

        if za == 0:
            p = Product('photon')
        elif za == 1:
            p = Product('neutron')
        elif za == 1000:
            p = Product('electron')
        else:
            Z, A = divmod(za, 1000)
            p = Product(f'{ATOMIC_SYMBOL[Z]}{A}')

        p.yield_ = yield_

        """
        # Set reference frame
        if reference_frame == 'laboratory':
            p.center_of_mass = False
        elif reference_frame == 'center-of-mass':
            p.center_of_mass = True
        elif reference_frame == 'light-heavy':
            p.center_of_mass = (awr <= 4.0)
        """

        if law == 0:
            # No distribution given
            pass
        if law == 1:
            # Continuum energy-angle distribution

            # Peak ahead to determine type of distribution
            position = file_obj.tell()
            params = get_cont_record(file_obj)
            file_obj.seek(position)

            lang = params[2]
            if lang == 1:
                p.distribution = [CorrelatedAngleEnergy.from_endf(file_obj)]
            elif lang == 2:
                # Products need to be described in the center-of-mass system
                product_center_of_mass = False
                if reference_frame == 'center-of-mass':
                    product_center_of_mass = True
                elif reference_frame == 'light-heavy':
                    product_center_of_mass = (awr <= 4.0)
                # TODO: 'breakup' logic not implemented

                if product_center_of_mass is False:
                    raise IOError(
                        "Kalbach-Mann representation must be defined in the "
                        "'center-of-mass' system"
                    )

                zat = ev.target["atomic_number"] * 1000 + ev.target["mass_number"]
                projectile_mass = ev.projectile["mass"]
                p.distribution = [KalbachMann.from_endf(file_obj,
                                                        za,
                                                        zat,
                                                        projectile_mass)]

        elif law == 2:
            # Discrete two-body scattering
            params, tab2 = get_tab2_record(file_obj)
            ne = params[5]
            energy = np.zeros(ne)
            mu = []
            for i in range(ne):
                items, values = get_list_record(file_obj)
                energy[i] = items[1]
                lang = items[2]
                if lang == 0:
                    mu.append(Legendre(values))
                elif lang == 12:
                    mu.append(Tabular(values[::2], values[1::2]))
                elif lang == 14:
                    mu.append(Tabular(values[::2], values[1::2],
                                      'log-linear'))

            angle_dist = AngleDistribution(energy, mu)
            dist = UncorrelatedAngleEnergy(angle_dist)
            p.distribution = [dist]
            # TODO: Add level-inelastic info?

        elif law == 3:
            # Isotropic discrete emission
            p.distribution = [UncorrelatedAngleEnergy()]
            # TODO: Add level-inelastic info?

        elif law == 4:
            # Discrete two-body recoil
            pass

        elif law == 5:
            # Charged particle elastic scattering
            pass

        elif law == 6:
            # N-body phase-space distribution
            p.distribution = [NBodyPhaseSpace.from_endf(file_obj)]

        elif law == 7:
            # Laboratory energy-angle distribution
            p.distribution = [LaboratoryAngleEnergy.from_endf(file_obj)]

        products.append(p)

    return products


def _get_fission_products_ace(ace):
    """Generate fission products from an ACE table

    Parameters
    ----------
    ace : openmc.data.ace.Table
        ACE table to read from

    Returns
    -------
    products : list of openmc.data.Product
        Prompt and delayed fission neutrons
    derived_products : list of openmc.data.Product
        "Total" fission neutron

    """
    # No NU block
    if ace.jxs[2] == 0:
        return None, None

    products = []
    derived_products = []

    # Either prompt nu or total nu is given
    if ace.xss[ace.jxs[2]] > 0:
        whichnu = 'prompt' if ace.jxs[24] > 0 else 'total'

        neutron = Product('neutron')
        neutron.emission_mode = whichnu

        idx = ace.jxs[2]
        LNU = int(ace.xss[idx])
        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
            neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        products.append(neutron)

    # Both prompt nu and total nu
    elif ace.xss[ace.jxs[2]] < 0:
        # Read prompt neutron yield
        prompt_neutron = Product('neutron')
        prompt_neutron.emission_mode = 'prompt'

        idx = ace.jxs[2] + 1
        LNU = int(ace.xss[idx])
        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
            prompt_neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            prompt_neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        # Read total neutron yield
        total_neutron = Product('neutron')
        total_neutron.emission_mode = 'total'

        idx = ace.jxs[2] + int(abs(ace.xss[ace.jxs[2]])) + 1
        LNU = int(ace.xss[idx])

        if LNU == 1:
            # Polynomial function form of nu
            NC = int(ace.xss[idx+1])
            coefficients = ace.xss[idx+2 : idx+2+NC].copy()
            for i in range(coefficients.size):
                coefficients[i] *= EV_PER_MEV**(-i)
            total_neutron.yield_ = Polynomial(coefficients)
        elif LNU == 2:
            # Tabular data form of nu
            total_neutron.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        products.append(prompt_neutron)
        derived_products.append(total_neutron)

    # Check for delayed nu data
    if ace.jxs[24] > 0:
        yield_delayed = Tabulated1D.from_ace(ace, ace.jxs[24] + 1)

        # Delayed neutron precursor distribution
        idx = ace.jxs[25]
        n_group = ace.nxs[8]
        total_group_probability = 0.
        for group in range(n_group):
            delayed_neutron = Product('neutron')
            delayed_neutron.emission_mode = 'delayed'

            # Convert units of inverse shakes to inverse seconds
            delayed_neutron.decay_rate = ace.xss[idx] * 1.e8

            group_probability = Tabulated1D.from_ace(ace, idx + 1)
            if np.all(group_probability.y == group_probability.y[0]):
                delayed_neutron.yield_ = deepcopy(yield_delayed)
                delayed_neutron.yield_.y *= group_probability.y[0]
                total_group_probability += group_probability.y[0]
            else:
                # Get union energy grid and ensure energies are within
                # interpolable range of both functions
                max_energy = min(yield_delayed.x[-1], group_probability.x[-1])
                energy = np.union1d(yield_delayed.x, group_probability.x)
                energy = energy[energy <= max_energy]

                # Calculate group yield
                group_yield = yield_delayed(energy) * group_probability(energy)
                delayed_neutron.yield_ = Tabulated1D(energy, group_yield)

            # Advance position
            nr = int(ace.xss[idx + 1])
            ne = int(ace.xss[idx + 2 + 2*nr])
            idx += 3 + 2*nr + 2*ne

            # Energy distribution for delayed fission neutrons
            location_start = int(ace.xss[ace.jxs[26] + group])
            delayed_neutron.distribution.append(
                AngleEnergy.from_ace(ace, ace.jxs[27], location_start))

            products.append(delayed_neutron)

        # Renormalize delayed neutron yields to reflect fact that in ACE
        # file, the sum of the group probabilities is not exactly one
        for product in products[1:]:
            if total_group_probability > 0.:
                product.yield_.y /= total_group_probability

    return products, derived_products


def _get_fission_products_endf(ev):
    """Generate fission products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation

    Returns
    -------
    products : list of openmc.data.Product
        Prompt and delayed fission neutrons
    derived_products : list of openmc.data.Product
        "Total" fission neutron

    """
    products = []
    derived_products = []

    if (1, 456) in ev.section:
        prompt_neutron = Product('neutron')
        prompt_neutron.emission_mode = 'prompt'

        # Prompt nu values
        file_obj = StringIO(ev.section[1, 456])
        lnu = get_head_record(file_obj)[3]
        if lnu == 1:
            # Polynomial representation
            items, coefficients = get_list_record(file_obj)
            prompt_neutron.yield_ = Polynomial(coefficients)
        elif lnu == 2:
            # Tabulated representation
            params, prompt_neutron.yield_ = get_tab1_record(file_obj)

        products.append(prompt_neutron)

    if (1, 452) in ev.section:
        total_neutron = Product('neutron')
        total_neutron.emission_mode = 'total'

        # Total nu values
        file_obj = StringIO(ev.section[1, 452])
        lnu = get_head_record(file_obj)[3]
        if lnu == 1:
            # Polynomial representation
            items, coefficients = get_list_record(file_obj)
            total_neutron.yield_ = Polynomial(coefficients)
        elif lnu == 2:
            # Tabulated representation
            params, total_neutron.yield_ = get_tab1_record(file_obj)

        if (1, 456) in ev.section:
            derived_products.append(total_neutron)
        else:
            products.append(total_neutron)

    if (1, 455) in ev.section:
        file_obj = StringIO(ev.section[1, 455])

        # Determine representation of delayed nu data
        items = get_head_record(file_obj)
        ldg = items[2]
        lnu = items[3]

        if ldg == 0:
            # Delayed-group constants energy independent
            items, decay_constants = get_list_record(file_obj)
            for constant in decay_constants:
                delayed_neutron = Product('neutron')
                delayed_neutron.emission_mode = 'delayed'
                delayed_neutron.decay_rate = constant
                products.append(delayed_neutron)
        elif ldg == 1:
            # Delayed-group constants energy dependent
            raise NotImplementedError('Delayed neutron with energy-dependent '
                                      'group constants.')

        # In MF=1, MT=455, the delayed-group abundances are actually not
        # specified if the group constants are energy-independent. In this case,
        # the abundances must be inferred from MF=5, MT=455 where multiple
        # energy distributions are given.
        if lnu == 1:
            # Nu represented as polynomial
            items, coefficients = get_list_record(file_obj)
            yield_ = Polynomial(coefficients)
            for neutron in products[-6:]:
                neutron.yield_ = deepcopy(yield_)
        elif lnu == 2:
            # Nu represented by tabulation
            params, yield_ = get_tab1_record(file_obj)
            for neutron in products[-6:]:
                neutron.yield_ = deepcopy(yield_)

        if (5, 455) in ev.section:
            file_obj = StringIO(ev.section[5, 455])
            items = get_head_record(file_obj)
            nk = items[4]
            if nk > 1 and len(decay_constants) == 1:
                # If only one precursor group is listed in MF=1, MT=455, use the
                # energy spectra from MF=5 to split them into different groups
                for _ in range(nk - 1):
                    products.append(deepcopy(products[1]))
            elif nk != len(decay_constants):
                raise ValueError(
                    'Number of delayed neutron fission spectra ({}) does not '
                    'match number of delayed neutron precursors ({}).'.format(
                        nk, len(decay_constants)))
            for i in range(nk):
                params, applicability = get_tab1_record(file_obj)
                dist = UncorrelatedAngleEnergy()
                dist.energy = EnergyDistribution.from_endf(file_obj, params)

                delayed_neutron = products[1 + i]
                yield_ = delayed_neutron.yield_

                # Here we handle the fact that the delayed neutron yield is the
                # product of the total delayed neutron yield and the
                # "applicability" of the energy distribution law in file 5.
                if isinstance(yield_, Tabulated1D):
                    if np.all(applicability.y == applicability.y[0]):
                        yield_.y *= applicability.y[0]
                    else:
                        # Get union energy grid and ensure energies are within
                        # interpolable range of both functions
                        max_energy = min(yield_.x[-1], applicability.x[-1])
                        energy = np.union1d(yield_.x, applicability.x)
                        energy = energy[energy <= max_energy]

                        # Calculate group yield
                        group_yield = yield_(energy) * applicability(energy)
                        delayed_neutron.yield_ = Tabulated1D(energy, group_yield)
                elif isinstance(yield_, Polynomial):
                    if len(yield_) == 1:
                        delayed_neutron.yield_ = deepcopy(applicability)
                        delayed_neutron.yield_.y *= yield_.coef[0]
                    else:
                        if np.all(applicability.y == applicability.y[0]):
                            yield_.coef[0] *= applicability.y[0]
                        else:
                            raise NotImplementedError(
                                'Total delayed neutron yield and delayed group '
                                'probability are both energy-dependent.')

                delayed_neutron.distribution.append(dist)

    return products, derived_products


def _get_activation_products(ev, rx):
    """Generate activation products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        The ENDF evaluation
    rx : openmc.data.Reaction
        Reaction which generates activation products

    Returns
    -------
    products : list of openmc.data.Product
        Activation products

    """
    file_obj = StringIO(ev.section[8, rx.mt])

    # Determine total number of states and whether decay chain is given in a
    # decay sublibrary
    items = get_head_record(file_obj)
    n_states = items[4]
    decay_sublib = (items[5] == 1)

    # Determine if file 9/10 are present
    present = {9: False, 10: False}
    for _ in range(n_states):
        if decay_sublib:
            items = get_cont_record(file_obj)
        else:
            items, values = get_list_record(file_obj)
        lmf = items[2]
        if lmf == 9:
            present[9] = True
        elif lmf == 10:
            present[10] = True

    products = []

    for mf in (9, 10):
        if not present[mf]:
            continue

        file_obj = StringIO(ev.section[mf, rx.mt])
        items = get_head_record(file_obj)
        n_states = items[4]
        for i in range(n_states):
            # Determine what the product is
            items, xs = get_tab1_record(file_obj)
            Z, A = divmod(items[2], 1000)
            excited_state = items[3]

            # Get GNDS name for product
            symbol = ATOMIC_SYMBOL[Z]
            if excited_state > 0:
                name = f'{symbol}{A}_e{excited_state}'
            else:
                name = f'{symbol}{A}'

            p = Product(name)
            if mf == 9:
                p.yield_ = xs
            else:
                # Re-interpolate production cross section and neutron cross
                # section to union energy grid
                energy = np.union1d(xs.x, rx.xs['0K'].x)
                prod_xs = xs(energy)
                neutron_xs = rx.xs['0K'](energy)
                idx = np.where(neutron_xs > 0)

                # Calculate yield as ratio
                yield_ = np.zeros_like(energy)
                yield_[idx] = prod_xs[idx] / neutron_xs[idx]
                p.yield_ = Tabulated1D(energy, yield_)

            # Check if product already exists from MF=6 and if it does, just
            # overwrite the existing yield.
            for product in rx.products:
                if name == product.particle:
                    product.yield_ = p.yield_
                    break
            else:
                products.append(p)

    return products


def _get_photon_products_ace(ace, rx):
    """Generate photon products from an ACE table

    Parameters
    ----------
    ace : openmc.data.ace.Table
        ACE table to read from
    rx : openmc.data.Reaction
        Reaction that generates photons

    Returns
    -------
    photons : list of openmc.Products
        Photons produced from reaction with given MT

    """
    n_photon_reactions = ace.nxs[6]
    photon_mts = ace.xss[ace.jxs[13]:ace.jxs[13] +
                         n_photon_reactions].astype(int)

    photons = []
    for i in range(n_photon_reactions):
        # Determine corresponding reaction
        neutron_mt = photon_mts[i] // 1000

        if neutron_mt != rx.mt:
            continue

        # Create photon product and assign to reactions
        photon = Product('photon')

        # ==================================================================
        # Photon yield / production cross section

        loca = int(ace.xss[ace.jxs[14] + i])
        idx = ace.jxs[15] + loca - 1
        mftype = int(ace.xss[idx])
        idx += 1

        if mftype in (12, 16):
            # Yield data taken from ENDF File 12 or 6
            mtmult = int(ace.xss[idx])
            assert mtmult == neutron_mt

            # Read photon yield as function of energy
            photon.yield_ = Tabulated1D.from_ace(ace, idx + 1)

        elif mftype == 13:
            # Cross section data from ENDF File 13

            # Energy grid index at which data starts
            threshold_idx = int(ace.xss[idx]) - 1
            n_energy = int(ace.xss[idx + 1])
            energy = ace.xss[ace.jxs[1] + threshold_idx:
                             ace.jxs[1] + threshold_idx + n_energy]*EV_PER_MEV

            # Get photon production cross section
            photon_prod_xs = ace.xss[idx + 2:idx + 2 + n_energy]
            neutron_xs = list(rx.xs.values())[0](energy)
            idx = np.where(neutron_xs > 0.)

            # Calculate photon yield
            yield_ = np.zeros_like(photon_prod_xs)
            yield_[idx] = photon_prod_xs[idx] / neutron_xs[idx]
            photon.yield_ = Tabulated1D(energy, yield_)

        else:
            raise ValueError(f"MFTYPE must be 12, 13, 16. Got {mftype}")

        # ==================================================================
        # Photon energy distribution

        location_start = int(ace.xss[ace.jxs[18] + i])
        distribution = AngleEnergy.from_ace(ace, ace.jxs[19], location_start)
        assert isinstance(distribution, UncorrelatedAngleEnergy)

        # ==================================================================
        # Photon angular distribution
        loc = int(ace.xss[ace.jxs[16] + i])

        if loc == 0:
            # No angular distribution data are given for this reaction,
            # isotropic scattering is asssumed in LAB
            energy = np.array([photon.yield_.x[0], photon.yield_.x[-1]])
            mu_isotropic = Uniform(-1., 1.)
            distribution.angle = AngleDistribution(
                energy, [mu_isotropic, mu_isotropic])
        else:
            distribution.angle = AngleDistribution.from_ace(ace, ace.jxs[17], loc)

        # Add to list of distributions
        photon.distribution.append(distribution)
        photons.append(photon)

    return photons


def _get_photon_products_endf(ev, rx):
    """Generate photon products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation
        ENDF evaluation to read from
    rx : openmc.data.Reaction
        Reaction that generates photons

    Returns
    -------
    products : list of openmc.Products
        Photons produced from reaction with given MT

    """
    products = []

    if (12, rx.mt) in ev.section:
        file_obj = StringIO(ev.section[12, rx.mt])

        items = get_head_record(file_obj)
        option = items[2]

        if option == 1:
            # Multiplicities given
            n_discrete_photon = items[4]
            if n_discrete_photon > 1:
                items, total_yield = get_tab1_record(file_obj)
            for k in range(n_discrete_photon):
                photon = Product('photon')

                # Get photon yield
                items, photon.yield_ = get_tab1_record(file_obj)

                # Get photon energy distribution
                law = items[3]
                dist = UncorrelatedAngleEnergy()
                if law == 1:
                    # TODO: Get file 15 distribution
                    pass
                elif law == 2:
                    energy = items[0]
                    primary_flag = items[2]
                    dist.energy = DiscretePhoton(primary_flag, energy,
                                                 ev.target['mass'])

                photon.distribution.append(dist)
                products.append(photon)

        elif option == 2:
            # Transition probability arrays given
            ppyield = {}
            ppyield['type'] = 'transition'
            ppyield['transition'] = transition = {}

            # Determine whether simple (LG=1) or complex (LG=2) transitions
            lg = items[3]

            # Get transition data
            items, values = get_list_record(file_obj)
            transition['energy_start'] = items[0]
            transition['energies'] = np.array(values[::lg + 1])
            transition['direct_probability'] = np.array(values[1::lg + 1])
            if lg == 2:
                # Complex case
                transition['conditional_probability'] = np.array(
                    values[2::lg + 1])

    elif (13, rx.mt) in ev.section:
        file_obj = StringIO(ev.section[13, rx.mt])

        # Determine option
        items = get_head_record(file_obj)
        n_discrete_photon = items[4]
        if n_discrete_photon > 1:
            items, total_xs = get_tab1_record(file_obj)
        for k in range(n_discrete_photon):
            photon = Product('photon')
            items, xs = get_tab1_record(file_obj)

            # Re-interpolate photon production cross section and neutron cross
            # section to union energy grid
            energy = np.union1d(xs.x, rx.xs['0K'].x)
            photon_prod_xs = xs(energy)
            neutron_xs = rx.xs['0K'](energy)
            idx = np.where(neutron_xs > 0)

            # Calculate yield as ratio
            yield_ = np.zeros_like(energy)
            yield_[idx] = photon_prod_xs[idx] / neutron_xs[idx]
            photon.yield_ = Tabulated1D(energy, yield_)

            # Get photon energy distribution
            law = items[3]
            dist = UncorrelatedAngleEnergy()
            if law == 1:
                # TODO: Get file 15 distribution
                pass
            elif law == 2:
                energy = items[1]
                primary_flag = items[2]
                dist.energy = DiscretePhoton(primary_flag, energy,
                                             ev.target['mass'])

            photon.distribution.append(dist)
            products.append(photon)

    return products


class Reaction(EqualityMixin):
    """A nuclear reaction

    A Reaction object represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution.

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction.

    Attributes
    ----------
    center_of_mass : bool
        Indicates whether scattering kinematics should be performed in the
        center-of-mass or laboratory reference frame.
        grid above the threshold value in barns.
    redundant : bool
        Indicates whether or not this is a redundant reaction
    mt : int
        The ENDF MT number for this reaction.
    q_value : float
        The Q-value of this reaction in eV.
    xs : dict of str to openmc.data.Function1D
        Microscopic cross section for this reaction as a function of incident
        energy; these cross sections are provided in a dictionary where the key
        is the temperature of the cross section set.
    products : Iterable of openmc.data.Product
        Reaction products
    derived_products : Iterable of openmc.data.Product
        Derived reaction products. Used for 'total' fission neutron data when
        prompt/delayed data also exists.

    """

    def __init__(self, mt):
        self._center_of_mass = True
        self._redundant = False
        self._q_value = 0.
        self._xs = {}
        self._products = []
        self._derived_products = []
        self.evaluated_product_metadata = None

        self.mt = mt

    def __repr__(self):
        if self.mt in REACTION_NAME:
            return f"<Reaction: MT={self.mt} {REACTION_NAME[self.mt]}>"
        else:
            return f"<Reaction: MT={self.mt}>"

    @property
    def center_of_mass(self):
        return self._center_of_mass

    @center_of_mass.setter
    def center_of_mass(self, center_of_mass):
        cv.check_type('center of mass', center_of_mass, (bool, np.bool_))
        self._center_of_mass = center_of_mass

    @property
    def redundant(self):
        return self._redundant

    @redundant.setter
    def redundant(self, redundant):
        cv.check_type('redundant', redundant, (bool, np.bool_))
        self._redundant = redundant

    @property
    def q_value(self):
        return self._q_value

    @q_value.setter
    def q_value(self, q_value):
        cv.check_type('Q value', q_value, Real)
        self._q_value = q_value

    @property
    def products(self):
        return self._products

    @products.setter
    def products(self, products):
        cv.check_type('reaction products', products, Iterable, Product)
        self._products = products

    @property
    def derived_products(self):
        return self._derived_products

    @derived_products.setter
    def derived_products(self, derived_products):
        cv.check_type('reaction derived products', derived_products,
                      Iterable, Product)
        self._derived_products = derived_products

    @property
    def xs(self):
        return self._xs

    @xs.setter
    def xs(self, xs):
        cv.check_type('reaction cross section dictionary', xs, MutableMapping)
        for key, value in xs.items():
            cv.check_type('reaction cross section temperature', key, str)
            cv.check_type('reaction cross section', value, Callable)
        self._xs = xs

    def to_hdf5(self, group):
        """Write reaction to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['mt'] = self.mt
        if self.mt in REACTION_NAME:
            group.attrs['label'] = np.bytes_(REACTION_NAME[self.mt])
        else:
            group.attrs['label'] = np.bytes_(self.mt)
        group.attrs['Q_value'] = self.q_value
        group.attrs['center_of_mass'] = 1 if self.center_of_mass else 0
        group.attrs['redundant'] = 1 if self.redundant else 0
        for T in self.xs:
            Tgroup = group.create_group(T)
            if self.xs[T] is not None:
                dset = Tgroup.create_dataset('xs', data=self.xs[T].y)
                threshold_idx = getattr(self.xs[T], '_threshold_idx', 0)
                dset.attrs['threshold_idx'] = threshold_idx
        for i, p in enumerate(self.products):
            pgroup = group.create_group(f'product_{i}')
            p.to_hdf5(pgroup)
        if self.evaluated_product_metadata is not None:
            _metadata_to_hdf5(
                self.evaluated_product_metadata, group, self.mt)

    @classmethod
    def from_hdf5(cls, group, energy):
        """Generate reaction from an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from
        energy : dict
            Dictionary whose keys are temperatures (e.g., '300K') and values are
            arrays of energies at which cross sections are tabulated at.

        Returns
        -------
        openmc.data.Reaction
            Reaction data

        """

        mt = group.attrs['mt']
        rx = cls(mt)
        rx.q_value = group.attrs['Q_value']
        rx.center_of_mass = bool(group.attrs['center_of_mass'])
        rx.redundant = bool(group.attrs.get('redundant', False))

        # Read cross section at each temperature
        for T, Tgroup in group.items():
            if T.endswith('K'):
                if 'xs' in Tgroup:
                    # Make sure temperature has associated energy grid
                    if T not in energy:
                        raise ValueError(
                            'Could not create reaction cross section for MT={} '
                            'at T={} because no corresponding energy grid '
                            'exists.'.format(mt, T))
                    xs = Tgroup['xs'][()]
                    threshold_idx = Tgroup['xs'].attrs['threshold_idx']
                    tabulated_xs = Tabulated1D(energy[T][threshold_idx:], xs)
                    tabulated_xs._threshold_idx = threshold_idx
                    rx.xs[T] = tabulated_xs

        # Determine number of products
        n_product = 0
        for name in group:
            if name.startswith('product_'):
                n_product += 1

        # Read reaction products
        for i in range(n_product):
            pgroup = group[f'product_{i}']
            rx.products.append(Product.from_hdf5(pgroup))

        if 'evaluated_product_metadata' in group:
            rx.evaluated_product_metadata = _metadata_from_hdf5(
                group['evaluated_product_metadata'], mt)

        return rx

    @classmethod
    def from_ace(cls, ace, i_reaction):
        # Get nuclide energy grid
        n_grid = ace.nxs[3]
        grid = ace.xss[ace.jxs[1]:ace.jxs[1] + n_grid]*EV_PER_MEV

        # Convert data temperature to a "300.0K" number for indexing
        # temperature data
        strT = str(int(round(ace.temperature*EV_PER_MEV / K_BOLTZMANN))) + "K"

        if i_reaction > 0:
            mt = int(ace.xss[ace.jxs[3] + i_reaction - 1])
            rx = cls(mt)

            # Get Q-value of reaction
            rx.q_value = ace.xss[ace.jxs[4] + i_reaction - 1]*EV_PER_MEV

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(ace.xss[ace.jxs[6] + i_reaction - 1])

            # Determine starting index on energy grid
            threshold_idx = int(ace.xss[ace.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(ace.xss[ace.jxs[7] + loc])
            energy = grid[threshold_idx:threshold_idx + n_energy]

            # Read reaction cross section
            xs = ace.xss[ace.jxs[7] + loc + 1:ace.jxs[7] + loc + 1 + n_energy]

            # For damage energy production, convert to eV
            if mt == 444:
                xs *= EV_PER_MEV

            # Fix negatives -- known issue for Y89 in JEFF 3.2
            if np.any(xs < 0.0):
                warn("Negative cross sections found for MT={} in {}. Setting "
                     "to zero.".format(rx.mt, ace.name))
                xs[xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(energy, xs)
            tabulated_xs._threshold_idx = threshold_idx
            rx.xs[strT] = tabulated_xs

            # ==================================================================
            # YIELD AND ANGLE-ENERGY DISTRIBUTION

            # Determine multiplicity
            ty = int(ace.xss[ace.jxs[5] + i_reaction - 1])
            rx.center_of_mass = (ty < 0)
            if i_reaction < ace.nxs[5] + 1:
                if ty != 19:
                    if abs(ty) > 100:
                        # Energy-dependent neutron yield
                        idx = ace.jxs[11] + abs(ty) - 101
                        yield_ = Tabulated1D.from_ace(ace, idx)
                    else:
                        # 0-order polynomial i.e. a constant
                        yield_ = Polynomial((abs(ty),))

                    neutron = Product('neutron')
                    neutron.yield_ = yield_
                    rx.products.append(neutron)
                else:
                    assert mt in FISSION_MTS
                    rx.products, rx.derived_products = _get_fission_products_ace(ace)

                    for p in rx.products:
                        if p.emission_mode in ('prompt', 'total'):
                            neutron = p
                            break
                    else:
                        raise Exception("Couldn't find prompt/total fission neutron")

                # Determine locator for ith energy distribution
                lnw = int(ace.xss[ace.jxs[10] + i_reaction - 1])
                while lnw > 0:
                    # Applicability of this distribution
                    neutron.applicability.append(Tabulated1D.from_ace(
                        ace, ace.jxs[11] + lnw + 2))

                    # Read energy distribution data
                    neutron.distribution.append(AngleEnergy.from_ace(
                        ace, ace.jxs[11], lnw, rx))

                    lnw = int(ace.xss[ace.jxs[11] + lnw - 1])

        else:
            # Elastic scattering
            mt = 2
            rx = cls(mt)

            # Get elastic cross section values
            elastic_xs = ace.xss[ace.jxs[1] + 3*n_grid:ace.jxs[1] + 4*n_grid]

            # Fix negatives -- known issue for Ti46,49,50 in JEFF 3.2
            if np.any(elastic_xs < 0.0):
                warn("Negative elastic scattering cross section found for {}. "
                     "Setting to zero.".format(ace.name))
                elastic_xs[elastic_xs < 0.0] = 0.0

            tabulated_xs = Tabulated1D(grid, elastic_xs)
            tabulated_xs._threshold_idx = 0
            rx.xs[strT] = tabulated_xs

            # No energy distribution for elastic scattering
            neutron = Product('neutron')
            neutron.distribution.append(UncorrelatedAngleEnergy())
            rx.products.append(neutron)

        # ======================================================================
        # ANGLE DISTRIBUTION (FOR UNCORRELATED)

        if i_reaction < ace.nxs[5] + 1:
            # Check if angular distribution data exist
            loc = int(ace.xss[ace.jxs[8] + i_reaction])
            if loc < 0:
                # Angular distribution is given as part of a product
                # angle-energy distribution
                angle_dist = None
            elif loc == 0:
                # Angular distribution is isotropic
                energy = [0.0, grid[-1]]
                mu = Uniform(-1., 1.)
                angle_dist = AngleDistribution(energy, [mu, mu])
            else:
                angle_dist = AngleDistribution.from_ace(ace, ace.jxs[9], loc)

            # Apply angular distribution to each uncorrelated angle-energy
            # distribution
            if angle_dist is not None:
                for d in neutron.distribution:
                    d.angle = angle_dist

        # ======================================================================
        # PHOTON PRODUCTION

        rx.products += _get_photon_products_ace(ace, rx)

        return rx

    @classmethod
    def from_endf(cls, ev, mt):
        """Generate a reaction from an ENDF evaluation

        Parameters
        ----------
        ev : openmc.data.endf.Evaluation or endf.Material
            ENDF evaluation
        mt : int
            The MT value of the reaction to get data for

        Returns
        -------
        rx : openmc.data.Reaction
            Reaction data

        """
        ev = as_evaluation(ev)
        rx = Reaction(mt)

        # Integrated cross section
        if (3, mt) in ev.section:
            file_obj = StringIO(ev.section[3, mt])
            get_head_record(file_obj)
            params, rx.xs['0K'] = get_tab1_record(file_obj)
            rx.q_value = params[1]

        # Get fission product yields (nu) as well as delayed neutron energy
        # distributions
        if mt in FISSION_MTS:
            rx.products, rx.derived_products = _get_fission_products_endf(ev)

        if (6, mt) in ev.section:
            # Product angle-energy distribution
            for product in _get_products(ev, mt):
                if mt in FISSION_MTS and product.particle == 'neutron':
                    rx.products[0].applicability = product.applicability
                    rx.products[0].distribution = product.distribution
                else:
                    rx.products.append(product)

        elif (4, mt) in ev.section or (5, mt) in ev.section:
            # Uncorrelated angle-energy distribution
            neutron = Product('neutron')

            # Note that the energy distribution for MT=455 is read in
            # _get_fission_products_endf rather than here
            if (5, mt) in ev.section:
                file_obj = StringIO(ev.section[5, mt])
                items = get_head_record(file_obj)
                nk = items[4]
                for i in range(nk):
                    params, applicability = get_tab1_record(file_obj)
                    dist = UncorrelatedAngleEnergy()
                    dist.energy = EnergyDistribution.from_endf(file_obj, params)

                    neutron.applicability.append(applicability)
                    neutron.distribution.append(dist)
            elif mt == 2:
                # Elastic scattering -- no energy distribution is given since it
                # can be calulcated analytically
                dist = UncorrelatedAngleEnergy()
                neutron.distribution.append(dist)
            elif mt >= 51 and mt < 91:
                # Level inelastic scattering -- no energy distribution is given
                # since it can be calculated analytically. Here we determine the
                # necessary parameters to create a LevelInelastic object
                dist = UncorrelatedAngleEnergy()

                A = ev.target['mass']
                threshold = (A + 1.)/A*abs(rx.q_value)
                mass_ratio = (A/(A + 1.))**2
                dist.energy = LevelInelastic(threshold, mass_ratio)

                neutron.distribution.append(dist)

            if (4, mt) in ev.section:
                for dist in neutron.distribution:
                    dist.angle = AngleDistribution.from_endf(ev, mt)

            if mt in FISSION_MTS and (5, mt) in ev.section:
                # For fission reactions,
                rx.products[0].applicability = neutron.applicability
                rx.products[0].distribution = neutron.distribution
            else:
                rx.products.append(neutron)

        if (8, mt) in ev.section:
            rx.products += _get_activation_products(ev, rx)

        if (12, mt) in ev.section or (13, mt) in ev.section:
            rx.products += _get_photon_products_endf(ev, rx)

        return rx
