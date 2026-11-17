"""Read-only inspection of native reaction secondary distributions."""
from ctypes import (Structure, POINTER, byref, cast, c_bool, c_char_p, c_double,
                    c_int, c_int64, c_uint64, c_size_t, c_void_p, string_at)
from numbers import Integral

import numpy as np
from numpy.ctypeslib import as_array

from . import _dll
from .error import _error_handler

__all__ = ['get_reaction_products', 'get_emission_distribution',
           'get_angular_distribution', 'sample_angular_distribution',
           'set_angular_distribution', 'sample_emission_energy']


class _Node(Structure):
    _fields_ = [('type', c_int), ('size', c_size_t), ('name', c_char_p),
                ('data', c_void_p)]


_dll.openmc_secondary_snapshot_create.argtypes = [
    c_int, c_int, c_int, c_int, c_bool, POINTER(c_void_p)]
_dll.openmc_secondary_snapshot_create.restype = c_int
_dll.openmc_secondary_snapshot_create.errcheck = _error_handler
_dll.openmc_secondary_snapshot_node.argtypes = [
    c_void_p, c_size_t, POINTER(_Node)]
_dll.openmc_secondary_snapshot_node.restype = c_int
_dll.openmc_secondary_snapshot_node.errcheck = _error_handler
_dll.openmc_secondary_snapshot_free.argtypes = [c_void_p]
_dll.openmc_secondary_snapshot_free.restype = None


def _index(value, name, optional=False):
    if optional and value is None:
        return -1
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f'{name} must be an integer')
    if not 0 <= value <= np.iinfo(np.int32).max:
        raise ValueError(f'{name} must be a nonnegative C integer')
    return int(value)


def _nuclide_index(nuclide):
    if isinstance(nuclide, str):
        from .nuclide import nuclides
        return nuclides[nuclide]._index
    return _index(nuclide, 'nuclide')


def _decode(snapshot, index=0):
    node = _Node()
    _dll.openmc_secondary_snapshot_node(snapshot, index, byref(node))
    tag, count = node.type, node.size
    name = node.name.decode()
    if tag == 0:
        value = None
    elif tag in (1, 2):
        children = cast(node.data, POINTER(c_size_t))
        entries = [_decode(snapshot, children[i]) for i in range(count)]
        value = dict(entries) if tag == 1 else [v for _, v in entries]
    elif tag == 3:
        value = string_at(node.data, count).decode()
    elif tag in (4, 5):
        value = cast(node.data, POINTER(c_int64))[0]
        if tag == 4:
            value = bool(value)
    elif tag == 6:
        value = cast(node.data, POINTER(c_double))[0]
    elif tag in (7, 8):
        dtype = np.int64 if tag == 7 else np.float64
        ctype = c_int64 if tag == 7 else c_double
        value = (as_array(cast(node.data, POINTER(ctype)), (count,)).copy()
                 if count else np.empty(0, dtype=dtype))
    else:
        raise RuntimeError(f'Unknown secondary snapshot value tag: {tag}')
    return name, value


def _snapshot(nuclide, mt, product, distribution, include_data):
    snapshot = c_void_p()
    _dll.openmc_secondary_snapshot_create(
        _nuclide_index(nuclide), _index(mt, 'mt', True),
        _index(product, 'product', True),
        _index(distribution, 'distribution', True), include_data,
        byref(snapshot))
    try:
        return _decode(snapshot)[1]
    finally:
        _dll.openmc_secondary_snapshot_free(snapshot)


def get_reaction_products(nuclide, mt=None, *, include_data=False):
    """Discover secondary distributions on a loaded nuclide's reaction products.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or library index. This does not load nuclear data.
    mt : int, optional
        Reaction MT. Omit to inspect all loaded reactions.
    include_data : bool, optional
        Include full native distribution snapshots as well as metadata.
        Defaults to False, avoiding copies of potentially large joint tables.

    Returns
    -------
    list of dict
        One entry per product, including nuclide, MT, product index, particle,
        emission mode, reference frame, yield, applicability, and distribution
        descriptors. Each descriptor identifies native angle-energy and energy
        types, distribution index, read/mutation capabilities and limitations.
        Empty product distribution lists are retained.

    Notes
    -----
    Requires loaded library data, but not ``simulation_init()``. Library calls
    must be serialized and reads are rejected during a batch. All returned
    arrays are independent copies. No cross sections, yields, applicability,
    frames, or distributions are changed.

    ``mutation_supported`` refers to the entire angle-energy law and is False
    in this read-only interface. ``energy_mutation_supported`` separately
    identifies existing ContinuousTabular editing via ``set_emission_spectrum``;
    that setter still requires simulation initialization and a batch boundary.

    Breakpoints use the one-based ENDF/HDF5 endpoint convention. Stored incident
    interpolation and actual sampling interpolation are reported separately.
    D1S spectra are tagged ``decay-photon`` and remain read-only because their
    storage is borrowed from a shared chain nuclide.
    """
    if not isinstance(include_data, bool):
        raise TypeError('include_data must be a bool')
    return _snapshot(nuclide, mt, None, None, include_data)


def get_emission_distribution(nuclide, mt, product=0, distribution=0):
    """Copy a native secondary angle-energy distribution into a tagged dictionary.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or library index.
    mt : int
        Reaction MT number.
    product : int, optional
        Zero-based reaction product index.
    distribution : int, optional
        Zero-based angle-energy distribution index within the product.

    Returns
    -------
    dict
        Native law identified by ``type``. Uncorrelated laws contain ``angle``
        and ``energy`` (None for an absent energy law). Correlated laws retain
        conditional angular tables inside outgoing-energy tables. Kalbach-Mann
        tables retain PDF/CDF, r, and a. Parametric laws retain native parameters
        and tabulated functions, including incident grids and interpolation.

    Notes
    -----
    Energies are in eV, angular abscissae are scattering cosines, and the frame
    is reported by :func:`get_reaction_products`. PDFs/CDFs are copied exactly
    as stored in C++, without renormalization or reconstruction. Changing this
    dictionary has no effect on transport. See :func:`get_reaction_products`
    for lifecycle requirements and capability discovery.
    """
    _index(mt, 'mt')
    _index(product, 'product')
    _index(distribution, 'distribution')
    products = _snapshot(nuclide, mt, product, distribution, True)
    descriptor = products[0]['distributions'][0]
    if not descriptor['read_supported']:
        raise NotImplementedError(descriptor['read_reason'])
    return descriptor['data']


def get_angular_distribution(nuclide, mt, product=0, distribution=0):
    """Copy an independently defined native angular law.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or library index.
    mt : int
        Reaction MT.
    product, distribution : int, optional
        Zero-based product and distribution indices.

    Returns
    -------
    dict
        ``isotropic`` or ``tabular-angle``. Tabular tables contain ``x`` (mu),
        ``pdf``, ``cdf``, native ``interpolation`` and the loader's ``integral``.
        Incident energies are in eV. Between incident energies the sampler
        selects adjacent tables with linear stochastic weights. No incident
        breakpoints are stored by AngleDistribution. All arrays are copies.

    Raises
    ------
    ValueError
        The law is conditional on outgoing energy. Obtain its native structure
        through :func:`get_emission_distribution` instead.

    Notes
    -----
    The product's reference frame is available through get_reaction_products.
    N-body and decay-photon samplers explicitly emit isotropically in the native
    frame. Library calls must be serialized and made outside transport batches.
    """
    law = get_emission_distribution(nuclide, mt, product, distribution)
    if law['type'] == 'uncorrelated':
        return law['angle']
    if law['type'] in ('nbody', 'decay-photon'):
        return law['angle']
    raise ValueError(
        f"{law['type']} angular law is conditional on outgoing energy; "
        'use get_emission_distribution() for the native joint structure')


_dll.openmc_secondary_sample_angle.argtypes = [
    c_int, c_int, c_int, c_int, c_double, POINTER(c_uint64), c_size_t,
    POINTER(c_double)]
_dll.openmc_secondary_sample_angle.restype = c_int
_dll.openmc_secondary_sample_angle.errcheck = _error_handler


def sample_angular_distribution(nuclide, mt, energy, n_samples, product=0,
                                distribution=0, *, seed=1):
    """Sample an independent angular law using the native transport sampler.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or library index.
    mt : int
        Reaction MT.
    energy : float
        Incident energy in eV (the energy supplied to the native sampler).
    n_samples : int
        Number of scattering cosines to draw.
    product, distribution : int, optional
        Zero-based product and distribution indices.
    seed : int, optional
        Caller-owned unsigned 64-bit RNG seed; transport streams are unaffected.

    Returns
    -------
    numpy.ndarray
        Scattering cosines in the native product frame. These precede transport
        frame conversion, target thermal motion, and azimuthal rotation.

    Notes
    -----
    Requires loaded data and serialized calls outside batches. Correlated and
    Kalbach-Mann laws are rejected because angle depends on outgoing energy.
    This selects a particular law, without sampling product applicability.
    """
    if isinstance(seed, bool) or not isinstance(seed, Integral):
        raise TypeError('seed must be an integer')
    if not 0 <= seed <= np.iinfo(np.uint64).max:
        raise ValueError('seed must fit in an unsigned 64-bit integer')
    mu = np.empty(_index(n_samples, 'n_samples'))
    rng = c_uint64(seed)
    _dll.openmc_secondary_sample_angle(
        _nuclide_index(nuclide), _index(mt, 'mt'), _index(product, 'product'),
        _index(distribution, 'distribution'), energy, byref(rng), mu.size,
        mu.ctypes.data_as(POINTER(c_double)))
    return mu


_dll.openmc_secondary_set_angle.argtypes = [
    c_int, c_int, c_int, c_int, POINTER(c_double), c_int, POINTER(c_int),
    POINTER(c_int), POINTER(c_double), POINTER(c_double), c_int]
_dll.openmc_secondary_set_angle.restype = c_int
_dll.openmc_secondary_set_angle.errcheck = _error_handler


def set_angular_distribution(nuclide, mt, angle, product=0, distribution=0):
    """Replace native uncorrelated angular PDFs on unchanged tabular grids.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide name or index.
    mt : int
        Reaction MT.
    angle : dict
        A tabular-angle snapshot. CDFs are ignored and rebuilt with the native
        histogram or linear-linear rule. All grids/interpolations must match.
    product, distribution : int, optional
        Zero-based indices.

    Notes
    -----
    Requires simulation initialization and serialized calls outside batches.
    Isotropic and conditional laws remain read-only; no implicit conversion.
    All tables are validated before replacing any native table.
    """
    if not isinstance(angle, dict) or angle.get('type') != 'tabular-angle':
        raise ValueError('Only existing tabular-angle laws can be replaced')
    incident = np.ascontiguousarray(angle['incident_energy'], dtype=float)
    tables = angle['tables']
    if incident.ndim != 1 or len(incident) != len(tables) or not len(tables):
        raise ValueError('Invalid angular incident grid/table count')
    offsets, interpolations, coordinates, probabilities = [0], [], [], []
    for table in tables:
        x = np.asarray(table['x'], dtype=float)
        pdf = np.asarray(table['pdf'], dtype=float)
        interp = _index(table['interpolation'], 'interpolation')
        if x.ndim != 1 or pdf.shape != x.shape or len(x) < 2:
            raise ValueError('Invalid angular table shape')
        coordinates.append(x)
        probabilities.append(pdf)
        interpolations.append(interp)
        offsets.append(offsets[-1] + len(x))
    if offsets[-1] > np.iinfo(np.int32).max:
        raise ValueError('Angular tables exceed C integer capacity')
    offsets = np.asarray(offsets, dtype=np.int32)
    interpolations = np.asarray(interpolations, dtype=np.int32)
    x, pdf = np.concatenate(coordinates), np.concatenate(probabilities)
    _dll.openmc_secondary_set_angle(
        _nuclide_index(nuclide), _index(mt, 'mt'), _index(product, 'product'),
        _index(distribution, 'distribution'),
        incident.ctypes.data_as(POINTER(c_double)), len(incident),
        offsets.ctypes.data_as(POINTER(c_int)),
        interpolations.ctypes.data_as(POINTER(c_int)),
        x.ctypes.data_as(POINTER(c_double)), pdf.ctypes.data_as(POINTER(c_double)),
        len(x))


_dll.openmc_secondary_sample_energy.argtypes = [
    c_int, c_int, c_int, c_int, c_double, POINTER(c_uint64), c_size_t,
    POINTER(c_double)]
_dll.openmc_secondary_sample_energy.restype = c_int
_dll.openmc_secondary_sample_energy.errcheck = _error_handler


def sample_emission_energy(nuclide, mt, energy, n_samples, product=0,
                         distribution=0, *, seed=1):
    """Sample an uncorrelated energy law with its native C++ sampler.

    Parameters
    ----------
    nuclide : str or int
        Loaded nuclide.
    mt : int
        Reaction MT.
    energy : float
        Incident energy in eV.
    n_samples : int
        Sample count.
    product, distribution : int, optional
        Zero-based indices.
    seed : int, optional
        Caller-owned unsigned 64-bit seed, independent of transport streams.

    Returns
    -------
    numpy.ndarray
        Outgoing energies in eV in the native product frame.
    """
    if isinstance(seed, bool) or not isinstance(seed, Integral):
        raise TypeError('seed must be an integer')
    if not 0 <= seed <= np.iinfo(np.uint64).max:
        raise ValueError('seed must fit unsigned 64-bit integer')
    result = np.empty(_index(n_samples, 'n_samples'))
    rng = c_uint64(seed)
    _dll.openmc_secondary_sample_energy(
        _nuclide_index(nuclide), _index(mt, 'mt'), _index(product, 'product'),
        _index(distribution, 'distribution'), energy, byref(rng), result.size,
        result.ctypes.data_as(POINTER(c_double)))
    return result
