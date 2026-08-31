from collections.abc import Iterable, Mapping
import copy
from numbers import Integral
import os
import pickle
from warnings import warn

import numpy as np

import openmc
import openmc.mgxs
import openmc.checkvalue as cv
from openmc.checkvalue import PathLike
from ..tallies import ESTIMATOR_TYPES

ROOM_TEMPERATURE_KELVIN = 294.0

class Library:
    """A multi-energy-group and multi-delayed-group cross section library for
    some energy group structure.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for deterministic neutronics calculations.

    This class helps automate the generation of MGXS and MDGXS objects for some
    energy group structure and domain type. The Library serves as a collection
    for MGXS and MDGXS objects with routines to automate the initialization of
    tallies for input files, the loading of tally data from statepoint files,
    data storage, energy group condensation and more.

    Parameters
    ----------
    geometry : openmc.Geometry
        A geometry which has been initialized with a root universe
    by_nuclide : bool
        If true, computes cross sections for each nuclide in each domain
    mgxs_types : Iterable of str
        The types of cross sections in the library (e.g., ['total', 'scatter'])
    name : str, optional
        Name of the multi-group cross section library. Used as a label to
        identify tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    geometry : openmc.Geometry
        An geometry which has been initialized with a root universe
    by_nuclide : bool
        If true, computes cross sections for each nuclide in each domain
    mgxs_types : Iterable of str
        The types of cross sections in the library (e.g., ['total', 'scatter'])
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    domains : Iterable of openmc.Material, openmc.Cell, openmc.Universe or openmc.RegularMesh
        The spatial domain(s) for which MGXS in the Library are computed
    correction : {'P0', None}
        Apply the P0 correction to scattering matrices if set to 'P0'
    transport_correction_ratios : dict or None
        The per-group transport correction ratios
        :math:`r_g = \\sigma_{tr,g} / \\sigma_{t,g}` applied when the library is
        written to an :class:`openmc.MGXSLibrary`. This is a nested dictionary
        keyed first by domain type (e.g., ``'material'``) and then by domain
        ID, whose values are iterables of per-group ratios (one ratio per
        energy group, ordered from fast to thermal to match the
        :class:`openmc.XSdata` group indexing). When a ``'transport'`` or
        ``'nu-transport'`` MGXS type has been tallied,
        :meth:`load_from_statepoint` fills this in automatically from the
        tallied data; entries supplied by the user beforehand are preserved and
        never overwritten. Explicit assignment may also be used to provide or
        edit the ratios directly. For each listed domain the transport-corrected
        total cross section is set to :math:`\\sigma_{tr,g} = r_g \\sigma_{t,g}`
        and the same change is applied to the in-group (diagonal) :math:`P_0`
        element of the scattering matrix, preserving the absorption balance.
        These stored ratios, rather than a correction derived directly from the
        tallies during data extraction, are what determine the transport
        correction, so editing an entry after :meth:`load_from_statepoint`
        changes the correction that is applied. When :attr:`correction` is
        ``None`` the plain :math:`\\sigma_{t,g}` is taken from the ``'total'``
        MGXS already written to the dataset; when :attr:`correction` is
        ``'P0'`` it is recovered from the transport MGXS. A ``'legendre'``
        :attr:`scatter_format` is required, and the ``'P0'`` recovery is only
        performed for isotropic data. Domains without an entry are left
        uncorrected. Defaults to ``None``.

        .. versionadded:: 0.16.1
    scatter_format : {'legendre', 'histogram'}
        Representation of the angular scattering distribution (default is
        'legendre')
    legendre_order : int
        The highest Legendre moment in the scattering matrix; this is used if
        :attr:`ScatterMatrixXS.scatter_format` is 'legendre'. (default is 0)
    histogram_bins : int
        The number of equally-spaced bins for the histogram representation of
        the angular scattering distribution; this is used if
        :attr:`ScatterMatrixXS.scatter_format` is 'histogram'. (default is 16)
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    num_delayed_groups : int
        Number of delayed groups
    num_polar : Integral
        Number of equi-width polar angle bins for angle discretization
    num_azimuthal : Integral
        Number of equi-width azimuthal angle bins for angle discretization
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the domain
        are included.
    estimator : str or None
        The tally estimator used to compute multi-group cross sections.
        If None, the default for each MGXS type is used.
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    all_mgxs : dict
        MGXS objects keyed by domain ID and cross section type
    sp_filename : str
        The filename of the statepoint with tally data used to the
        compute cross sections
    keff : Real or None
        The combined keff from the statepoint file with tally data used to
        compute cross sections (for eigenvalue calculations only)
    name : str, optional
        Name of the multi-group cross section library. Used as a label to
        identify tallies in OpenMC 'tallies.xml' file.
    sparse : bool
        Whether or not the Library's tallies use SciPy's LIL sparse matrix
        format for compressed data storage

    """

    def __init__(self, geometry, by_nuclide=False,
                 mgxs_types=None, name=''):

        self._name = ''
        self._geometry = None
        self._by_nuclide = None
        self._mgxs_types = []
        self._domain_type = None
        self._domains = 'all'
        self._energy_groups = None
        self._num_polar = 1
        self._num_azimuthal = 1
        self._nuclides = None
        self._num_delayed_groups = 0
        self._correction = 'P0'
        self._transport_correction_ratios = None
        self._scatter_format = 'legendre'
        self._legendre_order = 0
        self._histogram_bins = 16
        self._tally_trigger = None
        self._all_mgxs = {}
        self._sp_filename = None
        self._keff = None
        self._sparse = False
        self._estimator = None

        self.name = name
        self.geometry = geometry
        self.by_nuclide = by_nuclide

        if mgxs_types is not None:
            self.mgxs_types = mgxs_types

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self.name
            clone._geometry = self.geometry
            clone._by_nuclide = self.by_nuclide
            clone._mgxs_types = self.mgxs_types
            clone._domain_type = self.domain_type
            clone._domains = copy.deepcopy(self.domains)
            clone._correction = self.correction
            clone._transport_correction_ratios = \
                copy.deepcopy(self._transport_correction_ratios)
            clone._scatter_format = self.scatter_format
            clone._legendre_order = self.legendre_order
            clone._histogram_bins = self.histogram_bins
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._num_polar = self.num_polar
            clone._num_azimuthal = self.num_azimuthal
            clone._nuclides = self._nuclides
            clone._num_delayed_groups = self.num_delayed_groups
            clone._tally_trigger = copy.deepcopy(self.tally_trigger, memo)
            clone._all_mgxs = copy.deepcopy(self.all_mgxs)
            clone._sp_filename = self._sp_filename
            clone._keff = self._keff
            clone._sparse = self.sparse

            clone._all_mgxs = {}
            for domain in self.domains:
                clone.all_mgxs[domain.id] = {}
                for mgxs_type in self.mgxs_types:
                    mgxs = copy.deepcopy(self.all_mgxs[domain.id][mgxs_type])
                    clone.all_mgxs[domain.id][mgxs_type] = mgxs

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def geometry(self):
        return self._geometry

    @geometry.setter
    def geometry(self, geometry):
        cv.check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        cv.check_type('name', name, str)
        self._name = name

    @property
    def mgxs_types(self):
        return self._mgxs_types

    @mgxs_types.setter
    def mgxs_types(self, mgxs_types):
        all_mgxs_types = openmc.mgxs.MGXS_TYPES + openmc.mgxs.MDGXS_TYPES + \
            openmc.mgxs.ARBITRARY_VECTOR_TYPES + \
            openmc.mgxs.ARBITRARY_MATRIX_TYPES
        if mgxs_types == 'all':
            self._mgxs_types = all_mgxs_types
        else:
            cv.check_iterable_type('mgxs_types', mgxs_types, str)
            for mgxs_type in mgxs_types:
                cv.check_value('mgxs_type', mgxs_type, all_mgxs_types)
            self._mgxs_types = mgxs_types

    @property
    def by_nuclide(self):
        return self._by_nuclide

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)

        if by_nuclide and self.domain_type == 'mesh':
            raise ValueError('Unable to create MGXS library by nuclide with '
                             'mesh domain')

        self._by_nuclide = by_nuclide

    @property
    def domain_type(self):
        return self._domain_type

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, openmc.mgxs.DOMAIN_TYPES)

        if self.by_nuclide and domain_type == 'mesh':
            raise ValueError('Unable to create MGXS library by nuclide with '
                             'mesh domain')

        self._domain_type = domain_type

    @property
    def domains(self):
        if self._domains == 'all':
            if self.domain_type == 'material':
                return list(self.geometry.get_all_materials().values())
            elif self.domain_type == 'cell':
                return list(self.geometry.get_all_cells().values())
            elif self.domain_type in 'distribcell':
                return list(self.geometry.get_all_material_cells().values())
            elif self.domain_type == 'universe':
                return list(self.geometry.get_all_universes().values())
            elif self.domain_type == 'mesh':
                raise ValueError('Unable to get domains for Mesh domain type')
            else:
                raise ValueError('Unable to get domains without a domain type')
        else:
            return self._domains

    @domains.setter
    def domains(self, domains):

        # Use all materials, cells or universes in the geometry as domains
        if domains == 'all':
            self._domains = domains

        # User specified a list of material, cell or universe domains
        else:
            if self.domain_type == 'material':
                cv.check_type('domain', domains, Iterable, openmc.Material)
                all_domains = self.geometry.get_all_materials().values()
            elif self.domain_type == 'cell':
                cv.check_type('domain', domains, Iterable, openmc.Cell)
                all_domains = self.geometry.get_all_cells().values()
            elif self.domain_type == 'distribcell':
                cv.check_type('domain', domains, Iterable, openmc.Cell)
                all_domains = self.geometry.get_all_material_cells().values()
            elif self.domain_type == 'universe':
                cv.check_type('domain', domains, Iterable, openmc.Universe)
                all_domains = self.geometry.get_all_universes().values()
            elif self.domain_type == 'mesh':
                cv.check_type('domain', domains, Iterable, openmc.RegularMesh)

                # The mesh and geometry are independent, so set all_domains
                # to the input domains
                all_domains = domains
            else:
                raise ValueError('Unable to set domains with domain '
                                 'type "{}"'.format(self.domain_type))

            # Check that each domain can be found in the geometry
            for domain in domains:
                if domain not in all_domains:
                    raise ValueError('Domain "{}" could not be found in the '
                                     'geometry.'.format(domain))

            self._domains = list(domains)

    @property
    def nuclides(self):
        return self._nuclides

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_iterable_type('nuclides', nuclides, str)
        self._nuclides = nuclides

    @property
    def energy_groups(self):
        return self._energy_groups

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    @property
    def num_delayed_groups(self):
        return self._num_delayed_groups

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):

        cv.check_less_than('num delayed groups', num_delayed_groups,
                           openmc.mgxs.MAX_DELAYED_GROUPS, equality=True)
        cv.check_greater_than('num delayed groups', num_delayed_groups, 0,
                              equality=True)
        self._num_delayed_groups = num_delayed_groups

    @property
    def num_polar(self):
        return self._num_polar

    @num_polar.setter
    def num_polar(self, num_polar):
        cv.check_type('num_polar', num_polar, Integral)
        cv.check_greater_than('num_polar', num_polar, 0)
        self._num_polar = num_polar

    @property
    def num_azimuthal(self):
        return self._num_azimuthal

    @num_azimuthal.setter
    def num_azimuthal(self, num_azimuthal):
        cv.check_type('num_azimuthal', num_azimuthal, Integral)
        cv.check_greater_than('num_azimuthal', num_azimuthal, 0)
        self._num_azimuthal = num_azimuthal

    @property
    def correction(self):
        return self._correction

    @correction.setter
    def correction(self, correction):
        cv.check_value('correction', correction, ('P0', None))

        if self.scatter_format == 'legendre':
            if correction == 'P0' and self.legendre_order > 0:
                msg = 'The P0 correction will be ignored since the ' \
                      'scattering order {} is greater than '\
                      'zero'.format(self.legendre_order)
                warn(msg)
        elif self.scatter_format == 'histogram':
            msg = 'The P0 correction will be ignored since the ' \
                  'scatter format is set to histogram'
            warn(msg)

        self._correction = correction

    @property
    def transport_correction_ratios(self):
        return self._transport_correction_ratios

    @transport_correction_ratios.setter
    def transport_correction_ratios(self, ratios):
        if ratios is None:
            self._transport_correction_ratios = None
            return

        cv.check_type('transport_correction_ratios', ratios, Mapping)

        normalized = {}
        for domain_type, domain_ratios in ratios.items():
            cv.check_value('transport_correction_ratios domain type',
                           domain_type, openmc.mgxs.DOMAIN_TYPES)
            cv.check_type(f'transport_correction_ratios["{domain_type}"]',
                          domain_ratios, Mapping)

            normalized[domain_type] = {}
            for domain_id, group_ratios in domain_ratios.items():
                cv.check_type('transport correction ratio domain ID',
                              domain_id, Integral)

                try:
                    arr = np.asarray(group_ratios, dtype=float)
                except (ValueError, TypeError):
                    raise ValueError(
                        'Transport correction ratios for domain '
                        f'{domain_id} must be a 1-D iterable of real numbers.')

                if arr.ndim != 1:
                    raise ValueError(
                        'Transport correction ratios for domain '
                        f'{domain_id} must be a 1-D iterable (one ratio per '
                        'energy group).')

                if np.any(arr <= 0.0):
                    raise ValueError(
                        'Transport correction ratios for domain '
                        f'{domain_id} must be positive.')

                # Validate the number of groups when it is already known
                if self._energy_groups is not None:
                    cv.check_length(
                        'transport correction ratios for domain '
                        f'{domain_id}', arr, self.num_groups, self.num_groups)

                normalized[domain_type][int(domain_id)] = arr

        self._transport_correction_ratios = normalized

    @property
    def scatter_format(self):
        return self._scatter_format

    @scatter_format.setter
    def scatter_format(self, scatter_format):
        cv.check_value('scatter_format', scatter_format,
                       openmc.mgxs.MU_TREATMENTS)

        if scatter_format == 'histogram' and self.correction == 'P0':
            msg = 'The P0 correction will be ignored since the ' \
                  'scatter format is set to histogram'
            warn(msg)
            self.correction = None

        self._scatter_format = scatter_format

    @property
    def legendre_order(self):
        return self._legendre_order

    @legendre_order.setter
    def legendre_order(self, legendre_order):
        cv.check_type('legendre_order', legendre_order, Integral)
        cv.check_greater_than('legendre_order', legendre_order, 0,
                              equality=True)
        cv.check_less_than('legendre_order', legendre_order, 10, equality=True)

        if self.scatter_format == 'legendre':
            if self.correction == 'P0' and legendre_order > 0:
                msg = 'The P0 correction will be ignored since the ' \
                      'scattering order {} is greater than '\
                      'zero'.format(legendre_order)
                warn(msg, RuntimeWarning)
                self.correction = None
        elif self.scatter_format == 'histogram':
            msg = 'The legendre order will be ignored since the ' \
                  'scatter format is set to histogram'
            warn(msg)

        self._legendre_order = legendre_order

    @property
    def histogram_bins(self):
        return self._histogram_bins

    @histogram_bins.setter
    def histogram_bins(self, histogram_bins):
        cv.check_type('histogram_bins', histogram_bins, Integral)
        cv.check_greater_than('histogram_bins', histogram_bins, 0)

        if self.scatter_format == 'legendre':
            msg = 'The histogram bins will be ignored since the ' \
                  'scatter format is set to legendre'
            warn(msg)
        elif self.scatter_format == 'histogram':
            if self.correction == 'P0':
                msg = 'The P0 correction will be ignored since ' \
                      'a histogram representation of the scattering '\
                      'kernel is requested'
                warn(msg, RuntimeWarning)
                self.correction = None

        self._histogram_bins = histogram_bins

    @property
    def tally_trigger(self):
        return self._tally_trigger

    @tally_trigger.setter
    def tally_trigger(self, tally_trigger):
        cv.check_type('tally trigger', tally_trigger, openmc.Trigger)
        self._tally_trigger = tally_trigger

    @property
    def estimator(self):
        return self._estimator

    @estimator.setter
    def estimator(self, estimator):
        cv.check_value('estimator', estimator, ESTIMATOR_TYPES)
        self._estimator = estimator

    @property
    def num_groups(self):
        return self.energy_groups.num_groups

    @property
    def all_mgxs(self):
        return self._all_mgxs

    @property
    def sp_filename(self):
        return self._sp_filename

    @property
    def keff(self):
        return self._keff

    @property
    def sparse(self):
        return self._sparse

    @sparse.setter
    def sparse(self, sparse):
        """Convert tally data from NumPy arrays to SciPy list of lists (LIL)
        sparse matrices, and vice versa.

        This property may be used to reduce the amount of data in memory during
        tally data processing. The tally data will be stored as SciPy LIL
        matrices internally within the Tally object. All tally data access
        properties and methods will return data as a dense NumPy array.

        """

        cv.check_type('sparse', sparse, bool)

        # Sparsify or densify each MGXS in the Library
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                mgxs.sparse = self.sparse

        self._sparse = sparse

    def build_library(self):
        """Initialize MGXS objects in each domain and for each reaction type
        in the library.

        This routine will populate the all_mgxs instance attribute dictionary
        with MGXS subclass objects keyed by each domain ID (e.g., Material IDs)
        and cross section type (e.g., 'nu-fission', 'total', etc.).

        """

        # Initialize MGXS for each domain and mgxs type and store in dictionary
        for domain in self.domains:
            self.all_mgxs[domain.id] = {}
            for mgxs_type in self.mgxs_types:
                if mgxs_type in openmc.mgxs.MDGXS_TYPES:
                    mgxs = openmc.mgxs.MDGXS.get_mgxs(
                        mgxs_type, name=self.name, num_polar=self.num_polar,
                        num_azimuthal=self.num_azimuthal)
                else:
                    mgxs = openmc.mgxs.MGXS.get_mgxs(
                        mgxs_type, name=self.name, num_polar=self.num_polar,
                        num_azimuthal=self.num_azimuthal)

                mgxs.domain = domain
                mgxs.domain_type = self.domain_type
                mgxs.energy_groups = self.energy_groups
                mgxs.by_nuclide = self.by_nuclide
                if self.estimator is not None:
                    mgxs.estimator = self.estimator

                if mgxs_type in openmc.mgxs.MDGXS_TYPES:
                    if self.num_delayed_groups == 0:
                        mgxs.delayed_groups = None
                    else:
                        delayed_groups \
                            = list(range(1, self.num_delayed_groups + 1))
                        mgxs.delayed_groups = delayed_groups

                # If a tally trigger was specified, add it to the MGXS
                if self.tally_trigger is not None:
                    mgxs.tally_trigger = self.tally_trigger

                # Specify whether to use a transport ('P0') correction
                if isinstance(mgxs, openmc.mgxs.ScatterMatrixXS):
                    mgxs.correction = self.correction
                    mgxs.scatter_format = self.scatter_format
                    mgxs.legendre_order = self.legendre_order
                    mgxs.histogram_bins = self.histogram_bins

                if self.by_nuclide:
                    try:
                        domain_nuclides = domain.get_nuclides()
                    except AttributeError:
                        domain_nuclides = None
                    if self.nuclides:
                        if domain_nuclides:
                            mgxs.nuclides = [
                                nuclide for nuclide in self.nuclides
                                if nuclide in domain_nuclides
                            ] + ["total"]
                        else:
                            mgxs.nuclides = self.nuclides

                self.all_mgxs[domain.id][mgxs_type] = mgxs

    def add_to_tallies(self, tallies, merge=True):
        """Add tallies from all MGXS objects to a tallies object.

        NOTE: This assumes that :meth:`Library.build_library` has been called

        Parameters
        ----------
        tallies : openmc.Tallies
            A Tallies collection to add each MGXS' tallies to generate a
            'tallies.xml' input file for OpenMC
        merge : bool
            Indicate whether tallies should be merged when possible. Defaults
            to True.

        """

        cv.check_type('tallies', tallies, openmc.Tallies)

        # Add tallies from each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)

                if mgxs_type in openmc.mgxs.MDGXS_TYPES:
                    if self.num_delayed_groups == 0:
                        mgxs.delayed_groups = None
                    else:
                        mgxs.delayed_groups \
                            = list(range(1, self.num_delayed_groups + 1))

                for tally in mgxs.tallies.values():
                    tallies.append(tally, merge=merge)

    def add_to_tallies_file(self, tallies_file, merge=True):
        warn(
            "The Library.add_to_tallies_file(...) method has been renamed to"
            "add_to_tallies(...) and will be removed in a future version of "
            "OpenMC.", FutureWarning
        )
        self.add_to_tallies(tallies_file, merge=merge)

    def load_from_statepoint(self, statepoint):
        """Extracts tallies in an OpenMC StatePoint with the data needed to
        compute multi-group cross sections.

        This method is needed to compute cross section data from tallies
        in an OpenMC StatePoint object.

        NOTE: The statepoint must first be linked with an OpenMC Summary object.

        Parameters
        ----------
        statepoint : openmc.StatePoint
            An OpenMC StatePoint object with tally data

        Raises
        ------
        ValueError
            When this method is called with a statepoint that has not been
            linked with a summary object.

        """

        cv.check_type('statepoint', statepoint, openmc.StatePoint)

        if statepoint.summary is None:
            msg = 'Unable to load data from a statepoint which has not been ' \
                  'linked with a summary file'
            raise ValueError(msg)

        self._sp_filename = statepoint._f.filename
        self._geometry = statepoint.summary.geometry
        self._atomic_weight_ratios = statepoint.summary.nuclides

        if statepoint.run_mode == 'eigenvalue':
            self._keff = statepoint.keff.n

        # Load tallies for each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                mgxs.load_from_statepoint(statepoint)
                mgxs.sparse = self.sparse

        # Record the transport correction ratios implied by any tallied
        # transport-corrected total cross section (see the
        # transport_correction_ratios attribute)
        self._store_computed_transport_correction_ratios()

    def get_mgxs(self, domain, mgxs_type):
        """Return the MGXS object for some domain and reaction rate type.

        This routine searches the library for an MGXS object for the spatial
        domain and reaction rate type requested by the user.

        NOTE: This routine must be called after the build_library() routine.

        Parameters
        ----------
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh or Integral
            The material, cell, or universe object of interest (or its ID)
        mgxs_type : str
            The type of multi-group cross section object to return; allowable
            values are those MGXS to the Library and present in the
            mgxs_types attribute.

        Returns
        -------
        openmc.mgxs.MGXS
            The MGXS object for the requested domain and reaction rate type

        Raises
        ------
        ValueError
            If no MGXS object can be found for the requested domain or
            multi-group cross section type

        """

        if self.domain_type == 'material':
            cv.check_type('domain', domain, (openmc.Material, Integral))
        elif self.domain_type == 'cell' or self.domain_type == 'distribcell':
            cv.check_type('domain', domain, (openmc.Cell, Integral))
        elif self.domain_type == 'universe':
            cv.check_type('domain', domain, (openmc.Universe, Integral))
        elif self.domain_type == 'mesh':
            cv.check_type('domain', domain, (openmc.RegularMesh, Integral))

        # Check that requested domain is included in library
        if isinstance(domain, Integral):
            domain_id = domain
            for domain in self.domains:
                if domain_id == domain.id:
                    break
            else:
                msg = 'Unable to find MGXS for "{0}" "{1}" in ' \
                      'library'.format(self.domain_type, domain_id)
                raise ValueError(msg)
        else:
            domain_id = domain.id

        # Check that requested domain is included in library
        if mgxs_type not in self.mgxs_types:
            msg = f'Unable to find MGXS type "{mgxs_type}"'
            raise ValueError(msg)

        return self.all_mgxs[domain_id][mgxs_type]

    def get_condensed_library(self, coarse_groups):
        """Construct an energy-condensed version of this library.

        This routine condenses each of the multi-group cross sections in the
        library to a coarse energy group structure. NOTE: This routine must
        be called after the load_from_statepoint(...) routine loads the tallies
        from the statepoint into each of the cross sections.

        Parameters
        ----------
        coarse_groups : openmc.mgxs.EnergyGroups
            The coarse energy group structure of interest

        Returns
        -------
        openmc.mgxs.Library
            A new multi-group cross section library condensed to the group
            structure of interest

        Raises
        ------
        ValueError
            When this method is called before a statepoint has been loaded

        See also
        --------
        MGXS.get_condensed_xs(coarse_groups)

        """

        if self.sp_filename is None:
            msg = 'Unable to get a condensed coarse group cross section ' \
                  'library since the statepoint has not yet been loaded'
            raise ValueError(msg)

        cv.check_type('coarse_groups', coarse_groups, openmc.mgxs.EnergyGroups)
        cv.check_less_than('coarse groups', coarse_groups.num_groups,
                           self.num_groups, equality=True)
        cv.check_value('upper coarse energy', coarse_groups.group_edges[-1],
                       [self.energy_groups.group_edges[-1]])
        cv.check_value('lower coarse energy', coarse_groups.group_edges[0],
                       [self.energy_groups.group_edges[0]])

        # Clone this Library to initialize the condensed version
        condensed_library = copy.deepcopy(self)
        condensed_library.energy_groups = coarse_groups

        # Condense the MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = condensed_library.get_mgxs(domain, mgxs_type)
                condensed_mgxs = mgxs.get_condensed_xs(coarse_groups)
                condensed_library.all_mgxs[domain.id][mgxs_type] = condensed_mgxs

        return condensed_library

    def get_subdomain_avg_library(self):
        """Construct a subdomain-averaged version of this library.

        This routine averages each multi-group cross section across distribcell
        instances. The method performs spatial homogenization to compute the
        scalar flux-weighted average cross section across the subdomains.

        NOTE: This method is only relevant for distribcell domain types and
        simplys returns a deep copy of the library for all other domains types.

        Returns
        -------
        openmc.mgxs.Library
            A new multi-group cross section library averaged across subdomains

        Raises
        ------
        ValueError
            When this method is called before a statepoint has been loaded

        See also
        --------
        MGXS.get_subdomain_avg_xs(subdomains)

        """

        if self.sp_filename is None:
            msg = 'Unable to get a subdomain-averaged cross section ' \
                  'library since the statepoint has not yet been loaded'
            raise ValueError(msg)

        # Clone this Library to initialize the subdomain-averaged version
        subdomain_avg_library = copy.deepcopy(self)

        if subdomain_avg_library.domain_type == 'distribcell':
            subdomain_avg_library.domain_type = 'cell'
        else:
            return subdomain_avg_library

        # Subdomain average the MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = subdomain_avg_library.get_mgxs(domain, mgxs_type)
                if mgxs.domain_type == 'distribcell':
                    avg_mgxs = mgxs.get_subdomain_avg_xs()
                    subdomain_avg_library.all_mgxs[domain.id][mgxs_type] = avg_mgxs

        return subdomain_avg_library

    def build_hdf5_store(self, filename='mgxs.h5', directory='mgxs',
                         subdomains='all', nuclides='all', xs_type='macro',
                         row_column='inout', libver='earliest'):
        """Export the multi-group cross section library to an HDF5 binary file.

        This method constructs an HDF5 file which stores the library's
        multi-group cross section data. The data is stored in a hierarchy of
        HDF5 groups from the domain type, domain id, subdomain id (for
        distribcell domains), nuclides and cross section types. Two datasets for
        the mean and standard deviation are stored for each subdomain entry in
        the HDF5 file. The number of groups is stored as a file attribute.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file. Defaults to 'mgxs.h5'.
        directory : str
            Directory for the HDF5 file. Defaults to 'mgxs'.
        subdomains : {'all', 'avg'}
            Report all subdomains or the average of all subdomain cross sections
            in the report. Defaults to 'all'.
        nuclides : {'all', 'sum'}
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        row_column: {'inout', 'outin'}
            Store scattering matrices indexed first by incoming group and
            second by outgoing group ('inout'), or vice versa ('outin').
            Defaults to 'inout'.
        libver : {'earliest', 'latest'}
            Compatibility mode for the HDF5 file. 'latest' will produce files
            that are less backwards compatible but have performance benefits.

        Raises
        ------
        ValueError
            When this method is called before a statepoint has been loaded

        See also
        --------
        MGXS.build_hdf5_store(filename, directory, xs_type)

        """

        if self.sp_filename is None:
            msg = 'Unable to export multi-group cross section library ' \
                  'since a statepoint has not yet been loaded'
            raise ValueError(msg)

        cv.check_type('filename', filename, (str, PathLike))
        cv.check_type('directory', directory, str)

        import h5py

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Add an attribute for the number of energy groups to the HDF5 file
        full_filename = os.path.join(directory, filename)
        full_filename = full_filename.replace(' ', '-')
        f = h5py.File(full_filename, 'w', libver=libver)
        f.attrs['# groups'] = self.num_groups
        f.close()

        # Export MGXS for each domain and mgxs type to an HDF5 file
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.all_mgxs[domain.id][mgxs_type]

                if subdomains == 'avg':
                    mgxs = mgxs.get_subdomain_avg_xs()

                mgxs.build_hdf5_store(filename, directory, xs_type=xs_type,
                                      nuclides=nuclides, row_column=row_column)

    def dump_to_file(self, filename='mgxs', directory='mgxs'):
        """Store this Library object in a pickle binary file.

        Parameters
        ----------
        filename : str
            Filename for the pickle file. Defaults to 'mgxs'.
        directory : str
            Directory for the pickle file. Defaults to 'mgxs'.

        See also
        --------
        Library.load_from_file(filename, directory)

        """

        cv.check_type('filename', filename, (str, PathLike))
        cv.check_type('directory', directory, str)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_filename = os.path.join(directory, f'{filename}.pkl')
        full_filename = full_filename.replace(' ', '-')

        # Load and return pickled Library object
        with open(full_filename, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load_from_file(filename='mgxs', directory='mgxs'):
        """Load a Library object from a pickle binary file.

        Parameters
        ----------
        filename : str
            Filename for the pickle file. Defaults to 'mgxs'.
        directory : str
            Directory for the pickle file. Defaults to 'mgxs'.

        Returns
        -------
        openmc.mgxs.Library
            A Library object loaded from the pickle binary file

        See also
        --------
        Library.dump_to_file(mgxs_lib, filename, directory)

        """

        cv.check_type('filename', filename, (str, PathLike))
        cv.check_type('directory', directory, str)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_filename = os.path.join(directory, filename + '.pkl')
        full_filename = full_filename.replace(' ', '-')

        # Load and return pickled Library object
        with open(full_filename, 'rb') as f:
            return pickle.load(f)

    def _uncorrected_total_xs(self, transport_mgxs, nuclides, xs_type='macro',
                              subdomains='all'):
        """Recover the plain total cross section from a transport MGXS.

        A :class:`openmc.mgxs.TransportXS` tallies both the total reaction rate
        and the flux, so the un-corrected total cross section can be recovered
        as ``tallies['total'] / tallies['flux (tracklength)']``. The transport
        MGXS's own :meth:`~openmc.mgxs.MGXS.get_xs` machinery is reused (by
        temporarily overriding its cached cross-section tally) so that the
        energy-group ordering and nuclide handling match the transport-
        corrected cross section exactly.

        Parameters
        ----------
        transport_mgxs : openmc.mgxs.TransportXS
            The transport MGXS whose plain total cross section is wanted
        nuclides : str or Iterable of str
            The ``nuclides`` argument to forward to
            :meth:`~openmc.mgxs.MGXS.get_xs` (e.g., ``'total'``, ``'sum'`` or a
            list with a single nuclide name)
        xs_type : {'macro', 'micro'}
            The ``xs_type`` argument to forward to
            :meth:`~openmc.mgxs.MGXS.get_xs`
        subdomains : Iterable of int or 'all'
            The ``subdomains`` argument to forward to
            :meth:`~openmc.mgxs.MGXS.get_xs`

        Returns
        -------
        numpy.ndarray
            The plain total cross section (one value per energy group)

        """

        saved_xs_tally = transport_mgxs._xs_tally
        saved_rxn_rate_tally = transport_mgxs._rxn_rate_tally
        try:
            transport_mgxs._xs_tally = (
                transport_mgxs.tallies['total'] /
                transport_mgxs.tallies['flux (tracklength)'])
            transport_mgxs._compute_xs()
            return transport_mgxs.get_xs(nuclides=nuclides, xs_type=xs_type,
                                         subdomains=subdomains)
        finally:
            transport_mgxs._xs_tally = saved_xs_tally
            transport_mgxs._rxn_rate_tally = saved_rxn_rate_tally

    def _store_computed_transport_correction_ratios(self):
        """Populate transport_correction_ratios from tallied transport data.

        When a transport-type MGXS (``'transport'`` or ``'nu-transport'``) has
        been tallied, the per-group transport correction ratio
        :math:`r_g = \\sigma_{tr,g} / \\sigma_{t,g}` is computed for each domain
        and stored in :attr:`transport_correction_ratios`. The transport-
        corrected total :math:`\\sigma_{tr}` is the transport MGXS itself, while
        the plain total :math:`\\sigma_t` is recovered from the same MGXS's flux
        and total tallies (see :meth:`_uncorrected_total_xs`).

        Ratios that the user already supplied for a domain are never
        overwritten. Only isotropic, single-subdomain data (one ratio per
        energy group) is handled; angle-dependent or multi-subdomain domains
        are skipped.

        """

        # One ratio per group is only well defined for isotropic data
        if self.num_polar > 1 or self.num_azimuthal > 1:
            return

        # Use the transport-corrected total that matches the scattering
        # multiplicity treatment when both flavors are available
        if 'nu-transport' in self.mgxs_types:
            transport_type = 'nu-transport'
        elif 'transport' in self.mgxs_types:
            transport_type = 'transport'
        else:
            return

        nuclides = 'sum' if self.by_nuclide else 'total'

        # Start from any existing (e.g., user-provided) ratios so they are
        # preserved, and never overwrite an entry the user already set
        ratios = copy.deepcopy(self._transport_correction_ratios) or {}
        domain_ratios = ratios.get(self.domain_type, {})

        for domain in self.domains:
            if domain.id in domain_ratios:
                continue

            transport_mgxs = self.get_mgxs(domain, transport_type)
            sigma_tr = np.asarray(
                transport_mgxs.get_xs(nuclides=nuclides, xs_type='macro'),
                dtype=float)
            sigma_t = np.asarray(
                self._uncorrected_total_xs(transport_mgxs, nuclides),
                dtype=float)

            # Only isotropic, single-subdomain data yields one ratio per group
            if sigma_tr.shape != (self.num_groups,) or \
                    sigma_t.shape != (self.num_groups,):
                continue

            # Default to unity (no correction) where the total vanishes
            ratio = np.ones(self.num_groups)
            nonzero = sigma_t > 0.0
            ratio[nonzero] = sigma_tr[nonzero] / sigma_t[nonzero]
            domain_ratios[domain.id] = ratio

        if domain_ratios:
            ratios[self.domain_type] = domain_ratios
            self._transport_correction_ratios = ratios

    def _get_transport_correction_ratios(self, domain):
        """Return the per-group transport correction ratios for a domain.

        Parameters
        ----------
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            The domain of interest

        Returns
        -------
        numpy.ndarray or None
            The transport correction ratios (one per energy group) for the
            domain, or None if none were provided.

        """

        if not self._transport_correction_ratios:
            return None

        domain_ratios = self._transport_correction_ratios.get(self.domain_type)
        if not domain_ratios:
            return None

        return domain_ratios.get(domain.id)

    def _apply_transport_correction_ratios(self, xsdata, domain, temperature,
                                            nuclide='total', xs_type='macro',
                                            subdomain='all'):
        """Apply the stored transport correction ratios to an XSdata object.

        For each energy group ``g`` the transport-corrected total cross section
        is set to :math:`r_g \\sigma_{t,g}`, where :math:`\\sigma_{t,g}` is the
        plain total cross section and :math:`r_g` is the stored ratio. The same
        change is applied to the in-group (diagonal) :math:`P_0` element of the
        scattering matrix, so the absorption balance (total minus out-scatter)
        is preserved.

        Because the total and the scattering diagonal are shifted by the same
        amount, the stored ratios rather than a correction derived directly
        from the tallies during data extraction are what determine the
        transport correction. Editing an entry in
        :attr:`transport_correction_ratios` after
        :meth:`load_from_statepoint` therefore changes the correction applied
        here.

        When :attr:`correction` is ``None`` the total already in ``xsdata`` is
        the plain :math:`\\sigma_t`. When :attr:`correction` is ``'P0'`` the
        total is the tally-based transport-corrected :math:`\\sigma_{tr}`, so
        the plain :math:`\\sigma_t` is recovered from the transport MGXS (see
        :meth:`_uncorrected_total_xs`); this recovery is only performed for
        isotropic data, matching the ratios computed by
        :meth:`_store_computed_transport_correction_ratios`.

        Parameters
        ----------
        xsdata : openmc.XSdata
            The dataset to correct in place
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            The domain the dataset describes
        temperature : float
            Temperature (in Kelvin) of the data to correct
        nuclide : str
            The nuclide the dataset describes (or ``'total'`` for material-wise
            data); used to recover the plain total cross section
        xs_type : {'macro', 'micro'}
            Whether the dataset holds macroscopic or microscopic cross sections
        subdomain : Iterable of int or 'all'
            The subdomain the dataset describes (for mesh domains)

        """

        ratios = self._get_transport_correction_ratios(domain)
        if ratios is None:
            return

        # A legendre scattering matrix is needed to correct the P0 diagonal,
        # and exactly one ratio per energy group is required.
        if self.scatter_format != 'legendre':
            return
        if len(ratios) != self.num_groups:
            return

        i = xsdata._temperature_index(temperature)
        current_total = xsdata._total[i]
        if current_total is None:
            return

        # Determine the plain total cross section sigma_t.
        if self.correction is None:
            # The total already in xsdata is the plain total.
            sigma_t = current_total
        else:
            # The total in xsdata is the tally-corrected transport total, so
            # recover the plain total from the transport MGXS. The stored
            # ratios are only computed for isotropic data.
            if xsdata.representation == 'angle':
                return
            if 'nu-transport' in self.mgxs_types:
                transport_type = 'nu-transport'
            elif 'transport' in self.mgxs_types:
                transport_type = 'transport'
            else:
                return
            transport_mgxs = self.get_mgxs(domain, transport_type)
            sigma_t = np.asarray(
                self._uncorrected_total_xs(transport_mgxs, [nuclide], xs_type,
                                           subdomain), dtype=float)
            if sigma_t.size != np.size(current_total):
                return
            sigma_t = sigma_t.reshape(np.shape(current_total))

        # Set the transport-corrected total to r * sigma_t and shift the P0
        # scattering diagonal by the same amount. ratios has shape (G,) and
        # broadcasts against the trailing group axis for both the isotropic and
        # angle representations.
        target_total = ratios * sigma_t
        delta = target_total - current_total
        xsdata._total[i] = target_total

        scatter_matrix = xsdata._scatter_matrix[i]
        if scatter_matrix is not None:
            groups = np.arange(self.num_groups)
            if xsdata.representation == 'angle':
                scatter_matrix[:, :, groups, groups, 0] += delta
            else:
                scatter_matrix[groups, groups, 0] += delta

    def get_xsdata(self, domain, xsdata_name, nuclide='total', xs_type='macro',
                   subdomain=None, apply_domain_chi=False, temperature=ROOM_TEMPERATURE_KELVIN):
        """Generates an openmc.XSdata object describing a multi-group cross section
        dataset for writing to an openmc.MGXSLibrary object.

        Note that this method does not build an XSdata
        object with nested temperature tables.  The temperature of each
        XSdata object will be left at the default value of 294K.

        Parameters
        ----------
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.RegularMesh
            The domain for spatial homogenization
        xsdata_name : str
            Name to apply to the "xsdata" entry produced by this method
        nuclide : str
            A nuclide name string (e.g., 'U235').  Defaults to 'total' to
            obtain a material-wise macroscopic cross section.
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'. If the Library object is not tallied by
            nuclide this will be set to 'macro' regardless.
        subdomain : iterable of int
            This parameter is not used unless using a mesh domain. In that
            case, the subdomain is an [i,j,k] index (1-based indexing) of the
            mesh cell of interest in the openmc.RegularMesh object.  Note:
            this parameter currently only supports subdomains within a mesh,
            and not the subdomains of a distribcell.
        apply_domain_chi : bool
            This parameter sets whether (True) or not (False) the
            domain-averaged values of chi, chi-prompt, and chi-delayed are to
            be applied to each of the nuclide-dependent fission energy spectra
            of a domain. In effect, if this is True, then every nuclide in the
            domain receives the same flux-weighted Chi. This is useful for
            downstream multigroup solvers that precompute a material-specific
            chi before the transport solve provides group-wise fluxes. Defaults
            to False.
        temperature : float, optional
            The temperature to set in the XSdata object. Defaults to 294 K
            (room temperature).

        Returns
        -------
        xsdata : openmc.XSdata
            Multi-Group Cross Section dataset object.

        Raises
        ------
        ValueError
            When the Library object is initialized with insufficient types of
            cross sections for the Library.

        See also
        --------
        Library.create_mg_library()

        """

        cv.check_type('domain', domain, (openmc.Material, openmc.Cell,
                                         openmc.Universe, openmc.RegularMesh))
        cv.check_type('xsdata_name', xsdata_name, str)
        cv.check_type('nuclide', nuclide, str)
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])
        if subdomain is not None:
            cv.check_iterable_type('subdomain', subdomain, Integral,
                                   max_depth=3)

        # Make sure statepoint has been loaded
        if self._sp_filename is None:
            msg = 'A StatePoint must be loaded before calling ' \
                  'the create_mg_library() function'
            raise ValueError(msg)

        # If gathering material-specific data, set the xs_type to macro
        if not self.by_nuclide:
            xs_type = 'macro'

        # Build & add metadata to XSdata object
        name = xsdata_name
        if nuclide != 'total':
            name += '_' + nuclide
        if self.num_polar > 1 or self.num_azimuthal > 1:
            representation = 'angle'
        else:
            representation = 'isotropic'
        xsdata = openmc.XSdata(name, self.energy_groups,
                               temperatures=[temperature],
                               representation=representation)
        xsdata.num_delayed_groups = self.num_delayed_groups
        if self.num_polar > 1 or self.num_azimuthal > 1:
            xsdata.num_polar = self.num_polar
            xsdata.num_azimuthal = self.num_azimuthal

        if nuclide != 'total':
            xsdata.atomic_weight_ratio = self._atomic_weight_ratios[nuclide]

        if subdomain is None:
            subdomain = 'all'
        else:
            subdomain = [subdomain]

        # Now get xs data itself
        if 'nu-transport' in self.mgxs_types and self.correction == 'P0':
            mymgxs = self.get_mgxs(domain, 'nu-transport')
            xsdata.set_total_mgxs(mymgxs, temperature=temperature,
                                  xs_type=xs_type,
                                  nuclide=[nuclide],
                                  subdomain=subdomain)

        elif 'transport' in self.mgxs_types and self.correction == 'P0':
            mymgxs = self.get_mgxs(domain, 'transport')
            xsdata.set_total_mgxs(mymgxs, temperature=temperature,
                                  xs_type=xs_type,
                                  nuclide=[nuclide],
                                  subdomain=subdomain)

        elif 'total' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'total')
            xsdata.set_total_mgxs(mymgxs, temperature=temperature,
                                  xs_type=xs_type,
                                  nuclide=[nuclide],
                                  subdomain=subdomain)

        if 'absorption' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'absorption')
            xsdata.set_absorption_mgxs(mymgxs,
                                       temperature=temperature,
                                       xs_type=xs_type,
                                       nuclide=[nuclide],
                                       subdomain=subdomain)

        if 'fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'fission')
            xsdata.set_fission_mgxs(mymgxs, temperature=temperature,
                                    xs_type=xs_type,
                                    nuclide=[nuclide],
                                    subdomain=subdomain)

        if 'kappa-fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'kappa-fission')
            xsdata.set_kappa_fission_mgxs(mymgxs,
                                          temperature=temperature,
                                          xs_type=xs_type,
                                          nuclide=[nuclide],
                                          subdomain=subdomain)

        if 'inverse-velocity' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'inverse-velocity')
            xsdata.set_inverse_velocity_mgxs(mymgxs,
                                             temperature=temperature,
                                             xs_type=xs_type,
                                             nuclide=[nuclide],
                                             subdomain=subdomain)

        if 'nu-fission matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'nu-fission matrix')
            xsdata.set_nu_fission_mgxs(mymgxs,
                                       temperature=temperature,
                                       xs_type=xs_type,
                                       nuclide=[nuclide],
                                       subdomain=subdomain)

        if 'chi' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'chi')
            if apply_domain_chi and nuclide != "total":
                nuc = "sum"
            else:
                nuc = nuclide
            xsdata.set_chi_mgxs(mymgxs, temperature=temperature,
                                xs_type=xs_type,
                                nuclide=[nuc],
                                subdomain=subdomain)

        if 'chi-prompt' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'chi-prompt')
            if apply_domain_chi and nuclide != "total":
                nuc = "sum"
            else:
                nuc = nuclide
            xsdata.set_chi_prompt_mgxs(mymgxs, temperature=temperature,
                                       xs_type=xs_type,
                                       nuclide=[nuc],
                                       subdomain=subdomain)

        if 'chi-delayed' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'chi-delayed')
            if apply_domain_chi and nuclide != "total":
                nuc = "sum"
            else:
                nuc = nuclide
            xsdata.set_chi_delayed_mgxs(mymgxs, temperature=temperature,
                                        xs_type=xs_type,
                                        nuclide=[nuc],
                                        subdomain=subdomain)

        if 'nu-fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'nu-fission')
            xsdata.set_nu_fission_mgxs(mymgxs, temperature=temperature,
                                       xs_type=xs_type,
                                       nuclide=[nuclide],
                                       subdomain=subdomain)

        if 'prompt-nu-fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'prompt-nu-fission')
            xsdata.set_prompt_nu_fission_mgxs(mymgxs, temperature=temperature,
                                              xs_type=xs_type,
                                              nuclide=[nuclide],
                                              subdomain=subdomain)

        if 'prompt-nu-fission matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'prompt-nu-fission matrix')
            xsdata.set_prompt_nu_fission_mgxs(mymgxs, temperature=temperature,
                                              xs_type=xs_type,
                                              nuclide=[nuclide],
                                              subdomain=subdomain)

        if 'delayed-nu-fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'delayed-nu-fission')
            xsdata.set_delayed_nu_fission_mgxs(mymgxs, temperature=temperature,
                                               xs_type=xs_type,
                                               nuclide=[nuclide],
                                               subdomain=subdomain)

        if 'delayed-nu-fission matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'delayed-nu-fission matrix')
            xsdata.set_delayed_nu_fission_mgxs(mymgxs, temperature=temperature,
                                               xs_type=xs_type,
                                               nuclide=[nuclide],
                                               subdomain=subdomain)

        if 'beta' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'beta')
            xsdata.set_beta_mgxs(mymgxs, temperature=temperature, xs_type=xs_type,
                                 nuclide=[nuclide], subdomain=subdomain)

        if 'decay-rate' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'decay-rate')
            xsdata.set_decay_rate_mgxs(mymgxs, temperature=temperature, xs_type=xs_type,
                                       nuclide=[nuclide], subdomain=subdomain)

        # If multiplicity matrix is available, prefer that
        if 'multiplicity matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'multiplicity matrix')
            xsdata.set_multiplicity_matrix_mgxs(mymgxs, temperature=temperature,
                                                xs_type=xs_type,
                                                nuclide=[nuclide],
                                                subdomain=subdomain)
            using_multiplicity = True

        # multiplicity will fall back to using scatter and nu-scatter
        elif 'scatter matrix' in self.mgxs_types and \
             'nu-scatter matrix' in self.mgxs_types:
            scatt_mgxs = self.get_mgxs(domain, 'scatter matrix')
            nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
            xsdata.set_multiplicity_matrix_mgxs(nuscatt_mgxs, scatt_mgxs,
                                                temperature=temperature,
                                                xs_type=xs_type,
                                                nuclide=[nuclide],
                                                subdomain=subdomain)
            using_multiplicity = True

        # multiplicity will fall back to using scatter and nu-scatter
        elif 'consistent scatter matrix' in self.mgxs_types and \
             'consistent nu-scatter matrix' in self.mgxs_types:
            scatt_mgxs = self.get_mgxs(domain, 'consistent scatter matrix')
            nuscatt_mgxs = \
                self.get_mgxs(domain, 'consistent nu-scatter matrix')
            xsdata.set_multiplicity_matrix_mgxs(nuscatt_mgxs, scatt_mgxs,
                                                temperature=temperature,
                                                xs_type=xs_type,
                                                nuclide=[nuclide],
                                                subdomain=subdomain)
            using_multiplicity = True

        else:
            using_multiplicity = False

        if using_multiplicity:
            if 'nu-scatter matrix' in self.mgxs_types:
                nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
            else:
                nuscatt_mgxs = \
                    self.get_mgxs(domain, 'consistent nu-scatter matrix')
            xsdata.set_scatter_matrix_mgxs(nuscatt_mgxs, temperature=temperature,
                                           xs_type=xs_type,
                                           nuclide=[nuclide],
                                           subdomain=subdomain)
        else:
            if 'nu-scatter matrix' in self.mgxs_types or \
                    'consistent nu-scatter matrix' in self.mgxs_types:
                if 'nu-scatter matrix' in self.mgxs_types:
                    nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
                else:
                    nuscatt_mgxs = \
                        self.get_mgxs(domain, 'consistent nu-scatter matrix')
                xsdata.set_scatter_matrix_mgxs(nuscatt_mgxs, temperature=temperature,
                                               xs_type=xs_type,
                                               nuclide=[nuclide],
                                               subdomain=subdomain)

                # Since we are not using multiplicity, then
                # scattering multiplication (nu-scatter) must be
                # accounted for approximately by using an adjusted
                # absorption cross section.
                if 'total' in self.mgxs_types or 'transport' in self.mgxs_types:
                    if xsdata.scatter_format == 'legendre':
                        for i in range(len(xsdata.temperatures)):
                            if representation == 'isotropic':
                                xsdata._absorption[i] = \
                                    np.subtract(xsdata._total[i], np.sum(
                                        xsdata._scatter_matrix[i][:, :, 0],
                                        axis=1))
                            elif representation == 'angle':
                                xsdata._absorption[i] = \
                                    np.subtract(xsdata._total[i], np.sum(
                                        xsdata._scatter_matrix[i][:, :, :, :, 0],
                                        axis=3))
                    elif xsdata.scatter_format == 'histogram':
                        for i in range(len(xsdata.temperatures)):
                            if representation == 'isotropic':
                                xsdata._absorption[i] = \
                                    np.subtract(xsdata._total[i], np.sum(np.sum(
                                        xsdata._scatter_matrix[i][:, :, :],
                                        axis=2), axis=1))
                            elif representation == 'angle':
                                xsdata._absorption[i] = \
                                    np.subtract(xsdata._total[i], np.sum(np.sum(
                                        xsdata._scatter_matrix[i][:, :, :, :, :],
                                        axis=4), axis=3))
            # if only scatter matrices have been tallied, multiplicity cannot
            # be accounted for
            else:
                msg = 'Scatter multiplicity (such as (n,xn) reactions) '\
                      'are ignored since multiplicity or nu-scatter matrices '\
                      'were not tallied for ' + xsdata_name
                warn(msg, RuntimeWarning)

                if 'scatter matrix' in self.mgxs_types:
                    scatt_mgxs = self.get_mgxs(domain, 'scatter matrix')
                elif 'consistent scatter matrix' in self.mgxs_types:
                    scatt_mgxs = self.get_mgxs(domain, 'consistent scatter matrix')
                else:
                    raise ValueError(f'No scatter matrix found for {xsdata_name}.')

                xsdata.set_scatter_matrix_mgxs(scatt_mgxs, temperature=temperature,
                                               xs_type=xs_type,
                                               nuclide=[nuclide],
                                               subdomain=subdomain)

        # Apply the stored transport correction ratios to the total cross
        # section and the scattering matrix diagonal
        self._apply_transport_correction_ratios(
            xsdata, domain, temperature, nuclide=nuclide, xs_type=xs_type,
            subdomain=subdomain)

        return xsdata

    def create_mg_library(self, xs_type='macro', xsdata_names=None,
                          apply_domain_chi=False, temperature=ROOM_TEMPERATURE_KELVIN):
        """Creates an openmc.MGXSLibrary object to contain the MGXS data for the
        Multi-Group mode of OpenMC.

        Note that this library will not make use of nested temperature tables.
        Every dataset in the library will be treated as if it was at the same
        default temperature.

        Parameters
        ----------
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'. If the Library object is not tallied by
            nuclide this will be set to 'macro' regardless.
        xsdata_names : Iterable of str
            List of names to apply to the "xsdata" entries in the
            resultant mgxs data file. Defaults to 'set1', 'set2', ...
        apply_domain_chi : bool
            This parameter sets whether (True) or not (False) the
            domain-averaged values of chi, chi-prompt, and chi-delayed are to
            be applied to each of the nuclide-dependent fission energy spectra
            of a domain. In effect, if this is True, then every nuclide in the
            domain receives the same flux-weighted Chi. This is useful for
            downstream multigroup solvers that precompute a material-specific
            chi before the transport solve provides group-wise fluxes. Defaults
            to False.
        temperature : float, optional
            The temperature to set in the MGXSLibrary object. Defaults to 294 K
            (room temperature).

        Returns
        -------
        mgxs_file : openmc.MGXSLibrary
            Multi-Group Cross Section File that is ready to be printed to the
            file of choice by the user.

        Raises
        ------
        ValueError
            When the Library object is initialized with insufficient types of
            cross sections for the Library.

        See also
        --------
        Library.dump_to_file()
        Library.create_mg_mode()

        """

        # Check to ensure the Library contains the correct
        # multi-group cross section types
        self.check_library_for_openmc_mgxs()

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])
        if xsdata_names is not None:
            cv.check_iterable_type('xsdata_names', xsdata_names, str)

        # If gathering material-specific data, set the xs_type to macro
        if not self.by_nuclide:
            xs_type = 'macro'

        # Initialize file
        mgxs_file = openmc.MGXSLibrary(
            self.energy_groups, num_delayed_groups=self.num_delayed_groups)

        if self.domain_type == 'mesh':
            # Create the xsdata objects and add to the mgxs_file
            i = 0
            for domain in self.domains:
                if self.by_nuclide:
                    raise NotImplementedError("Mesh domains do not currently "
                                              "support nuclidic tallies")
                for subdomain in domain.indices:
                    # Build & add metadata to XSdata object
                    if xsdata_names is None:
                        xsdata_name = 'set' + str(i + 1)
                    else:
                        xsdata_name = xsdata_names[i]

                    # Create XSdata and Macroscopic for this domain
                    xsdata = self.get_xsdata(domain, xsdata_name,
                                             subdomain=subdomain)
                    mgxs_file.add_xsdata(xsdata)
                    i += 1

        else:
            # Create the xsdata object and add it to the mgxs_file
            for i, domain in enumerate(self.domains):
                if self.by_nuclide:
                    nuclides = domain.get_nuclides()
                else:
                    nuclides = ['total']
                for nuclide in nuclides:
                    # Build & add metadata to XSdata object
                    if xsdata_names is None:
                        xsdata_name = 'set' + str(i + 1)
                    else:
                        xsdata_name = xsdata_names[i]

                    xsdata = self.get_xsdata(domain, xsdata_name,
                                             nuclide=nuclide, xs_type=xs_type,
                                             apply_domain_chi=apply_domain_chi)

                    mgxs_file.add_xsdata(xsdata)

        return mgxs_file

    def create_mg_mode(self, xsdata_names=None, bc=['reflective'] * 6,
                       apply_domain_chi=False):
        """Creates an openmc.MGXSLibrary object to contain the MGXS data for the
        Multi-Group mode of OpenMC as well as the associated openmc.Materials
        and openmc.Geometry objects.

        The created Geometry is the same as that used to generate the MGXS
        data, with the only differences being modifications to point to
        newly-created Materials which point to the multi-group data. This
        method only creates a macroscopic MGXS Library even if nuclidic tallies
        are specified in the Library. Note that this library will not make
        use of nested temperature tables. Every dataset in the library will
        be treated as if it was at the same default temperature.

        Parameters
        ----------
        xsdata_names : Iterable of str
            List of names to apply to the "xsdata" entries in the
            resultant mgxs data file. Defaults to 'set1', 'set2', ...
        bc : iterable of {'reflective', 'periodic', 'transmission', or 'vacuum'}
            Boundary conditions for each of the four faces of a rectangle
            (if applying to a 2D mesh) or six faces of a parallelepiped
            (if applying to a 3D mesh) provided in the following order:
            [x min, x max, y min, y max, z min, z max].  2-D cells do not
            contain the z min and z max entries.
        apply_domain_chi : bool
            This parameter sets whether (True) or not (False) the
            domain-averaged values of chi, chi-prompt, and chi-delayed are to
            be applied to each of the nuclide-dependent fission energy spectra
            of a domain. In effect, if this is True, then every nuclide in the
            domain receives the same flux-weighted Chi. This is useful for
            downstream multigroup solvers that precompute a material-specific
            chi before the transport solve provides group-wise fluxes. Defaults
            to False.

        Returns
        -------
        mgxs_file : openmc.MGXSLibrary
            Multi-Group Cross Section File that is ready to be printed to the
            file of choice by the user.
        materials : openmc.Materials
            Materials file ready to be printed with all the macroscopic data
            present within this Library.
        geometry : openmc.Geometry
            Materials file ready to be printed with all the macroscopic data
            present within this Library.

        Raises
        ------
        ValueError
            When the Library object is initialized with insufficient types of
            cross sections for the Library.

        See also
        --------
        Library.create_mg_library()
        Library.dump_to_file()

        """

        # Check to ensure the Library contains the correct
        # multi-group cross section types
        self.check_library_for_openmc_mgxs()

        # If the domain type is a mesh, then there can only be one domain for
        # this method. This is because we can build a model automatically if
        # the user provided multiple mesh domains for library generation since
        # the multiple meshes could be overlapping or in disparate regions
        # of the continuous energy model. The next step makes sure there is
        # only one before continuing.
        if self.domain_type == 'mesh':
            cv.check_length("domains", self.domains, 1, 1)

        # Get the MGXS File Data
        mgxs_file = self.create_mg_library('macro', xsdata_names,
                                           apply_domain_chi=apply_domain_chi)

        # Now move on the creating the geometry and assigning materials
        if self.domain_type == 'mesh':
            root = openmc.Universe(name='root', universe_id=0)

            # Add cells representative of the mesh with reflective BC
            root_cell, cells = \
                self.domains[0].build_cells(bc)
            root.add_cell(root_cell)

            geometry = openmc.Geometry()
            geometry.root_universe = root
            materials = openmc.Materials()

            for i, subdomain in enumerate(self.domains[0].indices):
                xsdata = mgxs_file.xsdatas[i]

                # Build the macroscopic and assign it to the cell of
                # interest
                macroscopic = openmc.Macroscopic(name=xsdata.name)

                # Create Material and add to collection
                material = openmc.Material(name=xsdata.name)
                material.add_macroscopic(macroscopic)
                materials.append(material)

                # Set the materials for each of the universes
                cells[i].fill = materials[i]

        else:
            # Create a copy of the Geometry for these Macroscopics
            geometry = copy.deepcopy(self.geometry)
            materials = openmc.Materials()

            # Get all Cells from the Geometry for differentiation
            all_cells = geometry.get_all_cells().values()

            # Create the xsdata object and add it to the mgxs_file
            for i, domain in enumerate(self.domains):
                xsdata = mgxs_file.xsdatas[i]

                macroscopic = openmc.Macroscopic(name=xsdata.name)

                # Create Material and add to collection
                material = openmc.Material(name=xsdata.name)
                material.add_macroscopic(macroscopic)
                materials.append(material)

                # Differentiate Geometry with new Material
                if self.domain_type == 'material':
                    # Fill all appropriate Cells with new Material
                    for cell in all_cells:
                        if isinstance(cell.fill, openmc.Material) and cell.fill.id == domain.id:
                            cell.fill = material

                elif self.domain_type == 'cell':
                    for cell in all_cells:
                        if not isinstance(cell.fill, openmc.Material):
                            warn('If the library domain includes a lattice or universe cell '
                            'in conjunction with a consituent cell of that lattice/universe, '
                            'the multi-group simulation will fail')
                        if cell.id == domain.id:
                            cell.fill = material

        return mgxs_file, materials, geometry

    def check_library_for_openmc_mgxs(self):
        """This routine will check the MGXS Types within a Library
        to ensure the MGXS types provided can be used to create
        a MGXS Library for OpenMC's Multi-Group mode.

        The rules to check include:

        - Either total or transport must be present.

          - Both can be available if one wants, but we should
            use whatever corresponds to Library.correction (if P0: transport)

        - Absorption is required.
        - A nu-fission cross section and chi values are not required as a
          fixed source problem could be the target.
        - Fission and kappa-fission are not required as they are only
          needed to support tallies the user may wish to request.
        - Scattering multiplicity should have been tallied for increased model
          accuracy, either using a multiplicity or scatter and nu-scatter matrix
          tally.

        See also
        --------
        Library.create_mg_library()
        Library.create_mg_mode()

        """

        error_flag = False

        # if correction is 'P0', then transport must be provided
        # otherwise total must be provided
        if self.correction == 'P0':
            if ('transport' not in self.mgxs_types and
                'nu-transport' not in self.mgxs_types):
                error_flag = True
                warn('If the "correction" parameter is "P0", then a '
                     '"transport" or "nu-transport" MGXS type is required.')
        else:
            if 'total' not in self.mgxs_types:
                error_flag = True
                warn('If the "correction" parameter is None, then a '
                     '"total" MGXS type is required.')

        # Check consistency of "nu-transport" and "nu-scatter"
        if 'nu-transport' in self.mgxs_types:
            if not ('nu-scatter matrix' in self.mgxs_types or
                    'consistent nu-scatter matrix' in self.mgxs_types):
                error_flag = True
                warn('If a "nu-transport" MGXS type is used then a '
                     '"nu-scatter matrix" or "consistent nu-scatter matrix" '
                     'must also be used.')
        elif 'transport' in self.mgxs_types:
            if not ('scatter matrix' in self.mgxs_types or
                    'consistent scatter matrix' in self.mgxs_types):
                error_flag = True
                warn('If a "transport" MGXS type is used then a '
                     '"scatter matrix" or "consistent scatter matrix" '
                     'must also be used.')

        # Make sure there is some kind of a scattering matrix data
        if 'nu-scatter matrix' not in self.mgxs_types and \
            'consistent nu-scatter matrix' not in self.mgxs_types and \
            'scatter matrix' not in self.mgxs_types and \
            'consistent scatter matrix' not in self.mgxs_types:
            error_flag = True
            warn('A "nu-scatter matrix", "consistent nu-scatter matrix", '
                 '"scatter matrix", or "consistent scatter matrix" MGXS '
                 'type is required.')

        # Make sure there is some kind of a scattering multiplicity matrix data
        if 'multiplicity matrix' not in self.mgxs_types and \
            ('scatter matrix' not in self.mgxs_types or
             'nu-scatter matrix' not in self.mgxs_types) and\
            ('consistent scatter matrix' not in self.mgxs_types or
             'consistent nu-scatter matrix' not in self.mgxs_types):
            warn('A "multiplicity matrix" or both a "scatter" and "nu-scatter" '
                 'matrix MGXS type(s) should be provided.')

        # Ensure absorption is present
        if 'absorption' not in self.mgxs_types:
            error_flag = True
            warn('An "absorption" MGXS type is required but not provided.')

        # Validate the transport correction ratios that will be applied to the
        # exported data. The plain total cross section they require is already
        # guaranteed by the total/transport checks above (a 'total' MGXS when
        # correction is None, a 'transport'/'nu-transport' MGXS otherwise), so
        # only the ratio-specific requirements are checked here.
        if self._transport_correction_ratios:
            if self.scatter_format != 'legendre':
                error_flag = True
                warn('transport_correction_ratios require a "legendre" '
                     'scatter_format.')

            domain_ratios = \
                self._transport_correction_ratios.get(self.domain_type)
            if not domain_ratios:
                warn('The transport_correction_ratios do not contain any '
                     f'entries for the "{self.domain_type}" domain type, so '
                     'no transport correction will be applied.')
            else:
                domain_ids = [domain.id for domain in self.domains]
                for domain_id, ratios in domain_ratios.items():
                    if domain_id not in domain_ids:
                        warn(f'Domain {domain_id} in '
                             'transport_correction_ratios is not in the '
                             'Library and will be ignored.')
                    elif len(ratios) != self.num_groups:
                        error_flag = True
                        warn(f'Expected {self.num_groups} transport '
                             f'correction ratios for domain {domain_id} but '
                             f'got {len(ratios)}.')

        if error_flag:
            raise ValueError('Invalid MGXS configuration encountered.')
