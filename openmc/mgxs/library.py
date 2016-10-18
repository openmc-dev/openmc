import sys
import os
import copy
import pickle
from numbers import Integral
from collections import OrderedDict
from warnings import warn

import numpy as np

import openmc
import openmc.mgxs
import openmc.checkvalue as cv
from openmc.tallies import ESTIMATOR_TYPES


if sys.version_info[0] >= 3:
    basestring = str


class Library(object):
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
    openmc_geometry : openmc.Geometry
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
    openmc_geometry : openmc.Geometry
        An geometry which has been initialized with a root universe
    opencg_geometry : opencg.Geometry
        An OpenCG geometry object equivalent to the OpenMC geometry
        encapsulated by the summary file. Use of this attribute requires
        installation of the OpenCG Python module.
    by_nuclide : bool
        If true, computes cross sections for each nuclide in each domain
    mgxs_types : Iterable of str
        The types of cross sections in the library (e.g., ['total', 'scatter'])
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    domains : Iterable of openmc.Material, openmc.Cell, openmc.Universe or openmc.Mesh
        The spatial domain(s) for which MGXS in the Library are computed
    correction : {'P0', None}
        Apply the P0 correction to scattering matrices if set to 'P0'
    legendre_order : int
        The highest legendre moment in the scattering matrices (default is 0)
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    delayed_groups : list of int
        Delayed groups to filter out the xs
    estimator : str or None
        The tally estimator used to compute multi-group cross sections. If None,
        the default for each MGXS type is used.
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    all_mgxs : collections.OrderedDict
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

    def __init__(self, openmc_geometry, by_nuclide=False,
                 mgxs_types=None, name=''):

        self._name = ''
        self._openmc_geometry = None
        self._opencg_geometry = None
        self._by_nuclide = None
        self._mgxs_types = []
        self._domain_type = None
        self._domains = 'all'
        self._energy_groups = None
        self._delayed_groups = None
        self._correction = 'P0'
        self._legendre_order = 0
        self._tally_trigger = None
        self._all_mgxs = OrderedDict()
        self._sp_filename = None
        self._keff = None
        self._sparse = False
        self._estimator = None

        self.name = name
        self.openmc_geometry = openmc_geometry
        self.by_nuclide = by_nuclide

        if mgxs_types is not None:
            self.mgxs_types = mgxs_types

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self.name
            clone._openmc_geometry = self.openmc_geometry
            clone._opencg_geometry = None
            clone._by_nuclide = self.by_nuclide
            clone._mgxs_types = self.mgxs_types
            clone._domain_type = self.domain_type
            clone._domains = copy.deepcopy(self.domains)
            clone._correction = self.correction
            clone._legendre_order = self.legendre_order
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._delayed_groups = copy.deepcopy(self.delayed_groups, memo)
            clone._tally_trigger = copy.deepcopy(self.tally_trigger, memo)
            clone._all_mgxs = copy.deepcopy(self.all_mgxs)
            clone._sp_filename = self._sp_filename
            clone._keff = self._keff
            clone._sparse = self.sparse

            clone._all_mgxs = OrderedDict()
            for domain in self.domains:
                clone.all_mgxs[domain.id] = OrderedDict()
                for mgxs_type in self.mgxs_types:
                    mgxs = copy.deepcopy(self.all_mgxs[domain.id][mgxs_type])
                    clone.all_mgxs[domain.id][mgxs_type] = mgxs

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def openmc_geometry(self):
        return self._openmc_geometry

    @property
    def opencg_geometry(self):
        if self._opencg_geometry is None:
            from openmc.opencg_compatible import get_opencg_geometry
            self._opencg_geometry = get_opencg_geometry(self._openmc_geometry)
        return self._opencg_geometry

    @property
    def name(self):
        return self._name

    @property
    def mgxs_types(self):
        return self._mgxs_types

    @property
    def by_nuclide(self):
        return self._by_nuclide

    @property
    def domain_type(self):
        return self._domain_type

    @property
    def domains(self):
        if self._domains == 'all':
            if self.domain_type == 'material':
                return self.openmc_geometry.get_all_materials()
            elif self.domain_type in ['cell', 'distribcell']:
                return self.openmc_geometry.get_all_material_cells()
            elif self.domain_type == 'universe':
                return self.openmc_geometry.get_all_universes()
            elif self.domain_type == 'mesh':
                raise ValueError('Unable to get domains for Mesh domain type')
            else:
                raise ValueError('Unable to get domains without a domain type')
        else:
            return self._domains

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def delayed_groups(self):
        return self._delayed_groups

    @property
    def correction(self):
        return self._correction

    @property
    def legendre_order(self):
        return self._legendre_order

    @property
    def tally_trigger(self):
        return self._tally_trigger

    @property
    def estimator(self):
        return self._estimator

    @property
    def num_groups(self):
        return self.energy_groups.num_groups

    @property
    def num_delayed_groups(self):
        if self.delayed_groups == None:
            return 0
        else:
            return len(self.delayed_groups)

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

    @openmc_geometry.setter
    def openmc_geometry(self, openmc_geometry):
        cv.check_type('openmc_geometry', openmc_geometry, openmc.Geometry)
        self._openmc_geometry = openmc_geometry
        self._opencg_geometry = None

    @name.setter
    def name(self, name):
        cv.check_type('name', name, basestring)
        self._name = name

    @mgxs_types.setter
    def mgxs_types(self, mgxs_types):
        all_mgxs_types = openmc.mgxs.MGXS_TYPES + openmc.mgxs.MDGXS_TYPES
        if mgxs_types == 'all':
            self._mgxs_types = all_mgxs_types
        else:
            cv.check_iterable_type('mgxs_types', mgxs_types, basestring)
            for mgxs_type in mgxs_types:
                cv.check_value('mgxs_type', mgxs_type, all_mgxs_types)
            self._mgxs_types = mgxs_types

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)

        if by_nuclide == True and self.domain_type == 'mesh':
            raise ValueError('Unable to create MGXS library by nuclide with '
                             'mesh domain')

        self._by_nuclide = by_nuclide

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, openmc.mgxs.DOMAIN_TYPES)

        if self.by_nuclide == True and domain_type == 'mesh':
            raise ValueError('Unable to create MGXS library by nuclide with '
                             'mesh domain')

        self._domain_type = domain_type

    @domains.setter
    def domains(self, domains):

        # Use all materials, cells or universes in the geometry as domains
        if domains == 'all':
            self._domains = domains

        # User specified a list of material, cell or universe domains
        else:
            if self.domain_type == 'material':
                cv.check_iterable_type('domain', domains, openmc.Material)
                all_domains = self.openmc_geometry.get_all_materials()
            elif self.domain_type in ['cell', 'distribcell']:
                cv.check_iterable_type('domain', domains, openmc.Cell)
                all_domains = self.openmc_geometry.get_all_material_cells()
            elif self.domain_type == 'universe':
                cv.check_iterable_type('domain', domains, openmc.Universe)
                all_domains = self.openmc_geometry.get_all_universes()
            elif self.domain_type == 'mesh':
                cv.check_iterable_type('domain', domains, openmc.Mesh)

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

            self._domains = domains

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    @delayed_groups.setter
    def delayed_groups(self, delayed_groups):

        if delayed_groups != None:

            cv.check_type('delayed groups', delayed_groups, list, int)
            cv.check_greater_than('num delayed groups', len(delayed_groups), 0)

            # Check that the groups are within [1, MAX_DELAYED_GROUPS]
            for group in delayed_groups:
                cv.check_greater_than('delayed group', group, 0)
                cv.check_less_than('delayed group', group,
                                   openmc.mgxs.MAX_DELAYED_GROUPS,
                                   equality=True)

            self._delayed_groups = delayed_groups

    @correction.setter
    def correction(self, correction):
        cv.check_value('correction', correction, ('P0', None))

        if correction == 'P0' and self.legendre_order > 0:
            warn('The P0 correction will be ignored since the scattering '
                 'order "{}" is greater than zero'.format(self.legendre_order))

        self._correction = correction

    @legendre_order.setter
    def legendre_order(self, legendre_order):
        cv.check_type('legendre_order', legendre_order, Integral)
        cv.check_greater_than('legendre_order', legendre_order, 0, equality=True)
        cv.check_less_than('legendre_order', legendre_order, 10, equality=True)

        if self.correction == 'P0' and legendre_order > 0:
            msg = 'The P0 correction will be ignored since the scattering ' \
                  'order {} is greater than zero'.format(self.legendre_order)
            warn(msg, RuntimeWarning)
            self.correction = None

        self._legendre_order = legendre_order

    @tally_trigger.setter
    def tally_trigger(self, tally_trigger):
        cv.check_type('tally trigger', tally_trigger, openmc.Trigger)
        self._tally_trigger = tally_trigger

    @estimator.setter
    def estimator(self, estimator):
        cv.check_value('estimator', estimator, ESTIMATOR_TYPES)
        self._estimator = estimator

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
            self.all_mgxs[domain.id] = OrderedDict()
            for mgxs_type in self.mgxs_types:
                if mgxs_type in openmc.mgxs.MDGXS_TYPES:
                    mgxs = openmc.mgxs.MDGXS.get_mgxs(mgxs_type, name=self.name)
                else:
                    mgxs = openmc.mgxs.MGXS.get_mgxs(mgxs_type, name=self.name)

                mgxs.domain = domain
                mgxs.domain_type = self.domain_type
                mgxs.energy_groups = self.energy_groups
                mgxs.by_nuclide = self.by_nuclide
                if self.estimator is not None:
                    mgxs.estimator = self.estimator

                if mgxs_type in openmc.mgxs.MDGXS_TYPES:
                    mgxs.delayed_groups = self.delayed_groups

                # If a tally trigger was specified, add it to the MGXS
                if self.tally_trigger is not None:
                    mgxs.tally_trigger = self.tally_trigger

                # Specify whether to use a transport ('P0') correction
                if isinstance(mgxs, openmc.mgxs.ScatterMatrixXS):
                    mgxs.correction = self.correction
                    mgxs.legendre_order = self.legendre_order

                self.all_mgxs[domain.id][mgxs_type] = mgxs

    def add_to_tallies_file(self, tallies_file, merge=True):
        """Add all tallies from all MGXS objects to a tallies file.

        NOTE: This assumes that :meth:`Library.build_library` has been called

        Parameters
        ----------
        tallies_file : openmc.Tallies
            A Tallies collection to add each MGXS' tallies to generate a
            'tallies.xml' input file for OpenMC
        merge : bool
            Indicate whether tallies should be merged when possible. Defaults
            to True.

        """

        cv.check_type('tallies_file', tallies_file, openmc.Tallies)

        # Add tallies from each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                for tally in mgxs.tallies.values():
                    tallies_file.append(tally, merge=merge)

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
        self._openmc_geometry = statepoint.summary.openmc_geometry
        self._nuclides = statepoint.summary.nuclides

        if statepoint.run_mode == 'k-eigenvalue':
            self._keff = statepoint.k_combined[0]

        # Load tallies for each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                mgxs.load_from_statepoint(statepoint)
                mgxs.sparse = self.sparse

    def get_mgxs(self, domain, mgxs_type):
        """Return the MGXS object for some domain and reaction rate type.

        This routine searches the library for an MGXS object for the spatial
        domain and reaction rate type requested by the user.

        NOTE: This routine must be called after the build_library() routine.

        Parameters
        ----------
        domain : Material or Cell or Universe or Integral
            The material, cell, or universe object of interest (or its ID)
        mgxs_type : {'total', 'transport', 'nu-transport', 'absorption', 'capture', 'fission', 'nu-fission', 'kappa-fission', 'scatter', 'nu-scatter', 'scatter matrix', 'nu-scatter matrix', 'multiplicity matrix', 'nu-fission matrix', chi', 'chi-prompt', 'inverse-velocity', 'prompt-nu-fission', 'delayed-nu-fission', 'chi-delayed', 'beta'}
            The type of multi-group cross section object to return

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
            cv.check_type('domain', domain, (openmc.Mesh, Integral))

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
            msg = 'Unable to find MGXS type "{0}"'.format(mgxs_type)
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
        Library
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
                         row_column='inout'):
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

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        import h5py

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Add an attribute for the number of energy groups to the HDF5 file
        full_filename = os.path.join(directory, filename)
        full_filename = full_filename.replace(' ', '-')
        f = h5py.File(full_filename, 'w')
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

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_filename = os.path.join(directory, filename + '.pkl')
        full_filename = full_filename.replace(' ', '-')

        # Load and return pickled Library object
        pickle.dump(self, open(full_filename, 'wb'))

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
        Library
            A Library object loaded from the pickle binary file

        See also
        --------
        Library.dump_to_file(mgxs_lib, filename, directory)

        """

        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_filename = os.path.join(directory, filename + '.pkl')
        full_filename = full_filename.replace(' ', '-')

        # Load and return pickled Library object
        return pickle.load(open(full_filename, 'rb'))

    def get_xsdata(self, domain, xsdata_name, nuclide='total', xs_type='macro',
                   order=None, subdomain=None):
        """Generates an openmc.XSdata object describing a multi-group cross section
        dataset for writing to an openmc.MGXSLibrary object.

        Note that this method does not build an XSdata
        object with nested temperature tables.  The temperature of each
        XSdata object will be left at the default value of 300K.

        Parameters
        ----------
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
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
        order : int
            Scattering order for this data entry.  Default is None,
            which will set the XSdata object to use the order of the
            Library.
        subdomain : iterable of int
            This parameter is not used unless using a mesh domain. In that
            case, the subdomain is an [i,j,k] index (1-based indexing) of the
            mesh cell of interest in the openmc.Mesh object.  Note:
            this parameter currently only supports subdomains within a mesh,
            and not the subdomains of a distribcell.

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
                                         openmc.Universe, openmc.Mesh))
        cv.check_type('xsdata_name', xsdata_name, basestring)
        cv.check_type('nuclide', nuclide, basestring)
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])
        cv.check_type('order', order, (type(None), Integral))
        if order is not None:
            cv.check_greater_than('order', order, 0, equality=True)
            cv.check_less_than('order', order, 10, equality=True)
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
        xsdata = openmc.XSdata(name, self.energy_groups)

        if order is None:
            # Set the order to the Library's order (the defualt behavior)
            xsdata.order = self.legendre_order
        else:
            # Set the order of the xsdata object to the minimum of
            # the provided order or the Library's order.
            xsdata.order = min(order, self.legendre_order)

        # Right now only 'legendre' data and isotropic weighting is supported
        self.scatter_format = 'legendre'
        self.representation = 'isotropic'

        if nuclide != 'total':
            xsdata.atomic_weight_ratio = self._nuclides[nuclide][1]

        if subdomain is None:
            subdomain = 'all'
        else:
            subdomain = [subdomain]

        # Now get xs data itself
        if 'nu-transport' in self.mgxs_types and self.correction == 'P0':
            mymgxs = self.get_mgxs(domain, 'nu-transport')
            xsdata.set_total_mgxs(mymgxs, xs_type=xs_type, nuclide=[nuclide],
                                  subdomain=subdomain)
        elif 'total' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'total')
            xsdata.set_total_mgxs(mymgxs, xs_type=xs_type, nuclide=[nuclide],
                                  subdomain=subdomain)
        if 'absorption' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'absorption')
            xsdata.set_absorption_mgxs(mymgxs, xs_type=xs_type,
                                       nuclide=[nuclide],
                                       subdomain=subdomain)
        if 'fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'fission')
            xsdata.set_fission_mgxs(mymgxs, xs_type=xs_type,
                                    nuclide=[nuclide], subdomain=subdomain)
        if 'kappa-fission' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'kappa-fission')
            xsdata.set_kappa_fission_mgxs(mymgxs, xs_type=xs_type,
                                          nuclide=[nuclide],
                                          subdomain=subdomain)
        # For chi and nu-fission we can either have only a nu-fission matrix
        # provided, or vectors of chi and nu-fission provided
        if 'nu-fission matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'nu-fission matrix')
            xsdata.set_nu_fission_mgxs(mymgxs, xs_type=xs_type,
                                       nuclide=[nuclide],
                                       subdomain=subdomain)
        else:
            if 'chi' in self.mgxs_types:
                mymgxs = self.get_mgxs(domain, 'chi')
                xsdata.set_chi_mgxs(mymgxs, xs_type=xs_type, nuclide=[nuclide],
                                    subdomain=subdomain)
            if 'nu-fission' in self.mgxs_types:
                mymgxs = self.get_mgxs(domain, 'nu-fission')
                xsdata.set_nu_fission_mgxs(mymgxs, xs_type=xs_type,
                                           nuclide=[nuclide],
                                           subdomain=subdomain)
        # If multiplicity matrix is available, prefer that
        if 'multiplicity matrix' in self.mgxs_types:
            mymgxs = self.get_mgxs(domain, 'multiplicity matrix')
            xsdata.set_multiplicity_matrix_mgxs(mymgxs, xs_type=xs_type,
                                                nuclide=[nuclide],
                                                subdomain=subdomain)
            using_multiplicity = True
        # multiplicity will fall back to using scatter and nu-scatter
        elif ((('scatter matrix' in self.mgxs_types) and
               ('nu-scatter matrix' in self.mgxs_types))):
            scatt_mgxs = self.get_mgxs(domain, 'scatter matrix')
            nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
            xsdata.set_multiplicity_matrix_mgxs(nuscatt_mgxs, scatt_mgxs,
                                                xs_type=xs_type,
                                                nuclide=[nuclide],
                                                subdomain=subdomain)
            using_multiplicity = True
        else:
            using_multiplicity = False

        if using_multiplicity:
            nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
            xsdata.set_scatter_matrix_mgxs(nuscatt_mgxs, xs_type=xs_type,
                                           nuclide=[nuclide],
                                           subdomain=subdomain)
        else:
            if 'nu-scatter matrix' in self.mgxs_types:
                nuscatt_mgxs = self.get_mgxs(domain, 'nu-scatter matrix')
                xsdata.set_scatter_matrix_mgxs(nuscatt_mgxs, xs_type=xs_type,
                                               nuclide=[nuclide],
                                               subdomain=subdomain)

                # Since we are not using multiplicity, then
                # scattering multiplication (nu-scatter) must be
                # accounted for approximately by using an adjusted
                # absorption cross section.
                if 'total' in self.mgxs_types or 'transport' in self.mgxs_types:
                    for i in range(len(xsdata.temperatures)):
                        xsdata._absorption[i] = \
                            np.subtract(xsdata._total[i], np.sum(
                                xsdata._scatter_matrix[i][0, :, :], axis=1))

        return xsdata

    def create_mg_library(self, xs_type='macro', xsdata_names=None):
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
            cv.check_iterable_type('xsdata_names', xsdata_names, basestring)

        # If gathering material-specific data, set the xs_type to macro
        if not self.by_nuclide:
            xs_type = 'macro'

        # Initialize file
        mgxs_file = openmc.MGXSLibrary(self.energy_groups)

        if self.domain_type == 'mesh':
            # Create the xsdata objects and add to the mgxs_file
            i = 0
            for domain in self.domains:
                if self.by_nuclide:
                    raise NotImplementedError("Mesh domains do not currently "
                                              "support nuclidic tallies")
                for subdomain in domain.cell_generator():
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
                    if nuclide != 'total':
                        xsdata_name += '_' + nuclide

                    xsdata = self.get_xsdata(domain, xsdata_name,
                                             nuclide=nuclide, xs_type=xs_type)

                    mgxs_file.add_xsdata(xsdata)

        return mgxs_file

    def create_mg_mode(self, xsdata_names=None, bc=['reflective'] * 6):
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
        mgxs_file = self.create_mg_library('macro', xsdata_names)

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

            for i, subdomain in enumerate(self.domains[0].cell_generator()):
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
            geometry = copy.deepcopy(self.openmc_geometry)
            materials = openmc.Materials()

            # Get all Cells from the Geometry for differentiation
            all_cells = geometry.get_all_material_cells()

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
                        if cell.fill.id == domain.id:
                            cell.fill = material

                elif self.domain_type == 'cell':
                    for cell in all_cells:
                        if cell.id == domain.id:
                            cell.fill = material

        return mgxs_file, materials, geometry

    def check_library_for_openmc_mgxs(self):
        """This routine will check the MGXS Types within a Library
        to ensure the MGXS types provided can be used to create
        a MGXS Library for OpenMC's Multi-Group mode.

        The rules to check include:

        - Either total or transport should be present.

          - Both can be available if one wants, but we should
            use whatever corresponds to Library.correction (if P0: transport)

        - Absorption and total (or transport) are required.
        - A nu-fission cross section and chi values are not required as a
          fixed source problem could be the target.
        - Fission and kappa-fission are not required as they are only
          needed to support tallies the user may wish to request.
        - A nu-scatter matrix is required.

          - Having a multiplicity matrix is preferred.
          - Having both nu-scatter (of any order) and scatter
            (at least isotropic) matrices is the second choice.
          - If only nu-scatter, need total (not transport), to
            be used in adjusting absorption
            (i.e., reduced_abs = tot - nuscatt)

        See also
        --------
        Library.create_mg_library()
        Library.create_mg_mode()

        """

        error_flag = False
        # Ensure absorption is present
        if 'absorption' not in self.mgxs_types:
            error_flag = True
            warn('An "absorption" MGXS type is required but not provided.')
        # Ensure nu-scattering matrix is required
        if 'nu-scatter matrix' not in self.mgxs_types:
            error_flag = True
            warn('A "nu-scatter matrix" MGXS type is required but not provided.')
        else:
            # Ok, now see the status of scatter and/or multiplicity
            if ((('scatter matrix' not in self.mgxs_types) and
                 ('multiplicity matrix' not in self.mgxs_types))):
                # We dont have data needed for multiplicity matrix, therefore
                # we need total, and not transport.
                if 'total' not in self.mgxs_types:
                    error_flag = True
                    warn('A "total" MGXS type is required if a '
                         'scattering matrix is not provided.')
        # Total or transport can be present, but if using
        # self.correction=="P0", then we should use transport.
        if (((self.correction is "P0") and
             ('nu-transport' not in self.mgxs_types))):
            error_flag = True
            warn('A "nu-transport" MGXS type is required since a "P0" '
                 'correction is applied, but a "nu-transport" MGXS is '
                 'not provided.')
        elif (((self.correction is None) and
               ('total' not in self.mgxs_types))):
            error_flag = True
            warn('A "total" MGXS type is required, but not provided.')

        if error_flag:
            raise ValueError('Invalid MGXS configuration encountered.')
