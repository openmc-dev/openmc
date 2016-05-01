import sys
import os
import copy
import pickle
from numbers import Integral
from collections import OrderedDict
import numpy as np
import warnings.warn as warn

import openmc
import openmc.mgxs
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str


class Library(object):
    """A multi-group cross section library for some energy group structure.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for deterministic neutronics calculations.

    This class helps automate the generation of MGXS objects for some energy
    group structure and domain type. The Library serves as a collection for
    MGXS objects with routines to automate the initialization of tallies for
    input files, the loading of tally data from statepoint files, data storage,
    energy group condensation and more.

    Parameters
    ----------
    openmc_geometry : openmc.Geometry
        An geometry which has been initialized with a root universe
    by_nuclide : bool
        If true, computes cross sections for each nuclide in each domain
    mgxs_types : Iterable of str
        The types of cross sections in the library (e.g., ['total', 'scatter'])
    name : str, optional
        Name of the multi-group cross section. library Used as a label to
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
    domain_type : {'material', 'cell', 'distribcell', 'universe'}
        Domain type for spatial homogenization
    domains : Iterable of openmc.Material, openmc.Cell or openmc.Universe
        The spatial domain(s) for which MGXS in the Library are computed
    correction : {'P0', None}
        Apply the P0 correction to scattering matrices if set to 'P0'
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
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
        self._correction = 'P0'
        self._energy_groups = None
        self._tally_trigger = None
        self._all_mgxs = OrderedDict()
        self._sp_filename = None
        self._keff = None
        self._sparse = False

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
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
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
            else:
                raise ValueError('Unable to get domains without a domain type')
        else:
            return self._domains

    @property
    def correction(self):
        return self._correction

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def tally_trigger(self):
        return self._tally_trigger

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
        if mgxs_types == 'all':
            self._mgxs_types = openmc.mgxs.MGXS_TYPES
        else:
            cv.check_iterable_type('mgxs_types', mgxs_types, basestring)
            for mgxs_type in mgxs_types:
                cv.check_value('mgxs_type', mgxs_type, openmc.mgxs.MGXS_TYPES)
            self._mgxs_types = mgxs_types

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)
        self._by_nuclide = by_nuclide

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, tuple(openmc.mgxs.DOMAIN_TYPES))
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
            else:
                msg = 'Unable to set domains with ' \
                      'domain type "{}"'.format(self.domain_type)
                raise ValueError(msg)

            # Check that each domain can be found in the geometry
            for domain in domains:
                if domain not in all_domains:
                    msg = 'Domain "{}" could not be found in the ' \
                          'geometry.'.format(domain)
                    raise ValueError(msg)

            self._domains = domains

    @correction.setter
    def correction(self, correction):
        cv.check_value('correction', correction, ('P0', None))
        self._correction = correction

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    @tally_trigger.setter
    def tally_trigger(self, tally_trigger):
        cv.check_type('tally trigger', tally_trigger, openmc.Trigger)
        self._tally_trigger = tally_trigger

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
                mgxs = openmc.mgxs.MGXS.get_mgxs(mgxs_type, name=self.name)
                mgxs.domain = domain
                mgxs.domain_type = self.domain_type
                mgxs.energy_groups = self.energy_groups
                mgxs.by_nuclide = self.by_nuclide

                # If a tally trigger was specified, add it to the MGXS
                if self.tally_trigger:
                    mgxs.tally_trigger = self.tally_trigger

                # Specify whether to use a transport ('P0') correction
                if isinstance(mgxs, openmc.mgxs.ScatterMatrixXS):
                    mgxs.correction = self.correction

                self.all_mgxs[domain.id][mgxs_type] = mgxs

    def add_to_tallies_file(self, tallies_file, merge=True):
        """Add all tallies from all MGXS objects to a tallies file.

        NOTE: This assumes that :meth:`Library.build_library` has been called

        Parameters
        ----------
        tallies_file : openmc.Tallies
            A Tallies collection to add each MGXS' tallies to generate a
            "tallies.xml" input file for OpenMC
        merge : bool
            Indicate whether tallies should be merged when possible. Defaults
            to True.

        """

        cv.check_type('tallies_file', tallies_file, openmc.Tallies)

        # Add tallies from each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                for tally_id, tally in mgxs.tallies.items():
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
        mgxs_type : {'total', 'transport', 'absorption', 'capture', 'fission', 'nu-fission', 'kappa-fission', 'scatter', 'nu-scatter', 'scatter matrix', 'nu-scatter matrix', 'chi'}
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

        # Check that requested domain is included in library
        if cv._isinstance(domain, Integral):
            domain_id = domain
            for domain in self.domains:
                if domain_id == domain.id:
                    break
            else:
                msg = 'Unable to find MGXS for {0} "{1}" in ' \
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
                         subdomains='all', nuclides='all', xs_type='macro'):
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
            may be a list of nuclide name strings (e.g., ['U-235', 'U-238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

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
        f.attrs["# groups"] = self.num_groups
        f.close()

        # Export MGXS for each domain and mgxs type to an HDF5 file
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.all_mgxs[domain.id][mgxs_type]

                if subdomains == 'avg':
                    mgxs = mgxs.get_subdomain_avg_xs()

                mgxs.build_hdf5_store(filename, directory,
                                      xs_type=xs_type, nuclides=nuclides)

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

    def write_mg_library(self, xs_type='micro', domain_names=None, xs_ids=None,
                         filename='mg_cross_sections', directory='./'):
        """Create a cross-section data library file for the Multi-Group
        mode of OpenMC.

        Parameters
        ----------
        xs_type: {'macro', 'micro'}
            Provide the macro or micro cross section in units of cm^-1 or
            barns. Defaults to 'macro'. If the Library object is not tallied by
            nuclide this will be set to 'macro' regardless
        domain_names : Iterable of str
            List of names to apply to the xsdata entries in the
            resultant mgxs data file. Defaults to "set1", "set2", ...
        xs_ids : str or Iterable of str
            Cross section set identifier (i.e., "71c") for all
            data sets (if only str) or for each individual one
            (if iterable of str). Defaults to '1g'
        filename : str
            Filename for the pickle file. Defaults to 'mg_cross_sections'.
        directory : str
            Directory for the pickle file. Defaults to './' (the
            current working directory).

        See also
        --------
        Library.dump_to_file(mgxs_lib, filename, directory)

        """

        # Check data types provided

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])
        if not self.by_nuclide:
            xs_type = 'macro'
        if domain_names is not None:
            cv.check_iterable_type('domain_names', filename, basestring)
        if xs_ids is not None:
            if isinstance(xs_ids, basestring):
                # If we only have a string lets convert it now to a list
                # of strings.
                xs_ids = [xs_ids for i in range(len(self.domains))]
            else:
                cv.check_iterable_type('xs_ids', xs_ids, basestring)
        else:
            xs_ids = ['.1g']
        cv.check_type('filename', filename, basestring)
        cv.check_type('directory', directory, basestring)

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        full_filename = os.path.join(directory, filename + '.xml')
        full_filename = full_filename.replace(' ', '-')

        # Initialize file
        mgxs_file = openmc.MGXSLibrary(self.energy_groups)

        # Set the scattering order as isotropic until
        # support  for higher orders are included
        order = 0

        # Build XSdata objects
        xsdatas = []
        for i in range(len(self.domains)):
            id = self.domains[i].id
            if not self.by_nuclide:
                # Use k instead of i simply because k will be used for
                # the nuclide index in the else part of this conditional
                # and using k allows us to use the same code.
                k = i
                # Build & add metadata to XSdata object
                # (Use i here because k in nuclides will add chars to this)
                if domain_names is None:
                    name = 'set' + str(i + 1)
                else:
                    name = domain_names[i]
                name += xs_ids[k]
                xsdata = openmc.XSdata(name, self.energy_groups)
                xsdata.order = order

                # Now get xs data itself
                if 'total' in self.mgxs_types:
                    xsdata.set_total(self.all_mgxs[id]['total'],
                                     xs_type=xs_type, subdomains=(k + 1,))
                if 'absorption' in self.mgxs_types:
                    xsdata.set_absorption(self.all_mgxs[id]['absorption'],
                                          xs_type=xs_type, subdomains=(k + 1,))
                if 'fission' in self.mgxs_types:
                    xsdata.set_fission(self.all_mgxs[id]['fission'],
                                       xs_type=xs_type, subdomains=(k + 1,))
                if 'kappa-fission' in self.mgxs_types:
                    xsdata.set_k_fission(self.all_mgxs[id]['kappa-fission'],
                                         xs_type=xs_type, subdomains=(k + 1,))
                if 'chi' in self.mgxs_types:
                    xsdata.set_chi(self.all_mgxs[id]['chi'],
                                   xs_type=xs_type, subdomains=(k + 1,))
                if 'nu-fission' in self.mgxs_types:
                    xsdata.set_nu_fission(self.all_mgxs[id]['nu-fission'],
                                          xs_type=xs_type, subdomains=(k + 1,))
                # multiplicity requires scatter and nu-scatter
                if (('scatter' in self.mgxs_types) and ('nu-scatter' in
                                                        self.mgxs_types)):
                    xsdata.set_multiplicity(self.all_mgxs[id]['nu-scatter'],
                                            self.all_mgxs[id]['scatter'],
                                            xs_type=xs_type,
                                            subdomains=(k + 1,))
                    using_multiplicity = True
                else:
                    using_multiplicity = False

                if using_multiplicity:
                    xsdata.set_scatter(self.all_mgxs[id]['scatter'],
                                       xs_type=xs_type,
                                       subdomains=(k + 1,))
                else:
                    if 'nu-scatter' in self.mgxs_types:
                        xsdata.set_scatter(self.all_mgxs[id]['nu-scatter'],
                                           xs_type=xs_type,
                                           subdomains=(k + 1,))
                        # Since we are not using multiplicity, then
                        # scattering multiplication (nu-scatter) must be
                        # accounted for approximately by using an adjusted
                        # absorption cross section.
                        if self.total is not None:
                            xsdata.absorption = \
                                np.subtract(xsdata.total,
                                            np.sum(xsdata.scatter[0, :, :],
                                                   axis=1))
                        else:
                            # Total isnt included so we cant do the above
                            # approximation w/out changing absorption instead.
                            # That can be done with:
                            # SigA' = SigA - (nuSigS - SigS)
                            # Doing so would mean essentially duplicating
                            # set_scatter from MGXSLibrary to obtain the
                            # SigS, which would be big and ugly once
                            # angle filters are available.
                            # Instead, raise a warning about the
                            # lack of neutron balance and then use scatter
                            # instead of nu-scatter
                            if 'scatter' in self.mgxs_types:
                                msg = "To properly use the 'nu-scatter' " + \
                                      "MGXS type and maintain neutron " + \
                                      "balance, a 'total' MGXS type " + \
                                      "should be provided."
                                warn(msg)
                                xsdata.set_scatter(
                                    self.all_mgxs[id]['scatter'],
                                    xs_type=xs_type, subdomains=(k + 1,))
                            else:
                                # Welp, cant do that either. Quit while ahead.
                                msg = "Total X/S must be provided if using" + \
                                      " nu-scatter as the scattering data"
                                raise ValueError(msg)

                xsdatas.append(xsdata)
            else:
                pass

        # Add XSdatas to file

        # Finally, write the file
        mgxs_file.export_to_xml(full_filename)
