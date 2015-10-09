import sys
import copy
from numbers import Integral

import openmc
import openmc.mgxs
import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str


class Library(object):

    def __init__(self, openmc_geometry, by_nuclide=False,
                 mgxs_types=None, name=''):

        self._name = ''
        self._openmc_geometry = None
        self._by_nuclide = None
        self._mgxs_types = []
        self._domain_type = None
        self._energy_groups = None
        self._all_mgxs = {}

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
            clone._by_nuclide = self.by_nuclide
            clone._mgxs_types = self.mgxs_types
            clone._domain_type = self.domain_type
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._all_mgxs = self.all_mgxs

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
    def openmc_geometry(self):
        return self._openmc_geometry

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
    def domains(self):
        if self.domain_type is None:
            raise ValueError('Unable to get all domains without a domain type')

        if self.domain_type == 'material':
            return self.openmc_geometry.get_all_materials()
        elif self.domain_type == 'cell' or self.domain_type == 'distribcell':
            return self.openmc_geometry.get_all_material_cells()
        elif self.domain_type == 'universe':
            return self.openmc_geometry.get_all_universes()

    @property
    def domain_type(self):
        return self._domain_type

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def num_groups(self):
        return self.energy_groups.num_groups

    @property
    def all_mgxs(self):
        return self._all_mgxs

    @openmc_geometry.setter
    def openmc_geometry(self, openmc_geometry):
        cv.check_type('openmc_geometry', openmc_geometry, openmc.Geometry)
        self._openmc_geometry = openmc_geometry

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

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy groups', energy_groups, openmc.mgxs.EnergyGroups)
        self._energy_groups = energy_groups

    def build_library(self):
        """
        """

        # Initialize MGXS for each domain and mgxs type and store in dictionary
        for domain in self.domains:
            self.all_mgxs[domain.id] = {}
            for mgxs_type in self.mgxs_types:
                mgxs = openmc.mgxs.MGXS.get_mgxs(mgxs_type, name=self.name)
                mgxs.domain = domain
                mgxs.domain_type = self.domain_type
                mgxs.energy_groups = self.energy_groups
                mgxs.by_nuclide = self.by_nuclide
                mgxs.create_tallies()
                self.all_mgxs[domain.id][mgxs_type] = mgxs

    def add_to_tallies_file(self, tallies_file):
        """

        NOTE: This assumes that build_library() has been called

        :param tallies_file:
        :return:
        """

        cv.check_type('tallies_file', tallies_file, openmc.TalliesFile)

        # Add tallies from each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                for tally_id, tally in mgxs.tallies.items():
                    tallies_file.add_tally(tally, merge=True)

    def load_from_statepoint(self, statepoint):
        """

        :param statepoint:
        :return:
        """

        cv.check_type('statepoint', statepoint, openmc.StatePoint)

        # Load tallies for each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.get_mgxs(domain, mgxs_type)
                mgxs.load_from_statepoint(statepoint)
                mgxs.compute_xs()

    def get_mgxs(self, domain, mgxs_type):
        """

        :param domain:
        :param mgxs_type:
        :return:
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
                      'library'.format(self.domain_type, domain)
                raise ValueError(msg)
        else:
            domain_id = domain.id

        # Check that requested domain is included in library
        if mgxs_type not in self.mgxs_types:
                msg = 'Unable to find MGXS type "{0}"'.format(mgxs_type)
                raise ValueError(msg)

        return self.all_mgxs[domain_id][mgxs_type]

    def get_condensed_library(self, coarse_groups):
        """

        :param coarse_groups:
        :return:
        """

        if self.energy_groups is None:
            msg = 'Unable to get a condensed coarse group cross section ' \
                  'library since the fine energy groups have not yet been set'
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

    def build_hdf5_store(self, filename='mgxs', directory='mgxs', xs_type='macro'):
        """Export the multi-group cross section library to an HDF5 binary file.

        This method constructs an HDF5 file which stores the multi-group
        cross section data. The data is stored in a hierarchy of HDF5 groups
        from the domain type, domain id, subdomain id (for distribcell domains),
        nuclides and cross section types. Two datasets for the mean and standard
        deviation are stored for each subdomain entry in the HDF5 file.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file (default is 'mgxs')
        directory : str
            Directory for the HDF5 file (default is 'mgxs')
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns

        """

        # Load tallies for each MGXS for each domain and mgxs type
        for domain in self.domains:
            for mgxs_type in self.mgxs_types:
                mgxs = self.all_mgxs[domain.id][mgxs_type]
                mgxs.build_hdf5_store(filename, directory, xs_type)