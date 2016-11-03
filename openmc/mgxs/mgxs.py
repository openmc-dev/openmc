from __future__ import division

from collections import OrderedDict
from numbers import Integral
import warnings
import os
import sys
import copy
from abc import ABCMeta
import itertools

from six import add_metaclass, string_types
import numpy as np

import openmc
import openmc.checkvalue as cv
from openmc.tallies import ESTIMATOR_TYPES
from openmc.mgxs import EnergyGroups


# Supported cross section types
MGXS_TYPES = ['total',
              'transport',
              'nu-transport',
              'absorption',
              'capture',
              'fission',
              'nu-fission',
              'kappa-fission',
              'scatter',
              'nu-scatter',
              'scatter matrix',
              'nu-scatter matrix',
              'multiplicity matrix',
              'nu-fission matrix',
              'chi',
              'chi-prompt',
              'inverse-velocity',
              'prompt-nu-fission']

# Supported domain types
DOMAIN_TYPES = ['cell',
                'distribcell',
                'universe',
                'material',
                'mesh']

# Filter types corresponding to each domain
_DOMAIN_TO_FILTER = {'cell': openmc.CellFilter,
                     'distribcell': openmc.DistribcellFilter,
                     'universe': openmc.UniverseFilter,
                     'material': openmc.MaterialFilter,
                     'mesh': openmc.MeshFilter}

# Supported domain classes
_DOMAINS = (openmc.Cell,
            openmc.Universe,
            openmc.Material,
            openmc.Mesh)


@add_metaclass(ABCMeta)
class MGXS(object):
    """An abstract multi-group cross section for some energy group structure
    within some spatial domain.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations.

    NOTE: Users should instantiate the subclasses of this abstract class.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """
    def __init__(self, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name=''):
        self._name = ''
        self._rxn_type = None
        self._by_nuclide = None
        self._nuclides = None
        self._estimator = 'tracklength'
        self._domain = None
        self._domain_type = None
        self._energy_groups = None
        self._tally_trigger = None
        self._tallies = None
        self._rxn_rate_tally = None
        self._xs_tally = None
        self._sparse = False
        self._loaded_sp = False
        self._derived = False
        self._hdf5_key = None
        self._valid_estimators = ESTIMATOR_TYPES

        self.name = name
        self.by_nuclide = by_nuclide

        if domain_type is not None:
            self.domain_type = domain_type
        if domain is not None:
            self.domain = domain
        if energy_groups is not None:
            self.energy_groups = energy_groups

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, copy it
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._name = self.name
            clone._rxn_type = self.rxn_type
            clone._by_nuclide = self.by_nuclide
            clone._nuclides = copy.deepcopy(self._nuclides)
            clone._domain = self.domain
            clone._domain_type = self.domain_type
            clone._energy_groups = copy.deepcopy(self.energy_groups, memo)
            clone._tally_trigger = copy.deepcopy(self.tally_trigger, memo)
            clone._rxn_rate_tally = copy.deepcopy(self._rxn_rate_tally, memo)
            clone._xs_tally = copy.deepcopy(self._xs_tally, memo)
            clone._sparse = self.sparse
            clone._derived = self.derived

            clone._tallies = OrderedDict()
            for tally_type, tally in self.tallies.items():
                clone.tallies[tally_type] = copy.deepcopy(tally, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def name(self):
        return self._name

    @property
    def rxn_type(self):
        return self._rxn_type

    @property
    def by_nuclide(self):
        return self._by_nuclide

    @property
    def domain(self):
        return self._domain

    @property
    def domain_type(self):
        return self._domain_type

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
    def scores(self):
        return ['flux', self.rxn_type]

    @property
    def filters(self):
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.EnergyFilter(group_edges)
        return [[energy_filter]] * len(self.scores)

    @property
    def tally_keys(self):
        return self.scores

    @property
    def estimator(self):
        return self._estimator

    @property
    def tallies(self):

        # Instantiate tallies if they do not exist
        if self._tallies is None:

            # Initialize a collection of Tallies
            self._tallies = OrderedDict()

            # Create a domain Filter object
            filter_type = _DOMAIN_TO_FILTER[self.domain_type]
            if self.domain_type == 'mesh':
                domain_filter = filter_type(self.domain)
            else:
                domain_filter = filter_type(self.domain.id)

            # Create each Tally needed to compute the multi group cross section
            tally_metadata = zip(self.scores, self.tally_keys, self.filters)
            for score, key, filters in tally_metadata:
                self._tallies[key] = openmc.Tally(name=self.name)
                self._tallies[key].scores = [score]
                self._tallies[key].estimator = self.estimator
                self._tallies[key].filters = [domain_filter]

                # If a tally trigger was specified, add it to each tally
                if self.tally_trigger:
                    trigger_clone = copy.deepcopy(self.tally_trigger)
                    trigger_clone.scores = [score]
                    self._tallies[key].triggers.append(trigger_clone)

                # Add non-domain specific Filters (e.g., 'energy') to the Tally
                for add_filter in filters:
                    self._tallies[key].filters.append(add_filter)

                # If this is a by-nuclide cross-section, add nuclides to Tally
                if self.by_nuclide and score != 'flux':
                    self._tallies[key].nuclides += self.get_nuclides()
                else:
                    self._tallies[key].nuclides.append('total')

        return self._tallies

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies[self.rxn_type]
            self._rxn_rate_tally.sparse = self.sparse

        return self._rxn_rate_tally

    @property
    def xs_tally(self):
        if self._xs_tally is None:
            if self.tallies is None:
                msg = 'Unable to get xs_tally since tallies have ' \
                      'not been loaded from a statepoint'
                raise ValueError(msg)

            self._xs_tally = self.rxn_rate_tally / self.tallies['flux']
            self._compute_xs()

        return self._xs_tally

    @property
    def sparse(self):
        return self._sparse

    @property
    def num_subdomains(self):
        if self.domain_type.startswith('avg('):
            domain_type = self.domain_type[4:-1]
        else:
            domain_type = self.domain_type
        filter_type = _DOMAIN_TO_FILTER[domain_type]
        domain_filter = self.xs_tally.find_filter(filter_type)
        return domain_filter.num_bins

    @property
    def num_nuclides(self):
        if self.by_nuclide:
            return len(self.get_nuclides())
        else:
            return 1

    @property
    def nuclides(self):
        if self.by_nuclide:
            return self.get_nuclides()
        else:
            return ['sum']

    @property
    def loaded_sp(self):
        return self._loaded_sp

    @property
    def derived(self):
        return self._derived

    @property
    def hdf5_key(self):
        if self._hdf5_key is not None:
            return self._hdf5_key
        else:
            return self._rxn_type

    @name.setter
    def name(self, name):
        cv.check_type('name', name, string_types)
        self._name = name

    @by_nuclide.setter
    def by_nuclide(self, by_nuclide):
        cv.check_type('by_nuclide', by_nuclide, bool)
        self._by_nuclide = by_nuclide

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_iterable_type('nuclides', nuclides, string_types)
        self._nuclides = nuclides

    @estimator.setter
    def estimator(self, estimator):
        cv.check_value('estimator', estimator, self._valid_estimators)
        self._estimator = estimator

    @domain.setter
    def domain(self, domain):
        cv.check_type('domain', domain, _DOMAINS)
        self._domain = domain

        # Assign a domain type
        if self.domain_type is None:
            if isinstance(domain, openmc.Material):
                self._domain_type = 'material'
            elif isinstance(domain, openmc.Cell):
                self._domain_type = 'cell'
            elif isinstance(domain, openmc.Universe):
                self._domain_type = 'universe'
            elif isinstance(domain, openmc.Mesh):
                self._domain_type = 'mesh'

    @domain_type.setter
    def domain_type(self, domain_type):
        cv.check_value('domain type', domain_type, DOMAIN_TYPES)
        self._domain_type = domain_type

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

        # Sparsify or densify the derived MGXS tallies and the base tallies
        if self._xs_tally:
            self.xs_tally.sparse = sparse
        if self._rxn_rate_tally:
            self.rxn_rate_tally.sparse = sparse

        for tally_name in self.tallies:
            self.tallies[tally_name].sparse = sparse

        self._sparse = sparse

    @staticmethod
    def get_mgxs(mgxs_type, domain=None, domain_type=None,
                 energy_groups=None, by_nuclide=False, name=''):
        """Return a MGXS subclass object for some energy group structure within
        some spatial domain for some reaction type.

        This is a factory method which can be used to quickly create MGXS
        subclass objects for various reaction types.

        Parameters
        ----------
        mgxs_type : {'total', 'transport', 'nu-transport', 'absorption', 'capture', 'fission', 'nu-fission', 'kappa-fission', 'scatter', 'nu-scatter', 'scatter matrix', 'nu-scatter matrix', 'multiplicity matrix', 'nu-fission matrix', 'chi', 'chi-prompt', 'inverse-velocity', 'prompt-nu-fission'}
            The type of multi-group cross section object to return
        domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
            The domain for spatial homogenization
        domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
            The domain type for spatial homogenization
        energy_groups : openmc.mgxs.EnergyGroups
            The energy group structure for energy condensation
        by_nuclide : bool
            If true, computes cross sections for each nuclide in domain.
            Defaults to False
        name : str, optional
            Name of the multi-group cross section. Used as a label to identify
            tallies in OpenMC 'tallies.xml' file. Defaults to the empty string.

        Returns
        -------
        openmc.mgxs.MGXS
            A subclass of the abstract MGXS class for the multi-group cross
            section type requested by the user

        """

        cv.check_value('mgxs_type', mgxs_type, MGXS_TYPES)

        if mgxs_type == 'total':
            mgxs = TotalXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'transport':
            mgxs = TransportXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-transport':
            mgxs = NuTransportXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'absorption':
            mgxs = AbsorptionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'capture':
            mgxs = CaptureXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'fission':
            mgxs = FissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-fission':
            mgxs = NuFissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'kappa-fission':
            mgxs = KappaFissionXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'scatter':
            mgxs = ScatterXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-scatter':
            mgxs = NuScatterXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'scatter matrix':
            mgxs = ScatterMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-scatter matrix':
            mgxs = NuScatterMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'multiplicity matrix':
            mgxs = MultiplicityMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'nu-fission matrix':
            mgxs = NuFissionMatrixXS(domain, domain_type, energy_groups)
        elif mgxs_type == 'chi':
            mgxs = Chi(domain, domain_type, energy_groups)
        elif mgxs_type == 'chi-prompt':
            mgxs = ChiPrompt(domain, domain_type, energy_groups)
        elif mgxs_type == 'inverse-velocity':
            mgxs = InverseVelocity(domain, domain_type, energy_groups)
        elif mgxs_type == 'prompt-nu-fission':
            mgxs = PromptNuFissionXS(domain, domain_type, energy_groups)

        mgxs.by_nuclide = by_nuclide
        mgxs.name = name
        return mgxs

    def get_nuclides(self):
        """Get all nuclides in the cross section's spatial domain.

        Returns
        -------
        list of str
            A list of the string names for each nuclide in the spatial domain
            (e.g., ['U235', 'U238', 'O16'])

        Raises
        ------
        ValueError
            When this method is called before the spatial domain has been set.

        """

        if self.domain is None:
            raise ValueError('Unable to get all nuclides without a domain')

        # If the user defined nuclides, return them
        if self._nuclides:
            return self._nuclides

        # Otherwise, return all nuclides in the spatial domain
        else:
            return self.domain.get_nuclides()

    def get_nuclide_density(self, nuclide):
        """Get the atomic number density in units of atoms/b-cm for a nuclide
        in the cross section's spatial domain.

        Parameters
        ----------
        nuclide : str
            A nuclide name string (e.g., 'U235')

        Returns
        -------
        float
            The atomic number density (atom/b-cm) for the nuclide of interest

        """

        cv.check_type('nuclide', nuclide, string_types)

        # Get list of all nuclides in the spatial domain
        nuclides = self.domain.get_nuclide_densities()

        return nuclides[nuclide][1] if nuclide in nuclides else 0.0

    def get_nuclide_densities(self, nuclides='all'):
        """Get an array of atomic number densities in units of atom/b-cm for all
        nuclides in the cross section's spatial domain.

        Parameters
        ----------
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the atom densities for all nuclides
            in the spatial domain. The special string 'sum' will return the atom
            density summed across all nuclides in the spatial domain. Defaults
            to 'all'.

        Returns
        -------
        numpy.ndarray of float
            An array of the atomic number densities (atom/b-cm) for each of the
            nuclides in the spatial domain

        Raises
        ------
        ValueError
            When this method is called before the spatial domain has been set.

        """

        if self.domain is None:
            raise ValueError('Unable to get nuclide densities without a domain')

        # Sum the atomic number densities for all nuclides
        if nuclides == 'sum':
            nuclides = self.get_nuclides()
            densities = np.zeros(1, dtype=np.float)
            for nuclide in nuclides:
                densities[0] += self.get_nuclide_density(nuclide)

        # Tabulate the atomic number densities for all nuclides
        elif nuclides == 'all':
            nuclides = self.get_nuclides()
            densities = np.zeros(self.num_nuclides, dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] += self.get_nuclide_density(nuclide)

        # Tabulate the atomic number densities for each specified nuclide
        else:
            densities = np.zeros(len(nuclides), dtype=np.float)
            for i, nuclide in enumerate(nuclides):
                densities[i] = self.get_nuclide_density(nuclide)

        return densities

    def _compute_xs(self):
        """Performs generic cleanup after a subclass' uses tally arithmetic to
        compute a multi-group cross section as a derived tally.

        This method replaces CrossNuclides generated by tally arithmetic with
        the original Nuclide objects in the xs_tally instance attribute. The
        simple Nuclides allow for cleaner output through Pandas DataFrames as
        well as simpler data access through the get_xs(...) class method.

        In addition, this routine resets NaNs in the multi group cross section
        array to 0.0. This may be needed occur if no events were scored in
        certain tally bins, which will lead to a divide-by-zero situation.

        """

        # If computing xs for each nuclide, replace CrossNuclides with originals
        if self.by_nuclide:
            self.xs_tally._nuclides = []
            nuclides = self.get_nuclides()
            for nuclide in nuclides:
                self.xs_tally.nuclides.append(openmc.Nuclide(nuclide))

        # Remove NaNs which may have resulted from divide-by-zero operations
        self.xs_tally._mean = np.nan_to_num(self.xs_tally.mean)
        self.xs_tally._std_dev = np.nan_to_num(self.xs_tally.std_dev)
        self.xs_tally.sparse = self.sparse

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

        cv.check_type('statepoint', statepoint, openmc.statepoint.StatePoint)

        if statepoint.summary is None:
            msg = 'Unable to load data from a statepoint which has not been ' \
                  'linked with a summary file'
            raise ValueError(msg)

        # Override the domain object that loaded from an OpenMC summary file
        # NOTE: This is necessary for micro cross-sections which require
        # the isotopic number densities as computed by OpenMC
        if self.domain_type == 'cell' or self.domain_type == 'distribcell':
            self.domain = statepoint.summary.get_cell_by_id(self.domain.id)
        elif self.domain_type == 'universe':
            self.domain = statepoint.summary.get_universe_by_id(self.domain.id)
        elif self.domain_type == 'material':
            self.domain = statepoint.summary.get_material_by_id(self.domain.id)
        elif self.domain_type == 'mesh':
            self.domain = statepoint.meshes[self.domain.id]
        else:
            msg = 'Unable to load data from a statepoint for domain type {0} ' \
                  'which is not yet supported'.format(self.domain_type)
            raise ValueError(msg)

        # Use tally "slicing" to ensure that tallies correspond to our domain
        # NOTE: This is important if tally merging was used
        if self.domain_type == 'mesh':
            filters = [_DOMAIN_TO_FILTER[self.domain_type]]
            xyz = [range(1, x+1) for x in self.domain.dimension]
            filter_bins = [tuple(itertools.product(*xyz))]
        elif self.domain_type != 'distribcell':
            filters = [_DOMAIN_TO_FILTER[self.domain_type]]
            filter_bins = [(self.domain.id,)]
        # Distribcell filters only accept single cell - neglect it when slicing
        else:
            filters = []
            filter_bins = []

        # Clear any tallies previously loaded from a statepoint
        if self.loaded_sp:
            self._tallies = None
            self._xs_tally = None
            self._rxn_rate_tally = None
            self._loaded_sp = False

        # Find, slice and store Tallies from StatePoint
        # The tally slicing is needed if tally merging was used
        for tally_type, tally in self.tallies.items():
            sp_tally = statepoint.get_tally(
                tally.scores, tally.filters, tally.nuclides,
                estimator=tally.estimator, exact_filters=True)
            sp_tally = sp_tally.get_slice(
                tally.scores, filters, filter_bins, tally.nuclides)
            sp_tally.sparse = self.sparse
            self.tallies[tally_type] = sp_tally

        self._loaded_sp = True

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               value='mean', squeeze=True, **kwargs):
        r"""Returns an array of multi-group cross sections.

        This method constructs a 3D NumPy array for the requested
        multi-group cross section data for one or more subdomains
        (1st dimension), energy groups (2nd dimension), and nuclides
        (3rd dimension).

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        numpy.ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, string_types):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append(openmc.EnergyFilter)
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # If user requested the sum for all nuclides, use tally summation
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(filters=filters,
                                     filter_bins=filter_bins, value=value)
        else:
            xs = self.xs_tally.get_values(filters=filters, filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        if groups == 'all':
            num_groups = self.num_groups
        else:
            num_groups = len(groups)

        # Reshape tally data array with separate axes for domain and energy
        num_subdomains = int(xs.shape[0] / num_groups)
        new_shape = (num_subdomains, num_groups) + xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[:, ::-1, :]

        if squeeze:
            xs = np.squeeze(xs)
            xs = np.atleast_1d(xs)

        return xs

    def get_condensed_xs(self, coarse_groups):
        """Construct an energy-condensed version of this cross section.

        Parameters
        ----------
        coarse_groups : openmc.mgxs.EnergyGroups
            The coarse energy group structure of interest

        Returns
        -------
        MGXS
            A new MGXS condensed to the group structure of interest

        """

        cv.check_type('coarse_groups', coarse_groups, EnergyGroups)
        cv.check_less_than('coarse groups', coarse_groups.num_groups,
                           self.num_groups, equality=True)
        cv.check_value('upper coarse energy', coarse_groups.group_edges[-1],
                       [self.energy_groups.group_edges[-1]])
        cv.check_value('lower coarse energy', coarse_groups.group_edges[0],
                       [self.energy_groups.group_edges[0]])

        # Clone this MGXS to initialize the condensed version
        condensed_xs = copy.deepcopy(self)
        condensed_xs._rxn_rate_tally = None
        condensed_xs._xs_tally = None
        condensed_xs._sparse = False
        condensed_xs._energy_groups = coarse_groups

        # Build energy indices to sum across
        energy_indices = []
        for group in range(coarse_groups.num_groups, 0, -1):
            low, high = coarse_groups.get_group_bounds(group)
            low_index = np.where(self.energy_groups.group_edges == low)[0][0]
            energy_indices.append(low_index)

        fine_edges = self.energy_groups.group_edges

        # Condense each of the tallies to the coarse group structure
        for tally in condensed_xs.tallies.values():

            # Make condensed tally derived and null out sum, sum_sq
            tally._derived = True
            tally._sum = None
            tally._sum_sq = None

            # Get tally data arrays reshaped with one dimension per filter
            mean = tally.get_reshaped_data(value='mean')
            std_dev = tally.get_reshaped_data(value='std_dev')

            # Sum across all applicable fine energy group filters
            for i, tally_filter in enumerate(tally.filters):
                if not isinstance(tally_filter, (openmc.EnergyFilter,
                                                 openmc.EnergyoutFilter)):
                    continue
                elif len(tally_filter.bins) != len(fine_edges):
                    continue
                elif not np.allclose(tally_filter.bins, fine_edges):
                    continue
                else:
                    tally_filter.bins = coarse_groups.group_edges
                    mean = np.add.reduceat(mean, energy_indices, axis=i)
                    std_dev = np.add.reduceat(std_dev**2, energy_indices, axis=i)
                    std_dev = np.sqrt(std_dev)

            # Reshape condensed data arrays with one dimension for all filters
            mean = np.reshape(mean, tally.shape)
            std_dev = np.reshape(std_dev, tally.shape)

            # Override tally's data with the new condensed data
            tally._mean = mean
            tally._std_dev = std_dev

        # Compute the energy condensed multi-group cross section
        condensed_xs.sparse = self.sparse
        return condensed_xs

    def get_subdomain_avg_xs(self, subdomains='all'):
        """Construct a subdomain-averaged version of this cross section.

        This method is useful for averaging cross sections across distribcell
        instances. The method performs spatial homogenization to compute the
        scalar flux-weighted average cross section across the subdomains.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs to average across. Defaults to 'all'.

        Returns
        -------
        openmc.mgxs.MGXS
            A new MGXS averaged across the subdomains of interest

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        # Construct a collection of the subdomain filter bins to average across
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains)
        else:
            subdomains = None

        # Clone this MGXS to initialize the subdomain-averaged version
        avg_xs = copy.deepcopy(self)

        if self.derived:
            avg_xs._rxn_rate_tally = avg_xs.rxn_rate_tally.average(
                filter_type=_DOMAIN_TO_FILTER[self.domain_type],
                filter_bins=subdomains)
        else:
            avg_xs._rxn_rate_tally = None
            avg_xs._xs_tally = None

            # Average each of the tallies across subdomains
            for tally_type, tally in avg_xs.tallies.items():
                filt_type = _DOMAIN_TO_FILTER[self.domain_type]
                tally_avg = tally.average(filter_type=filt_type,
                                          filter_bins=subdomains)
                avg_xs.tallies[tally_type] = tally_avg

        avg_xs._domain_type = 'avg({0})'.format(self.domain_type)
        avg_xs.sparse = self.sparse
        return avg_xs

    def get_slice(self, nuclides=[], groups=[]):
        """Build a sliced MGXS for the specified nuclides and energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        groups : list of int
            A list of energy group indices starting at 1 for the high energies
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MGXS
            A new MGXS object which encapsulates the subset of data requested
            for the nuclide(s) and/or energy group(s) requested in the
            parameters.

        """

        cv.check_iterable_type('nuclides', nuclides, string_types)
        cv.check_iterable_type('energy_groups', groups, Integral)

        # Build lists of filters and filter bins to slice
        filters = []
        filter_bins = []

        if len(groups) != 0:
            energy_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                energy_bins.append(group_bounds)
            filter_bins.append(tuple(energy_bins))
            filters.append(openmc.EnergyFilter)

        # Clone this MGXS to initialize the sliced version
        slice_xs = copy.deepcopy(self)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice each of the tallies across nuclides and energy groups
        for tally_type, tally in slice_xs.tallies.items():
            slice_nuclides = [nuc for nuc in nuclides if nuc in tally.nuclides]
            if len(groups) != 0 and tally.contains_filter(openmc.EnergyFilter):
                tally_slice = tally.get_slice(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=slice_nuclides)
            else:
                tally_slice = tally.get_slice(nuclides=slice_nuclides)
            slice_xs.tallies[tally_type] = tally_slice

        # Assign sliced energy group structure to sliced MGXS
        if groups:
            new_group_edges = []
            for group in groups:
                group_edges = self.energy_groups.get_group_bounds(group)
                new_group_edges.extend(group_edges)
            new_group_edges = np.unique(new_group_edges)
            slice_xs.energy_groups.group_edges = sorted(new_group_edges)

        # Assign sliced nuclides to sliced MGXS
        if nuclides:
            slice_xs.nuclides = nuclides

        slice_xs.sparse = self.sparse
        return slice_xs

    def can_merge(self, other):
        """Determine if another MGXS can be merged with this one

        If results have been loaded from a statepoint, then MGXS are only
        mergeable along one and only one of enegy groups or nuclides.

        Parameters
        ----------
        other : openmc.mgxs.MGXS
            MGXS to check for merging

        """

        if not isinstance(other, type(self)):
            return False

        # Compare reaction type, energy groups, nuclides, domain type
        if self.rxn_type != other.rxn_type:
            return False
        elif not self.energy_groups.can_merge(other.energy_groups):
            return False
        elif self.by_nuclide != other.by_nuclide:
            return False
        elif self.domain_type != other.domain_type:
            return False
        elif 'distribcell' not in self.domain_type and self.domain != other.domain:
            return False
        elif not self.xs_tally.can_merge(other.xs_tally):
            return False
        elif not self.rxn_rate_tally.can_merge(other.rxn_rate_tally):
            return False

        # If all conditionals pass then MGXS are mergeable
        return True

    def merge(self, other):
        """Merge another MGXS with this one

        MGXS are only mergeable if their energy groups and nuclides are either
        identical or mutually exclusive. If results have been loaded from a
        statepoint, then MGXS are only mergeable along one and only one of
        energy groups or nuclides.

        Parameters
        ----------
        other : openmc.mgxs.MGXS
            MGXS to merge with this one

        Returns
        -------
        merged_mgxs : openmc.mgxs.MGXS
            Merged MGXS

        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge MGXS')

        # Create deep copy of tally to return as merged tally
        merged_mgxs = copy.deepcopy(self)
        merged_mgxs._derived = True

        # Merge energy groups
        if self.energy_groups != other.energy_groups:
            merged_groups = self.energy_groups.merge(other.energy_groups)
            merged_mgxs.energy_groups = merged_groups

        # Merge nuclides
        if self.nuclides != other.nuclides:

            # The nuclides must be mutually exclusive
            for nuclide in self.nuclides:
                if nuclide in other.nuclides:
                    msg = 'Unable to merge MGXS with shared nuclides'
                    raise ValueError(msg)

            # Concatenate lists of nuclides for the merged MGXS
            merged_mgxs.nuclides = self.nuclides + other.nuclides

        # Null base tallies but merge reaction rate and cross section tallies
        merged_mgxs._tallies = OrderedDict()
        merged_mgxs._rxn_rate_tally = self.rxn_rate_tally.merge(other.rxn_rate_tally)
        merged_mgxs._xs_tally = self.xs_tally.merge(other.xs_tally)

        return merged_mgxs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Print a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x+1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
            elif nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, string_types)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Generate the header for an individual XS
        xs_header = '\tCross Sections [{0}]:'.format(self.get_units(xs_type))

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell' or self.domain_type == 'mesh':
                string += '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if nuclide != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                string += '{0: <16}\n'.format(xs_header)
                template = '{0: <12}Group {1} [{2: <10} - {3: <10}eV]:\t'

                # Loop over energy groups ranges
                for group in range(1, self.num_groups+1):
                    bounds = self.energy_groups.get_group_bounds(group)
                    string += template.format('', group, bounds[0], bounds[1])
                    average = self.get_xs([group], [subdomain], [nuclide],
                                          xs_type=xs_type, value='mean')
                    rel_err = self.get_xs([group], [subdomain], [nuclide],
                                          xs_type=xs_type, value='rel_err')
                    average = average.flatten()[0]
                    rel_err = rel_err.flatten()[0] * 100.
                    string += '{:.2e} +/- {:1.2e}%'.format(average, rel_err)
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)

    def build_hdf5_store(self, filename='mgxs.h5', directory='mgxs',
                         subdomains='all', nuclides='all',
                         xs_type='macro', row_column='inout', append=True):
        """Export the multi-group cross section data to an HDF5 binary file.

        This method constructs an HDF5 file which stores the multi-group
        cross section data. The data is stored in a hierarchy of HDF5 groups
        from the domain type, domain id, subdomain id (for distribcell domains),
        nuclides and cross section type. Two datasets for the mean and standard
        deviation are stored for each subdomain entry in the HDF5 file.

        NOTE: This requires the h5py Python package.

        Parameters
        ----------
        filename : str
            Filename for the HDF5 file. Defaults to 'mgxs.h5'.
        directory : str
            Directory for the HDF5 file. Defaults to 'mgxs'.
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
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
        append : bool
            If true, appends to an existing HDF5 file with the same filename
            directory (if one exists). Defaults to True.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.
        ImportError
            When h5py is not installed.

        """

        import h5py

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)
        filename = filename.replace(' ', '-')

        if append and os.path.isfile(filename):
            xs_results = h5py.File(filename, 'a')
        else:
            xs_results = h5py.File(filename, 'w')

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'avg(distribcell)':
            domain_filter = self.xs_tally.find_filter('avg(distribcell)')
            subdomains = domain_filter.bins
        elif self.domain_type == 'mesh':
            xyz = [range(1, x+1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
                densities = np.zeros(len(nuclides), dtype=np.float)
            elif nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, string_types)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Create an HDF5 group within the file for the domain
        domain_type_group = xs_results.require_group(self.domain_type)
        domain_group = domain_type_group.require_group(str(self.domain.id))

        # Determine number of digits to pad subdomain group keys
        num_digits = len(str(self.num_subdomains))

        # Create a separate HDF5 group for each subdomain
        for subdomain in subdomains:

            # Create an HDF5 group for the subdomain
            if self.domain_type == 'distribcell':
                group_name = ''.zfill(num_digits)
                subdomain_group = domain_group.require_group(group_name)
            else:
                subdomain_group = domain_group

            # Create a separate HDF5 group for this cross section
            rxn_group = subdomain_group.require_group(self.hdf5_key)

            # Create a separate HDF5 group for each nuclide
            for j, nuclide in enumerate(nuclides):

                if nuclide != 'sum':
                    density = densities[j]
                    nuclide_group = rxn_group.require_group(nuclide)
                    nuclide_group.require_dataset('density', dtype=np.float64,
                                                  data=[density], shape=(1,))
                else:
                    nuclide_group = rxn_group

                # Extract the cross section for this subdomain and nuclide
                average = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='mean',
                                      row_column=row_column)
                std_dev = self.get_xs(subdomains=[subdomain], nuclides=[nuclide],
                                      xs_type=xs_type, value='std_dev',
                                      row_column=row_column)

                # Add MGXS results data to the HDF5 group
                nuclide_group.require_dataset('average', dtype=np.float64,
                                              shape=average.shape, data=average)
                nuclide_group.require_dataset('std. dev.', dtype=np.float64,
                                              shape=std_dev.shape, data=std_dev)

        # Close the results HDF5 file
        xs_results.close()

    def export_xs_data(self, filename='mgxs', directory='mgxs',
                       format='csv', groups='all', xs_type='macro'):
        """Export the multi-group cross section data to a file.

        This method leverages the functionality in the Pandas library to export
        the multi-group cross section data in a variety of output file formats
        for storage and/or post-processing.

        Parameters
        ----------
        filename : str
            Filename for the exported file. Defaults to 'mgxs'.
        directory : str
            Directory for the exported file. Defaults to 'mgxs'.
        format : {'csv', 'excel', 'pickle', 'latex'}
            The format for the exported data file. Defaults to 'csv'.
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Store the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        cv.check_type('filename', filename, string_types)
        cv.check_type('directory', directory, string_types)
        cv.check_value('format', format, ['csv', 'excel', 'pickle', 'latex'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Make directory if it does not exist
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename = os.path.join(directory, filename)
        filename = filename.replace(' ', '-')

        # Get a Pandas DataFrame for the data
        df = self.get_pandas_dataframe(groups=groups, xs_type=xs_type)

        # Export the data using Pandas IO API
        if format == 'csv':
            df.to_csv(filename + '.csv', index=False)
        elif format == 'excel':
            if self.domain_type == 'mesh':
                df.to_excel(filename + '.xls')
            else:
                df.to_excel(filename + '.xls', index=False)
        elif format == 'pickle':
            df.to_pickle(filename + '.pkl')
        elif format == 'latex':
            if self.domain_type == 'distribcell':
                msg = 'Unable to export distribcell multi-group cross section' \
                      'data to a LaTeX table'
                raise NotImplementedError(msg)

            df.to_latex(filename + '.tex', bold_rows=True,
                        longtable=True, index=False)

            # Surround LaTeX table with code needed to run pdflatex
            with open(filename + '.tex', 'r') as original:
                data = original.read()
            with open(filename + '.tex', 'w') as modified:
                modified.write(
                    '\\documentclass[preview, 12pt, border=1mm]{standalone}\n')
                modified.write('\\usepackage{caption}\n')
                modified.write('\\usepackage{longtable}\n')
                modified.write('\\usepackage{booktabs}\n')
                modified.write('\\begin{document}\n\n')
                modified.write(data)
                modified.write('\n\\end{document}')

    def get_pandas_dataframe(self, groups='all', nuclides='all',
                             xs_type='macro', distribcell_paths=True):
        """Build a Pandas DataFrame for the MGXS data.

        This method leverages :meth:`openmc.Tally.get_pandas_dataframe`, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will include the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            include the cross sections summed over all nuclides. Defaults
            to 'all'.
        xs_type: {'macro', 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        distribcell_paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into
            a Multi-index column with a geometric "path" to each distribcell
            instance.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        if not isinstance(groups, string_types):
            cv.check_iterable_type('groups', groups, Integral)
        if nuclides != 'all' and nuclides != 'sum':
            cv.check_iterable_type('nuclides', nuclides, string_types)
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Get a Pandas DataFrame from the derived xs tally
        if self.by_nuclide and nuclides == 'sum':

            # Use tally summation to sum across all nuclides
            query_nuclides = [nuclides]
            xs_tally = self.xs_tally.summation(nuclides=self.get_nuclides())
            df = xs_tally.get_pandas_dataframe(
                distribcell_paths=distribcell_paths)

            # Remove nuclide column since it is homogeneous and redundant
            if self.domain_type == 'mesh':
                df.drop('sum(nuclide)', axis=1, level=0, inplace=True)
            else:
                df.drop('sum(nuclide)', axis=1, inplace=True)

        # If the user requested a specific set of nuclides
        elif self.by_nuclide and nuclides != 'all':
            query_nuclides = nuclides
            xs_tally = self.xs_tally.get_slice(nuclides=nuclides)
            df = xs_tally.get_pandas_dataframe(
                distribcell_paths=distribcell_paths)

        # If the user requested all nuclides, keep nuclide column in dataframe
        else:
            query_nuclides = self.nuclides
            df = self.xs_tally.get_pandas_dataframe(
                distribcell_paths=distribcell_paths)

        # Remove the score column since it is homogeneous and redundant
        if self.domain_type == 'mesh':
            df = df.drop('score', axis=1, level=0)
        else:
            df = df.drop('score', axis=1)

        # Override energy groups bounds with indices
        all_groups = np.arange(self.num_groups, 0, -1, dtype=np.int)
        all_groups = np.repeat(all_groups, len(query_nuclides))
        if 'energy low [eV]' in df and 'energyout low [eV]' in df:
            df.rename(columns={'energy low [eV]': 'group in'},
                      inplace=True)
            in_groups = np.tile(all_groups, int(self.num_subdomains))
            in_groups = np.repeat(in_groups, int(df.shape[0] / in_groups.size))
            df['group in'] = in_groups
            del df['energy high [eV]']

            df.rename(columns={'energyout low [eV]': 'group out'},
                      inplace=True)
            out_groups = np.repeat(all_groups, self.xs_tally.num_scores)
            out_groups = np.tile(out_groups, int(df.shape[0] / out_groups.size))
            df['group out'] = out_groups
            del df['energyout high [eV]']
            columns = ['group in', 'group out']

        elif 'energyout low [eV]' in df:
            df.rename(columns={'energyout low [eV]': 'group out'},
                      inplace=True)
            in_groups = np.tile(all_groups, int(df.shape[0] / all_groups.size))
            df['group out'] = in_groups
            del df['energyout high [eV]']
            columns = ['group out']

        elif 'energy low [eV]' in df:
            df.rename(columns={'energy low [eV]': 'group in'}, inplace=True)
            in_groups = np.tile(all_groups, int(df.shape[0] / all_groups.size))
            df['group in'] = in_groups
            del df['energy high [eV]']
            columns = ['group in']

        # Select out those groups the user requested
        if not isinstance(groups, string_types):
            if 'group in' in df:
                df = df[df['group in'].isin(groups)]
            if 'group out' in df:
                df = df[df['group out'].isin(groups)]

        # If user requested micro cross sections, divide out the atom densities
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            densities = np.repeat(densities, len(self.rxn_rate_tally.scores))
            tile_factor = df.shape[0] / len(densities)
            df['mean'] /= np.tile(densities, tile_factor)
            df['std. dev.'] /= np.tile(densities, tile_factor)

            # Replace NaNs by zeros (happens if nuclide density is zero)
            df['mean'].replace(np.nan, 0.0, inplace=True)
            df['std. dev.'].replace(np.nan, 0.0, inplace=True)

        # Sort the dataframe by domain type id (e.g., distribcell id) and
        # energy groups such that data is from fast to thermal
        if self.domain_type == 'mesh':
            mesh_str = 'mesh {0}'.format(self.domain.id)
            df.sort_values(by=[(mesh_str, 'x'), (mesh_str, 'y'), \
                               (mesh_str, 'z')] + columns, inplace=True)
        else:
            df.sort_values(by=[self.domain_type] + columns, inplace=True)

        return df

    def get_units(self, xs_type='macro'):
        """This method returns the units of a MGXS based on a desired xs_type.

        Parameters
        ----------
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section units.
            Defaults to 'macro'.

        Returns
        -------
        str
            A string representing the units of the MGXS.

        """

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        return 'cm^-1' if xs_type == 'macro' else 'barns'


@add_metaclass(ABCMeta)
class MatrixMGXS(MGXS):
    """An abstract multi-group cross section for some energy group structure
    within some spatial domain. This class is specifically intended for
    cross sections which depend on both the incoming and outgoing energy groups
    and are therefore represented by matrices. Examples of this include the
    scattering and nu-fission matrices.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations.

    NOTE: Users should instantiate the subclasses of this abstract class.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """
    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy = openmc.EnergyFilter(group_edges)
        energyout = openmc.EnergyoutFilter(group_edges)

        return [[energy], [energy, energyout]]

    def get_xs(self, in_groups='all', out_groups='all',
               subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               row_column='inout', value='mean', squeeze=True, **kwargs):
        """Returns an array of multi-group cross sections.

        This method constructs a 4D NumPy array for the requested
        multi-group cross section data for one or more subdomains
        (1st dimension), energy groups in (2nd dimension), energy groups out
        (3rd dimension), and nuclides (4th dimension).

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest. Defaults to 'all'.
        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            return the cross section summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        row_column: {'inout', 'outin'}
            Return the cross section indexed first by incoming group and
            second by outgoing group ('inout'), or vice versa ('outin').
            Defaults to 'inout'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral,
                                   max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(in_groups, string_types):
            cv.check_iterable_type('groups', in_groups, Integral)
            for group in in_groups:
                filters.append(openmc.EnergyFilter)
                filter_bins.append((
                    self.energy_groups.get_group_bounds(group),))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(out_groups, string_types):
            cv.check_iterable_type('groups', out_groups, Integral)
            for group in out_groups:
                filters.append(openmc.EnergyoutFilter)
                filter_bins.append((
                    self.energy_groups.get_group_bounds(group),))

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # Use tally summation if user requested the sum for all nuclides
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(filters=filters, filter_bins=filter_bins,
                                     value=value)
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        if in_groups == 'all':
            num_in_groups = self.num_groups
        else:
            num_in_groups = len(in_groups)

        if out_groups == 'all':
            num_out_groups = self.num_groups
        else:
            num_out_groups = len(out_groups)

        # Reshape tally data array with separate axes for domain and energy
        num_subdomains = int(xs.shape[0] / (num_in_groups * num_out_groups))
        new_shape = (num_subdomains, num_in_groups, num_out_groups)
        new_shape += xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Transpose the matrix if requested by user
        if row_column == 'outin':
            xs = np.swapaxes(xs, 1, 2)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[:, ::-1, ::-1, :]

        if squeeze:
            xs = np.squeeze(xs)
            xs = np.atleast_2d(xs)

        return xs

    def get_slice(self, nuclides=[], in_groups=[], out_groups=[]):
        """Build a sliced MatrixMGXS object for the specified nuclides and
        energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        in_groups : list of int
            A list of incoming energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        out_groups : list of int
            A list of outgoing energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MatrixMGXS
            A new MatrixMGXS object which encapsulates the subset of data
            requested for the nuclide(s) and/or energy group(s) requested in
            the parameters.

        """

        # Call super class method and null out derived tallies
        slice_xs = super(MatrixMGXS, self).get_slice(nuclides, in_groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice outgoing energy groups if needed
        if len(out_groups) != 0:
            filter_bins = []
            for group in out_groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]

            # Slice each of the tallies across energyout groups
            for tally_type, tally in slice_xs.tallies.items():
                if tally.contains_filter(openmc.EnergyoutFilter):
                    tally_slice = tally.get_slice(
                        filters=[openmc.EnergyoutFilter],
                        filter_bins=filter_bins)
                    slice_xs.tallies[tally_type] = tally_slice

        slice_xs.sparse = self.sparse
        return slice_xs

    def print_xs(self, subdomains='all', nuclides='all', xs_type='macro'):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            report the cross sections summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x+1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
            if nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, string_types)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', self.rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Generate the header for an individual XS
        xs_header = '\tCross Sections [{0}]:'.format(self.get_units(xs_type))

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        string += '{0: <16}\n'.format('\tEnergy Groups:')
        template = '{0: <12}Group {1} [{2: <10} - {3: <10}eV]\n'

        # Loop over energy groups ranges
        for group in range(1, self.num_groups + 1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += template.format('', group, bounds[0], bounds[1])

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += \
                    '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if xs_type != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                string += '{0: <16}\n'.format(xs_header)
                template = '{0: <12}Group {1} -> Group {2}:\t\t'

                # Loop over incoming/outgoing energy groups ranges
                for in_group in range(1, self.num_groups + 1):
                    for out_group in range(1, self.num_groups + 1):
                        string += template.format('', in_group, out_group)
                        average = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide],
                                        xs_type=xs_type, value='mean')
                        rel_err = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide],
                                        xs_type=xs_type, value='rel_err')
                        average = average.flatten()[0]
                        rel_err = rel_err.flatten()[0] * 100.
                        string += '{:1.2e} +/- {:1.2e}%'.format(average,
                                                                rel_err)
                        string += '\n'
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)


class TotalXS(MGXS):
    r"""A total multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group total cross sections for multi-group neutronics calculations. At
    a minimum, one needs to set the :attr:`TotalXS.energy_groups` and
    :attr:`TotalXS.domain` properties. Tallies for the flux and appropriate
    reaction rates over the specified domain are generated automatically via the
    :attr:`TotalXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`TotalXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    total cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \sigma_t (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`TotalXS.tally_keys` property and values
        are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TotalXS, self).__init__(domain, domain_type,
                                      groups, by_nuclide, name)
        self._rxn_type = 'total'


class TransportXS(MGXS):
    r"""A transport-corrected total multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`TransportXS.energy_groups` and
    :attr:`TransportXS.domain` properties. Tallies for the flux and appropriate
    reaction rates over the specified domain are generated automatically via the
    :attr:`TransportXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`TransportXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    transport-corrected total cross section is calculated as:

    .. math::

       \langle \sigma_t \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \sigma_t (r, E) \psi
       (r, E, \Omega) \\
       \langle \sigma_{s1} \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_{-1}^1 d\mu \; \mu \sigma_s
       (r, E' \rightarrow E, \Omega' \cdot \Omega)
       \phi (r, E', \Omega) \\
       \langle \phi \rangle &= \int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega) \\
       \sigma_{tr} &= \frac{\langle \sigma_t \phi \rangle - \langle \sigma_{s1}
       \phi \rangle}{\langle \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`TransportXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(TransportXS, self).__init__(domain, domain_type,
                                          groups, by_nuclide, name)
        self._rxn_type = 'transport'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']

    @property
    def scores(self):
        return ['flux', 'total', 'scatter-1']

    @property
    def filters(self):
        group_edges = self.energy_groups.group_edges
        energy_filter = openmc.EnergyFilter(group_edges)
        energyout_filter = openmc.EnergyoutFilter(group_edges)
        return [[energy_filter], [energy_filter], [energyout_filter]]

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            # Switch EnergyoutFilter to EnergyFilter.
            old_filt = self.tallies['scatter-1'].filters[-1]
            new_filt = openmc.EnergyFilter(old_filt.bins)
            new_filt.stride = old_filt.stride
            self.tallies['scatter-1'].filters[-1] = new_filt

            self._rxn_rate_tally = \
                self.tallies['total'] - self.tallies['scatter-1']
            self._rxn_rate_tally.sparse = self.sparse

        return self._rxn_rate_tally


class NuTransportXS(TransportXS):
    r"""A transport-corrected total multi-group cross section which
    accounts for neutron multiplicity in scattering reactions.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`NuTransportXS.energy_groups` and
    :attr:`NuTransportXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`NuTransportXS.tallies` property, which can then
    be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`NuTransportXS.xs_tally` property.

    The calculation of the transport-corrected cross section is the same as that
    for :class:`TransportXS` except that the scattering multiplicity is
    accounted for.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`NuTransportXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuTransportXS, self).__init__(domain, domain_type,
                                            groups, by_nuclide, name)
        self._rxn_type = 'nu-transport'

    @property
    def scores(self):
        return ['flux', 'total', 'nu-scatter-1']

    @property
    def tally_keys(self):
        return ['flux', 'total', 'scatter-1']


class AbsorptionXS(MGXS):
    r"""An absorption multi-group cross section.

    Absorption is defined as all reactions that do not produce secondary
    neutrons (disappearance) plus fission reactions.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group absorption cross sections for multi-group neutronics
    calculations. At a minimum, one needs to set the
    :attr:`AbsorptionXS.energy_groups` and :attr:`AbsorptionXS.domain`
    properties. Tallies for the flux and appropriate reaction rates over the
    specified domain are generated automatically via the
    :attr:`AbsorptionXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`AbsorptionXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    absorption cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \sigma_a (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`AbsorptionXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(AbsorptionXS, self).__init__(domain, domain_type,
                                           groups, by_nuclide, name)
        self._rxn_type = 'absorption'


class CaptureXS(MGXS):
    r"""A capture multi-group cross section.

    The neutron capture reaction rate is defined as the difference between
    OpenMC's 'absorption' and 'fission' reaction rate score types. This includes
    not only radiative capture, but all forms of neutron disappearance aside
    from fission (i.e., MT > 100).

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group capture cross sections for multi-group neutronics
    calculations. At a minimum, one needs to set the
    :attr:`CaptureXS.energy_groups` and :attr:`CaptureXS.domain`
    properties. Tallies for the flux and appropriate reaction rates over the
    specified domain are generated automatically via the
    :attr:`CaptureXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`CaptureXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    capture cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \left [ \sigma_a (r, E) \psi (r, E, \Omega) - \sigma_f (r, E) \psi (r, E,
       \Omega) \right ]}{\int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`CaptureXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(CaptureXS, self).__init__(domain, domain_type,
                                        groups, by_nuclide, name)
        self._rxn_type = 'capture'

    @property
    def scores(self):
        return ['flux', 'absorption', 'fission']

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = \
                self.tallies['absorption'] - self.tallies['fission']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally


class FissionXS(MGXS):
    r"""A fission multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group fission cross sections for multi-group neutronics
    calculations. At a minimum, one needs to set the
    :attr:`FissionXS.energy_groups` and :attr:`FissionXS.domain`
    properties. Tallies for the flux and appropriate reaction rates over the
    specified domain are generated automatically via the
    :attr:`FissionXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`FissionXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    fission cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \sigma_f (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`FissionXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(FissionXS, self).__init__(domain, domain_type,
                                        groups, by_nuclide, name)
        self._rxn_type = 'fission'


class NuFissionXS(MGXS):
    r"""A fission neutron production multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group fission neutron production cross sections for multi-group
    neutronics calculations. At a minimum, one needs to set the
    :attr:`NuFissionXS.energy_groups` and :attr:`NuFissionXS.domain`
    properties. Tallies for the flux and appropriate reaction rates over the
    specified domain are generated automatically via the
    :attr:`NuFissionXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`NuFissionXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    fission neutron production cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \nu\sigma_f (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.


    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`NuFissionXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuFissionXS, self).__init__(domain, domain_type,
                                          groups, by_nuclide, name)
        self._rxn_type = 'nu-fission'

class KappaFissionXS(MGXS):
    r"""A recoverable fission energy production rate multi-group cross section.

    The recoverable energy per fission, :math:`\kappa`, is defined as the
    fission product kinetic energy, prompt and delayed neutron kinetic energies,
    prompt and delayed :math:`\gamma`-ray total energies, and the total energy
    released by the delayed :math:`\beta` particles. The neutrino energy does
    not contribute to this response. The prompt and delayed :math:`\gamma`-rays
    are assumed to deposit their energy locally.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`KappaFissionXS.energy_groups` and
    :attr:`KappaFissionXS.domain` properties. Tallies for the flux and appropriate
    reaction rates over the specified domain are generated automatically via the
    :attr:`KappaFissionXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`KappaFissionXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    recoverable fission energy production rate cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \kappa\sigma_f (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`KappaFissionXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(KappaFissionXS, self).__init__(domain, domain_type,
                                             groups, by_nuclide, name)
        self._rxn_type = 'kappa-fission'


class ScatterXS(MGXS):
    r"""A scattering multi-group cross section.

    The scattering cross section is defined as the difference between the total
    and absorption cross sections.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`ScatterXS.energy_groups` and
    :attr:`ScatterXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`ScatterXS.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`ScatterXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    scattering cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \left [ \sigma_t (r, E) \psi (r, E, \Omega) - \sigma_a (r, E) \psi (r, E,
       \Omega) \right ]}{\int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`ScatterXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterXS, self).__init__(domain, domain_type,
                                        groups, by_nuclide, name)
        self._rxn_type = 'scatter'


class NuScatterXS(MGXS):
    r"""A scattering neutron production multi-group cross section.

    The neutron production from scattering is defined as the average number of
    neutrons produced from all neutron-producing reactions except for fission.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`NuScatterXS.energy_groups` and
    :attr:`NuScatterXS.domain` properties. Tallies for the flux and appropriate
    reaction rates over the specified domain are generated automatically via the
    :attr:`NuScatterXS.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`NuScatterXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    scattering neutron production cross section is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \sum_i \upsilon_i \sigma_i (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr
       \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    where :math:`\upsilon_i` is the multiplicity of the :math:`i`-th scattering
    reaction.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`NuScatterXS.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterXS, self).__init__(domain, domain_type,
                                          groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']


class ScatterMatrixXS(MatrixMGXS):
    r"""A scattering matrix multi-group cross section for one or more Legendre
    moments.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`ScatterMatrixXS.energy_groups` and
    :attr:`ScatterMatrixXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`ScatterMatrixXS.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`ScatterMatrixXS.xs_tally` property.

    For a spatial domain :math:`V`, incoming energy group
    :math:`[E_{g'},E_{g'-1}]`, and outgoing energy group :math:`[E_g,E_{g-1}]`,
    the scattering moments are calculated as:

    .. math::

       \langle \sigma_{s,\ell,g'\rightarrow g} \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega' \int_{E_{g'}}^{E_{g'-1}} dE' \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; P_\ell (\Omega \cdot \Omega') \sigma_s (r, E'
       \rightarrow E, \Omega' \cdot \Omega) \psi(r, E', \Omega')\\
       \langle \phi \rangle &= \int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega) \\
       \sigma_{s,\ell,g'\rightarrow g} &= \frac{\langle
       \sigma_{s,\ell,g'\rightarrow g} \phi \rangle}{\langle \phi \rangle}

    If the order is zero and a :math:`P_0` transport-correction is applied
    (default), the scattering matrix elements are:

    .. math::

       \sigma_{s,g'\rightarrow g} = \frac{\langle \sigma_{s,0,g'\rightarrow g}
       \phi \rangle - \delta_{gg'} \sum_{g''} \langle \sigma_{s,1,g''\rightarrow
       g} \phi \rangle}{\langle \phi \rangle}


    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    correction : 'P0' or None
        Apply the P0 correction to scattering matrices if set to 'P0'
    legendre_order : int
        The highest Legendre moment in the scattering matrix (default is 0)
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`ScatterMatrixXS.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ScatterMatrixXS, self).__init__(domain, domain_type,
                                              groups, by_nuclide, name)
        self._rxn_type = 'scatter'
        self._correction = 'P0'
        self._legendre_order = 0
        self._hdf5_key = 'scatter matrix'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']

    def __deepcopy__(self, memo):
        clone = super(ScatterMatrixXS, self).__deepcopy__(memo)
        clone._correction = self.correction
        clone._legendre_order = self.legendre_order
        return clone

    @property
    def correction(self):
        return self._correction

    @property
    def legendre_order(self):
        return self._legendre_order

    @property
    def scores(self):
        scores = ['flux']

        if self.correction == 'P0' and self.legendre_order == 0:
            scores += ['{}-0'.format(self.rxn_type),
                       '{}-1'.format(self.rxn_type)]
        else:
            scores += ['{}-P{}'.format(self.rxn_type, self.legendre_order)]

        return scores

    @property
    def filters(self):
        group_edges = self.energy_groups.group_edges
        energy = openmc.EnergyFilter(group_edges)
        energyout = openmc.EnergyoutFilter(group_edges)

        if self.correction == 'P0' and self.legendre_order == 0:
            filters = [[energy], [energy, energyout], [energyout]]
        else:
            filters = [[energy], [energy, energyout]]

        return filters

    @property
    def rxn_rate_tally(self):

        if self._rxn_rate_tally is None:

            # If using P0 correction subtract scatter-1 from the diagonal
            if self.correction == 'P0' and self.legendre_order == 0:
                scatter_p0 = self.tallies['{}-0'.format(self.rxn_type)]
                scatter_p1 = self.tallies['{}-1'.format(self.rxn_type)]
                energy_filter = scatter_p0.find_filter(openmc.EnergyFilter)
                energy_filter = copy.deepcopy(energy_filter)
                scatter_p1 = scatter_p1.diagonalize_filter(energy_filter)
                self._rxn_rate_tally = scatter_p0 - scatter_p1

            # Extract scattering moment reaction rate Tally
            else:
                tally_key = '{}-P{}'.format(self.rxn_type, self.legendre_order)
                self._rxn_rate_tally = self.tallies[tally_key]

            self._rxn_rate_tally.sparse = self.sparse

        return self._rxn_rate_tally

    @correction.setter
    def correction(self, correction):
        cv.check_value('correction', correction, ('P0', None))

        if correction == 'P0' and self.legendre_order > 0:
            msg = 'The P0 correction will be ignored since the scattering ' \
                  'order {} is greater than zero'.format(self.legendre_order)
            warnings.warn(msg)

        self._correction = correction

    @legendre_order.setter
    def legendre_order(self, legendre_order):
        cv.check_type('legendre_order', legendre_order, Integral)
        cv.check_greater_than('legendre_order', legendre_order, 0, equality=True)
        cv.check_less_than('legendre_order', legendre_order, 10, equality=True)

        if self.correction == 'P0' and legendre_order > 0:
            msg = 'The P0 correction will be ignored since the scattering ' \
                  'order {} is greater than zero'.format(self.legendre_order)
            warnings.warn(msg, RuntimeWarning)
            self.correction = None

        self._legendre_order = legendre_order

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

        # Clear any tallies previously loaded from a statepoint
        if self.loaded_sp:
            self._tallies = None
            self._xs_tally = None
            self._rxn_rate_tally = None
            self._loaded_sp = False

        # Expand scores to match the format in the statepoint
        # e.g., "scatter-P2" -> "scatter-0", "scatter-1", "scatter-2"
        if self.correction != 'P0' or self.legendre_order != 0:
            tally_key = '{}-P{}'.format(self.rxn_type, self.legendre_order)
            self.tallies[tally_key].scores = \
                [self.rxn_type + '-{}'.format(i) for i in range(self.legendre_order+1)]

        super(ScatterMatrixXS, self).load_from_statepoint(statepoint)

    def get_slice(self, nuclides=[], in_groups=[], out_groups=[],
                  legendre_order='same'):
        """Build a sliced ScatterMatrix for the specified nuclides and
        energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        in_groups : list of int
            A list of incoming energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        out_groups : list of int
            A list of outgoing energy group indices starting at 1 for the high
            energies (e.g., [1, 2, 3]; default is [])
        legendre_order : int or 'same'
            The highest Legendre moment in the sliced MGXS. If order is 'same'
            then the sliced MGXS will have the same Legendre moments as the
            original MGXS (default). If order is an integer less than the
            original MGXS' order, then only those Legendre moments up to that
            order will be included in the sliced MGXS.

        Returns
        -------
        openmc.mgxs.MatrixMGXS
            A new MatrixMGXS which encapsulates the subset of data requested
            for the nuclide(s) and/or energy group(s) requested in the
            parameters.

        """

        # Call super class method and null out derived tallies
        slice_xs = super(ScatterMatrixXS, self).get_slice(nuclides, in_groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice the Legendre order if needed
        if legendre_order != 'same':
            cv.check_type('legendre_order', legendre_order, Integral)
            cv.check_less_than('legendre_order', legendre_order,
                               self.legendre_order, equality=True)
            slice_xs.legendre_order = legendre_order

            # Slice the scattering tally
            tally_key = '{}-P{}'.format(self.rxn_type, self.legendre_order)
            expand_scores = \
                [self.rxn_type + '-{}'.format(i) for i in range(self.legendre_order+1)]
            slice_xs.tallies[tally_key] = \
                slice_xs.tallies[tally_key].get_slice(scores=expand_scores)

        # Slice outgoing energy groups if needed
        if len(out_groups) != 0:
            filter_bins = []
            for group in out_groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]

            # Slice each of the tallies across energyout groups
            for tally_type, tally in slice_xs.tallies.items():
                if tally.contains_filter(openmc.EnergyoutFilter):
                    tally_slice = tally.get_slice(
                        filters=[openmc.EnergyoutFilter], filter_bins=filter_bins)
                    slice_xs.tallies[tally_type] = tally_slice

        slice_xs.sparse = self.sparse
        return slice_xs

    def get_xs(self, in_groups='all', out_groups='all',
               subdomains='all', nuclides='all', moment='all',
               xs_type='macro', order_groups='increasing',
               row_column='inout', value='mean', squeeze=True):
        r"""Returns an array of multi-group cross sections.

        This method constructs a 5D NumPy array for the requested
        multi-group cross section data for one or more subdomains
        (1st dimension), energy groups in (2nd dimension), energy groups out
        (3rd dimension), nuclides (4th dimension), and moments (5th dimension).

        NOTE: The scattering moments are not multiplied by the :math:`(2l+1)/2`
        prefactor in the expansion of the scattering source into Legendre
        moments in the neutron transport equation.

        Parameters
        ----------
        in_groups : Iterable of Integral or 'all'
            Incoming energy groups of interest. Defaults to 'all'.
        out_groups : Iterable of Integral or 'all'
            Outgoing energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        moment : int or 'all'
            The scattering matrix moment to return. All moments will be
             returned if the moment is 'all' (default); otherwise, a specific
             moment will be returned.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        row_column: {'inout', 'outin'}
            Return the cross section indexed first by incoming group and
            second by outgoing group ('inout'), or vice versa ('outin').
            Defaults to 'inout'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group and subdomain is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral, max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(in_groups, string_types):
            cv.check_iterable_type('groups', in_groups, Integral)
            for group in in_groups:
                filters.append(openmc.EnergyFilter)
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(out_groups, string_types):
            cv.check_iterable_type('groups', out_groups, Integral)
            for group in out_groups:
                filters.append(openmc.EnergyoutFilter)
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # Construct CrossScore for requested scattering moment
        if moment != 'all':
            cv.check_type('moment', moment, Integral)
            cv.check_greater_than('moment', moment, 0, equality=True)
            cv.check_less_than(
                'moment', moment, self.legendre_order, equality=True)
            scores = [self.xs_tally.scores[moment]]
        else:
            scores = []

        # Construct a collection of the nuclides to retrieve from the xs tally
        if self.by_nuclide:
            if nuclides == 'all' or nuclides == 'sum' or nuclides == ['sum']:
                query_nuclides = self.get_nuclides()
            else:
                query_nuclides = nuclides
        else:
            query_nuclides = ['total']

        # Use tally summation if user requested the sum for all nuclides
        if nuclides == 'sum' or nuclides == ['sum']:
            xs_tally = self.xs_tally.summation(nuclides=query_nuclides)
            xs = xs_tally.get_values(scores=scores, filters=filters,
                                     filter_bins=filter_bins, value=value)
        else:
            xs = self.xs_tally.get_values(scores=scores, filters=filters,
                                          filter_bins=filter_bins,
                                          nuclides=query_nuclides, value=value)

        # Divide by atom number densities for microscopic cross sections
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            if value == 'mean' or value == 'std_dev':
                xs /= densities[np.newaxis, :, np.newaxis]

        # Convert and nans to zero
        xs = np.nan_to_num(xs)

        if in_groups == 'all':
            num_in_groups = self.num_groups
        else:
            num_in_groups = len(in_groups)

        if out_groups == 'all':
            num_out_groups = self.num_groups
        else:
            num_out_groups = len(out_groups)

        # Reshape tally data array with separate axes for domain and energy
        num_subdomains = int(xs.shape[0] / (num_in_groups * num_out_groups))
        new_shape = (num_subdomains, num_in_groups, num_out_groups)
        new_shape += xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Transpose the scattering matrix if requested by user
        if row_column == 'outin':
            xs = np.swapaxes(xs, 1, 2)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[:, ::-1, ::-1, :]

        if squeeze:
            xs = np.squeeze(xs)
            xs = np.atleast_2d(xs)

        return xs

    def get_pandas_dataframe(self, groups='all', nuclides='all', moment='all',
                             xs_type='macro', distribcell_paths=True):
        """Build a Pandas DataFrame for the MGXS data.

        This method leverages :meth:`openmc.Tally.get_pandas_dataframe`, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will include the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            include the cross sections summed over all nuclides. Defaults
            to 'all'.
        moment : int or 'all'
            The scattering matrix moment to return. All moments will be
             returned if the moment is 'all' (default); otherwise, a specific
             moment will be returned.
        xs_type: {'macro', 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        distribcell_paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into a
            Multi-index column with a geometric "path" to each distribcell
            instance.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        df = super(ScatterMatrixXS, self).get_pandas_dataframe(
            groups, nuclides, xs_type, distribcell_paths)

        # Add a moment column to dataframe
        if self.legendre_order > 0:
            # Insert a column corresponding to the Legendre moments
            moments = ['P{}'.format(i) for i in range(self.legendre_order+1)]
            moments = np.tile(moments, int(df.shape[0] / len(moments)))
            df['moment'] = moments

            # Place the moment column before the mean column
            columns = df.columns.tolist()
            mean_index = [i for i, s in enumerate(columns) if 'mean' in s][0]
            if self.domain_type == 'mesh':
                df = df[columns[:mean_index] + [('moment', '')] + columns[mean_index:-1]]
            else:
                df = df[columns[:mean_index] + ['moment'] + columns[mean_index:-1]]

        # Select rows corresponding to requested scattering moment
        if moment != 'all':
            cv.check_type('moment', moment, Integral)
            cv.check_greater_than('moment', moment, 0, equality=True)
            cv.check_less_than(
                'moment', moment, self.legendre_order, equality=True)
            df = df[df['moment'] == 'P{}'.format(moment)]

        return df

    def print_xs(self, subdomains='all', nuclides='all',
                 xs_type='macro', moment=0):
        """Prints a string representation for the multi-group cross section.

        Parameters
        ----------
        subdomains : Iterable of Integral or 'all'
            The subdomain IDs of the cross sections to include in the report.
            Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the report. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will report the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will report
            the cross sections summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        moment : int
            The scattering moment to print (default is 0)

        """

        # Construct a collection of the subdomains to report
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral)
        elif self.domain_type == 'distribcell':
            subdomains = np.arange(self.num_subdomains, dtype=np.int)
        elif self.domain_type == 'mesh':
            xyz = [range(1, x+1) for x in self.domain.dimension]
            subdomains = list(itertools.product(*xyz))
        else:
            subdomains = [self.domain.id]

        # Construct a collection of the nuclides to report
        if self.by_nuclide:
            if nuclides == 'all':
                nuclides = self.get_nuclides()
            if nuclides == 'sum':
                nuclides = ['sum']
            else:
                cv.check_iterable_type('nuclides', nuclides, string_types)
        else:
            nuclides = ['sum']

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        if self.correction != 'P0':
            rxn_type = '{0} (P{1})'.format(self.rxn_type, moment)
        else:
            rxn_type = self.rxn_type

        # Build header for string with type and domain info
        string = 'Multi-Group XS\n'
        string += '{0: <16}=\t{1}\n'.format('\tReaction Type', rxn_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain Type', self.domain_type)
        string += '{0: <16}=\t{1}\n'.format('\tDomain ID', self.domain.id)

        # Generate the header for an individual XS
        xs_header = '\tCross Sections [{0}]:'.format(self.get_units(xs_type))

        # If cross section data has not been computed, only print string header
        if self.tallies is None:
            print(string)
            return

        string += '{0: <16}\n'.format('\tEnergy Groups:')
        template = '{0: <12}Group {1} [{2: <10} - {3: <10}eV]\n'

        # Loop over energy groups ranges
        for group in range(1, self.num_groups + 1):
            bounds = self.energy_groups.get_group_bounds(group)
            string += template.format('', group, bounds[0], bounds[1])

        # Loop over all subdomains
        for subdomain in subdomains:

            if self.domain_type == 'distribcell':
                string += \
                    '{0: <16}=\t{1}\n'.format('\tSubdomain', subdomain)

            # Loop over all Nuclides
            for nuclide in nuclides:

                # Build header for nuclide type
                if xs_type != 'sum':
                    string += '{0: <16}=\t{1}\n'.format('\tNuclide', nuclide)

                # Build header for cross section type
                string += '{0: <16}\n'.format(xs_header)
                template = '{0: <12}Group {1} -> Group {2}:\t\t'

                # Loop over incoming/outgoing energy groups ranges
                for in_group in range(1, self.num_groups + 1):
                    for out_group in range(1, self.num_groups + 1):
                        string += template.format('', in_group, out_group)
                        average = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide], moment=moment,
                                        xs_type=xs_type, value='mean')
                        rel_err = \
                            self.get_xs([in_group], [out_group],
                                        [subdomain], [nuclide], moment=moment,
                                        xs_type=xs_type, value='rel_err')
                        average = average.flatten()[0]
                        rel_err = rel_err.flatten()[0] * 100.
                        string += '{:1.2e} +/- {:1.2e}%'.format(average,
                                                                rel_err)
                        string += '\n'
                    string += '\n'
                string += '\n'
            string += '\n'

        print(string)


class NuScatterMatrixXS(ScatterMatrixXS):
    """A scattering production matrix multi-group cross section for one or
    more Legendre moments.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`NuScatterMatrixXS.energy_groups` and
    :attr:`NuScatterMatrixXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`NuScatterMatrixXS.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`NuScatterMatrixXS.xs_tally` property.

    The calculation of the scattering-production matrix is the same as that for
    :class:`ScatterMatrixXS` except that the scattering multiplicity is
    accounted for.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    correction : 'P0' or None
        Apply the P0 correction to scattering matrices if set to 'P0'
    legendre_order : int
        The highest legendre moment in the scattering matrix (default is 0)
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`NuScatterMatrixXS.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuScatterMatrixXS, self).__init__(domain, domain_type,
                                                groups, by_nuclide, name)
        self._rxn_type = 'nu-scatter'
        self._hdf5_key = 'nu-scatter matrix'


class MultiplicityMatrixXS(MatrixMGXS):
    r"""The scattering multiplicity matrix.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`MultiplicityMatrixXS.energy_groups` and
    :attr:`MultiplicityMatrixXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`MultiplicityMatrixXS.tallies` property, which
    can then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`MultiplicityMatrixXS.xs_tally`
    property.

    For a spatial domain :math:`V`, incoming energy group
    :math:`[E_{g'},E_{g'-1}]`, and outgoing energy group :math:`[E_g,E_{g-1}]`,
    the multiplicity is calculated as:

    .. math::

       \langle \upsilon \sigma_{s,g'\rightarrow g} \phi \rangle &= \int_{r \in
       D} dr \int_{4\pi} d\Omega' \int_{E_{g'}}^{E_{g'-1}} dE' \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \sum_i \upsilon_i \sigma_i (r, E' \rightarrow
       E, \Omega' \cdot \Omega) \psi(r, E', \Omega') \\
       \langle \sigma_{s,g'\rightarrow g} \phi \rangle &= \int_{r \in
       D} dr \int_{4\pi} d\Omega' \int_{E_{g'}}^{E_{g'-1}} dE' \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \sum_i \upsilon_i \sigma_i (r, E' \rightarrow
       E, \Omega' \cdot \Omega) \psi(r, E', \Omega') \\
       \upsilon_{g'\rightarrow g} &= \frac{\langle \upsilon
       \sigma_{s,g'\rightarrow g} \rangle}{\langle \sigma_{s,g'\rightarrow g}
       \rangle}

    where :math:`\upsilon_i` is the multiplicity for the :math:`i`-th reaction.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`MultiplicityMatrixXS.tally_keys`
        property and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(MultiplicityMatrixXS, self).__init__(domain, domain_type, groups,
                                                   by_nuclide, name)
        self._rxn_type = 'multiplicity matrix'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']

    @property
    def scores(self):
        scores = ['nu-scatter', 'scatter']
        return scores

    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energy = openmc.EnergyFilter(group_edges)
        energyout = openmc.EnergyoutFilter(group_edges)

        return [[energy, energyout], [energy, energyout]]

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['nu-scatter']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

    @property
    def xs_tally(self):

        if self._xs_tally is None:
            scatter = self.tallies['scatter']

            # Compute the multiplicity
            self._xs_tally = self.rxn_rate_tally / scatter
            super(MultiplicityMatrixXS, self)._compute_xs()

        return self._xs_tally


class NuFissionMatrixXS(MatrixMGXS):
    r"""A fission production matrix multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`NuFissionMatrixXS.energy_groups` and
    :attr:`NuFissionMatrixXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`NuFissionMatrixXS.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`NuFissionMatrixXS.xs_tally` property.

    For a spatial domain :math:`V`, incoming energy group
    :math:`[E_{g'},E_{g'-1}]`, and outgoing energy group :math:`[E_g,E_{g-1}]`,
    the fission production is calculated as:

    .. math::

       \langle \nu\sigma_{f,g'\rightarrow g} \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega' \int_{E_{g'}}^{E_{g'-1}} dE' \int_{E_g}^{E_{g-1}} dE
       \; \chi(E) \nu\sigma_f (r, E') \psi(r, E', \Omega')\\
       \langle \phi \rangle &= \int_{r \in V} dr \int_{4\pi} d\Omega
       \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega) \\
       \nu\sigma_{f,g'\rightarrow g} &= \frac{\langle \nu\sigma_{f,g'\rightarrow
       g} \phi \rangle}{\langle \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`NuFissionMatrixXS.tally_keys`
        property and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(NuFissionMatrixXS, self).__init__(domain, domain_type,
                                                groups, by_nuclide, name)
        self._rxn_type = 'nu-fission'
        self._hdf5_key = 'nu-fission matrix'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']


class Chi(MGXS):
    r"""The fission spectrum.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`Chi.energy_groups` and
    :attr:`Chi.domain` properties. Tallies for the flux and appropriate reaction
    rates over the specified domain are generated automatically via the
    :attr:`Chi.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`Chi.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    fission spectrum is calculated as:

    .. math::

       \langle \nu\sigma_{f,g' \rightarrow g} \phi \rangle &= \int_{r \in V} dr
       \int_{4\pi} d\Omega' \int_0^\infty dE' \int_{E_g}^{E_{g-1}} dE \; \chi(E)
       \nu\sigma_f (r, E') \psi(r, E', \Omega')\\
       \langle \nu\sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu\sigma_f (r,
       E') \psi(r, E', \Omega') \\
       \chi_g &= \frac{\langle \nu\sigma_{f,g' \rightarrow g} \phi \rangle}
       {\langle \nu\sigma_f \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`Chi.tally_keys` property and values are
        instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U238', 'O16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(Chi, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'chi'
        self._estimator = 'analog'
        self._valid_estimators = ['analog']

    @property
    def scores(self):
        return ['nu-fission', 'nu-fission']

    @property
    def filters(self):
        # Create the non-domain specific Filters for the Tallies
        group_edges = self.energy_groups.group_edges
        energyout = openmc.EnergyoutFilter(group_edges)
        energyin = openmc.EnergyFilter([group_edges[0], group_edges[-1]])
        return [[energyin], [energyout]]

    @property
    def tally_keys(self):
        return ['nu-fission-in', 'nu-fission-out']

    @property
    def rxn_rate_tally(self):
        if self._rxn_rate_tally is None:
            self._rxn_rate_tally = self.tallies['nu-fission-out']
            self._rxn_rate_tally.sparse = self.sparse
        return self._rxn_rate_tally

    @property
    def xs_tally(self):

        if self._xs_tally is None:
            nu_fission_in = self.tallies['nu-fission-in']

            # Remove coarse energy filter to keep it out of tally arithmetic
            energy_filter = nu_fission_in.find_filter(openmc.EnergyFilter)
            nu_fission_in.remove_filter(energy_filter)

            # Compute chi
            self._xs_tally = self.rxn_rate_tally / nu_fission_in

            # Add the coarse energy filter back to the nu-fission tally
            nu_fission_in.filters.append(energy_filter)

        return self._xs_tally

    def get_slice(self, nuclides=[], groups=[]):
        """Build a sliced Chi for the specified nuclides and energy groups.

        This method constructs a new MGXS to encapsulate a subset of the data
        represented by this MGXS. The subset of data to include in the tally
        slice is determined by the nuclides and energy groups specified in
        the input parameters.

        Parameters
        ----------
        nuclides : list of str
            A list of nuclide name strings
            (e.g., ['U235', 'U238']; default is [])
        groups : list of Integral
            A list of energy group indices starting at 1 for the high energies
            (e.g., [1, 2, 3]; default is [])

        Returns
        -------
        openmc.mgxs.MGXS
            A new MGXS which encapsulates the subset of data requested
            for the nuclide(s) and/or energy group(s) requested in the
            parameters.

        """

        # Temporarily remove energy filter from nu-fission-in since its
        # group structure will work in super MGXS.get_slice(...) method
        nu_fission_in = self.tallies['nu-fission-in']
        energy_filter = nu_fission_in.find_filter(openmc.EnergyFilter)
        nu_fission_in.remove_filter(energy_filter)

        # Call super class method and null out derived tallies
        slice_xs = super(Chi, self).get_slice(nuclides, groups)
        slice_xs._rxn_rate_tally = None
        slice_xs._xs_tally = None

        # Slice energy groups if needed
        if len(groups) != 0:
            filter_bins = []
            for group in groups:
                group_bounds = self.energy_groups.get_group_bounds(group)
                filter_bins.append(group_bounds)
            filter_bins = [tuple(filter_bins)]

            # Slice nu-fission-out tally along energyout filter
            nu_fission_out = slice_xs.tallies['nu-fission-out']
            tally_slice = nu_fission_out.get_slice(
                filters=[openmc.EnergyoutFilter], filter_bins=filter_bins)
            slice_xs._tallies['nu-fission-out'] = tally_slice

        # Add energy filter back to nu-fission-in tallies
        self.tallies['nu-fission-in'].add_filter(energy_filter)
        slice_xs._tallies['nu-fission-in'].add_filter(energy_filter)

        slice_xs.sparse = self.sparse
        return slice_xs

    def merge(self, other):
        """Merge another Chi with this one

        If results have been loaded from a statepoint, then Chi are only
        mergeable along one and only one of energy groups or nuclides.

        Parameters
        ----------
        other : openmc.mgxs.MGXS
            MGXS to merge with this one

        Returns
        -------
        merged_mgxs : openmc.mgxs.MGXS
            Merged MGXS
        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge a Chi MGXS')

        # Create deep copy of tally to return as merged tally
        merged_mgxs = copy.deepcopy(self)
        merged_mgxs._derived = True
        merged_mgxs._rxn_rate_tally = None
        merged_mgxs._xs_tally = None

        # Merge energy groups
        if self.energy_groups != other.energy_groups:
            merged_groups = self.energy_groups.merge(other.energy_groups)
            merged_mgxs.energy_groups = merged_groups

        # Merge nuclides
        if self.nuclides != other.nuclides:

            # The nuclides must be mutually exclusive
            for nuclide in self.nuclides:
                if nuclide in other.nuclides:
                    msg = 'Unable to merge a Chi MGXS with shared nuclides'
                    raise ValueError(msg)

            # Concatenate lists of nuclides for the merged MGXS
            merged_mgxs.nuclides = self.nuclides + other.nuclides

        # Merge tallies
        for tally_key in self.tallies:
            merged_tally = self.tallies[tally_key].merge(other.tallies[tally_key])
            merged_mgxs.tallies[tally_key] = merged_tally

        return merged_mgxs

    def get_xs(self, groups='all', subdomains='all', nuclides='all',
               xs_type='macro', order_groups='increasing',
               value='mean', squeeze=True, **kwargs):
        """Returns an array of the fission spectrum.

        This method constructs a 3D NumPy array for the requested
        multi-group cross section data for one or more subdomains
        (1st dimension), energy groups (2nd dimension), and nuclides
        (3rd dimension).

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        subdomains : Iterable of Integral or 'all'
            Subdomain IDs of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            A list of nuclide name strings (e.g., ['U235', 'U238']). The
            special string 'all' will return the cross sections for all nuclides
            in the spatial domain. The special string 'sum' will return the
            cross section summed over all nuclides. Defaults to 'all'.
        xs_type: {'macro', 'micro'}
            This parameter is not relevant for chi but is included here to
            mirror the parent MGXS.get_xs(...) class method
        order_groups: {'increasing', 'decreasing'}
            Return the cross section indexed according to increasing or
            decreasing energy groups (decreasing or increasing energies).
            Defaults to 'increasing'.
        value : {'mean', 'std_dev', 'rel_err'}
            A string for the type of value to return. Defaults to 'mean'.
        squeeze : bool
            A boolean representing whether to eliminate the extra dimensions
            of the multi-dimensional array to be returned. Defaults to True.

        Returns
        -------
        numpy.ndarray
            A NumPy array of the multi-group cross section indexed in the order
            each group, subdomain and nuclide is listed in the parameters.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        cv.check_value('value', value, ['mean', 'std_dev', 'rel_err'])
        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # FIXME: Unable to get microscopic xs for mesh domain because the mesh
        # cells do not know the nuclide densities in each mesh cell.
        if self.domain_type == 'mesh' and xs_type == 'micro':
            msg = 'Unable to get micro xs for mesh domain since the mesh ' \
                  'cells do not know the nuclide densities in each mesh cell.'
            raise ValueError(msg)

        filters = []
        filter_bins = []

        # Construct a collection of the domain filter bins
        if not isinstance(subdomains, string_types):
            cv.check_iterable_type('subdomains', subdomains, Integral, max_depth=3)
            for subdomain in subdomains:
                filters.append(_DOMAIN_TO_FILTER[self.domain_type])
                filter_bins.append((subdomain,))

        # Construct list of energy group bounds tuples for all requested groups
        if not isinstance(groups, string_types):
            cv.check_iterable_type('groups', groups, Integral)
            for group in groups:
                filters.append(openmc.EnergyoutFilter)
                filter_bins.append((self.energy_groups.get_group_bounds(group),))

        # If chi was computed for each nuclide in the domain
        if self.by_nuclide:

            # Get the sum as the fission source weighted average chi for all
            # nuclides in the domain
            if nuclides == 'sum' or nuclides == ['sum']:

                # Retrieve the fission production tallies
                nu_fission_in = self.tallies['nu-fission-in']
                nu_fission_out = self.tallies['nu-fission-out']

                # Sum out all nuclides
                nuclides = self.get_nuclides()
                nu_fission_in = nu_fission_in.summation(nuclides=nuclides)
                nu_fission_out = nu_fission_out.summation(nuclides=nuclides)

                # Remove coarse energy filter to keep it out of tally arithmetic
                energy_filter = nu_fission_in.find_filter(openmc.EnergyFilter)
                nu_fission_in.remove_filter(energy_filter)

                # Compute chi and store it as the xs_tally attribute so we can
                # use the generic get_xs(...) method
                xs_tally = nu_fission_out / nu_fission_in

                # Add the coarse energy filter back to the nu-fission tally
                nu_fission_in.filters.append(energy_filter)

                xs = xs_tally.get_values(filters=filters,
                                         filter_bins=filter_bins, value=value)

            # Get chi for all nuclides in the domain
            elif nuclides == 'all':
                nuclides = self.get_nuclides()
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

            # Get chi for user-specified nuclides in the domain
            else:
                cv.check_iterable_type('nuclides', nuclides, string_types)
                xs = self.xs_tally.get_values(filters=filters,
                                              filter_bins=filter_bins,
                                              nuclides=nuclides, value=value)

        # If chi was computed as an average of nuclides in the domain
        else:
            xs = self.xs_tally.get_values(filters=filters,
                                          filter_bins=filter_bins, value=value)

        # Eliminate the trivial score dimension
        xs = np.squeeze(xs, axis=len(xs.shape) - 1)
        xs = np.nan_to_num(xs)

        # Reshape tally data array with separate axes for domain and energy
        if groups == 'all':
            num_groups = self.num_groups
        else:
            num_groups = len(groups)

        num_subdomains = int(xs.shape[0] / num_groups)
        new_shape = (num_subdomains, num_groups) + xs.shape[1:]
        xs = np.reshape(xs, new_shape)

        # Reverse data if user requested increasing energy groups since
        # tally data is stored in order of increasing energies
        if order_groups == 'increasing':
            xs = xs[:, ::-1, :]

        if squeeze:
            xs = np.squeeze(xs)
            xs = np.atleast_1d(xs)

        return xs

    def get_pandas_dataframe(self, groups='all', nuclides='all',
                             xs_type='macro', distribcell_paths=False):
        """Build a Pandas DataFrame for the MGXS data.

        This method leverages :meth:`openmc.Tally.get_pandas_dataframe`, but
        renames the columns with terminology appropriate for cross section data.

        Parameters
        ----------
        groups : Iterable of Integral or 'all'
            Energy groups of interest. Defaults to 'all'.
        nuclides : Iterable of str or 'all' or 'sum'
            The nuclides of the cross-sections to include in the dataframe. This
            may be a list of nuclide name strings (e.g., ['U235', 'U238']).
            The special string 'all' will include the cross sections for all
            nuclides in the spatial domain. The special string 'sum' will
            include the cross sections summed over all nuclides. Defaults to
            'all'.
        xs_type: {'macro', 'micro'}
            Return macro or micro cross section in units of cm^-1 or barns.
            Defaults to 'macro'.
        distribcell_paths : bool, optional
            Construct columns for distribcell tally filters (default is True).
            The geometric information in the Summary object is embedded into
            a Multi-index column with a geometric "path" to each distribcell
            instance.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame for the cross section data.

        Raises
        ------
        ValueError
            When this method is called before the multi-group cross section is
            computed from tally data.

        """

        # Build the dataframe using the parent class method
        df = super(Chi, self).get_pandas_dataframe(
            groups, nuclides, xs_type, distribcell_paths=distribcell_paths)

        # If user requested micro cross sections, multiply by the atom
        # densities to cancel out division made by the parent class method
        if xs_type == 'micro':
            if self.by_nuclide:
                densities = self.get_nuclide_densities(nuclides)
            else:
                densities = self.get_nuclide_densities('sum')
            tile_factor = df.shape[0] / len(densities)
            df['mean'] *= np.tile(densities, tile_factor)
            df['std. dev.'] *= np.tile(densities, tile_factor)

        return df

    def get_units(self, xs_type='macro'):
        """Returns the units of Chi.

        This method returns the units of Chi, which is "%" for both macro
        and micro xs types.

        Parameters
        ----------
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section units.
            Defaults to 'macro'.

        Returns
        -------
        str
            A string representing the units of Chi.

        """

        cv.check_value('xs_type', xs_type, ['macro', 'micro'])

        # Chi has the same units (%) for both macro and micro
        return '%'


class ChiPrompt(Chi):
    r"""The prompt fission spectrum.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`ChiPrompt.energy_groups` and
    :attr:`ChiPrompt.domain` properties. Tallies for the flux and appropriate
    reaction rates over the specified domain are generated automatically via the
    :attr:`ChiPrompt.tallies` property, which can then be appended to a
    :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`ChiPrompt.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    fission spectrum is calculated as:

    .. math::

       \langle \nu^p \sigma_{f,g' \rightarrow g} \phi \rangle &= \int_{r \in V}
       dr \int_{4\pi} d\Omega' \int_0^\infty dE' \int_{E_g}^{E_{g-1}} dE \;
       \chi(E)^p \nu^p \sigma_f (r, E') \psi(r, E', \Omega')\\
       \langle \nu^p \sigma_f \phi \rangle &= \int_{r \in V} dr \int_{4\pi}
       d\Omega' \int_0^\infty dE' \int_0^\infty dE \; \chi(E) \nu^p \sigma_f (r,
       E') \psi(r, E', \Omega') \\
       \chi_g^p &= \frac{\langle \nu^p \sigma_{f,g' \rightarrow g} \phi \rangle}
       {\langle \nu^p \sigma_f \phi \rangle}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : 'analog'
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`ChiPrompt.tally_keys` property and
        values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(ChiPrompt, self).__init__(domain, domain_type, groups, by_nuclide, name)
        self._rxn_type = 'chi-prompt'

    @property
    def scores(self):
        return ['prompt-nu-fission', 'prompt-nu-fission']


class InverseVelocity(MGXS):
    r"""An inverse velocity multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group neutron inverse velocities for multi-group neutronics
    calculations. The units of inverse velocity are seconds per centimeter. At a
    minimum, one needs to set the :attr:`InverseVelocity.energy_groups` and
    :attr:`InverseVelocity.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`InverseVelocity.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`InverseVelocity.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    neutron inverse velocities are calculated by tallying the flux-weighted
    inverse velocity and the flux. The inverse velocity is then the
    flux-weighted inverse velocity divided by the flux:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \frac{\psi (r, E, \Omega)}{v (r, E)}}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`InverseVelocity.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file) and the number of mesh cells for
        'mesh' domain types.
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(InverseVelocity, self).__init__(domain, domain_type,
                                              groups, by_nuclide, name)
        self._rxn_type = 'inverse-velocity'

    def get_units(self, xs_type='macro'):
        """Returns the units of InverseVelocity.

        This method returns the units of an InverseVelocity based on a desired
        xs_type.

        Parameters
        ----------
        xs_type: {'macro', 'micro'}
            Return the macro or micro cross section units.
            Defaults to 'macro'.

        Returns
        -------
        str
            A string representing the units of the InverseVelocity.

        """

        if xs_type == 'macro':
            return 'second/cm'
        else:
            raise ValueError('Unable to return the units of InverseVelocity for'
                             ' xs_type other than "macro"')


class PromptNuFissionXS(MGXS):
    r"""A prompt fission neutron production multi-group cross section.

    This class can be used for both OpenMC input generation and tally data
    post-processing to compute spatially-homogenized and energy-integrated
    multi-group cross sections for multi-group neutronics calculations. At a
    minimum, one needs to set the :attr:`PromptNuFissionXS.energy_groups` and
    :attr:`PromptNuFissionXS.domain` properties. Tallies for the flux and
    appropriate reaction rates over the specified domain are generated
    automatically via the :attr:`PromptNuFissionXS.tallies` property, which can
    then be appended to a :class:`openmc.Tallies` instance.

    For post-processing, the :meth:`MGXS.load_from_statepoint` will pull in the
    necessary data to compute multi-group cross sections from a
    :class:`openmc.StatePoint` instance. The derived multi-group cross section
    can then be obtained from the :attr:`PromptNuFissionXS.xs_tally` property.

    For a spatial domain :math:`V` and energy group :math:`[E_g,E_{g-1}]`, the
    fission spectrum is calculated as:

    .. math::

       \frac{\int_{r \in V} dr \int_{4\pi} d\Omega \int_{E_g}^{E_{g-1}} dE \;
       \nu\sigma_f^p (r, E) \psi (r, E, \Omega)}{\int_{r \in V} dr \int_{4\pi}
       d\Omega \int_{E_g}^{E_{g-1}} dE \; \psi (r, E, \Omega)}.

    Parameters
    ----------
    domain : openmc.Material or openmc.Cell or openmc.Universe or openmc.Mesh
        The domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        The domain type for spatial homogenization
    groups : openmc.mgxs.EnergyGroups
        The energy group structure for energy condensation
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    name : str, optional
        Name of the multi-group cross section. Used as a label to identify
        tallies in OpenMC 'tallies.xml' file.

    Attributes
    ----------
    name : str, optional
        Name of the multi-group cross section
    rxn_type : str
        Reaction type (e.g., 'total', 'nu-fission', etc.)
    by_nuclide : bool
        If true, computes cross sections for each nuclide in domain
    domain : Material or Cell or Universe or Mesh
        Domain for spatial homogenization
    domain_type : {'material', 'cell', 'distribcell', 'universe', 'mesh'}
        Domain type for spatial homogenization
    energy_groups : openmc.mgxs.EnergyGroups
        Energy group structure for energy condensation
    tally_trigger : openmc.Trigger
        An (optional) tally precision trigger given to each tally used to
        compute the cross section
    scores : list of str
        The scores in each tally used to compute the multi-group cross section
    filters : list of openmc.Filter
        The filters in each tally used to compute the multi-group cross section
    tally_keys : list of str
        The keys into the tallies dictionary for each tally used to compute
        the multi-group cross section
    estimator : {'tracklength', 'collision', 'analog'}
        The tally estimator used to compute the multi-group cross section
    tallies : collections.OrderedDict
        OpenMC tallies needed to compute the multi-group cross section. The keys
        are strings listed in the :attr:`PromptNuFissionXS.tally_keys` property
        and values are instances of :class:`openmc.Tally`.
    rxn_rate_tally : openmc.Tally
        Derived tally for the reaction rate tally used in the numerator to
        compute the multi-group cross section. This attribute is None
        unless the multi-group cross section has been computed.
    xs_tally : openmc.Tally
        Derived tally for the multi-group cross section. This attribute
        is None unless the multi-group cross section has been computed.
    num_subdomains : int
        The number of subdomains is unity for 'material', 'cell' and 'universe'
        domain types. This is equal to the number of cell instances
        for 'distribcell' domain types (it is equal to unity prior to loading
        tally data from a statepoint file).
    num_nuclides : int
        The number of nuclides for which the multi-group cross section is
        being tracked. This is unity if the by_nuclide attribute is False.
    nuclides : Iterable of str or 'sum'
        The optional user-specified nuclides for which to compute cross
        sections (e.g., 'U-238', 'O-16'). If by_nuclide is True but nuclides
        are not specified by the user, all nuclides in the spatial domain
        are included. This attribute is 'sum' if by_nuclide is false.
    sparse : bool
        Whether or not the MGXS' tallies use SciPy's LIL sparse matrix format
        for compressed data storage
    loaded_sp : bool
        Whether or not a statepoint file has been loaded with tally data
    derived : bool
        Whether or not the MGXS is merged from one or more other MGXS
    hdf5_key : str
        The key used to index multi-group cross sections in an HDF5 data store

    """

    def __init__(self, domain=None, domain_type=None,
                 groups=None, by_nuclide=False, name=''):
        super(PromptNuFissionXS, self).__init__(domain, domain_type, groups,
                                                by_nuclide, name)
        self._rxn_type = 'prompt-nu-fission'
