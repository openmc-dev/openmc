from collections import Iterable
from numbers import Real, Integral
import copy
import sys

import numpy as np

import openmc.checkvalue as cv


if sys.version_info[0] >= 3:
    basestring = str

class EnergyGroups(object):
    """An energy groups structure used for multi-group cross-sections.

    Parameters
    ----------
    group_edges : NumPy array
        The energy group boundaries [MeV]
    num_groups : Integral
        The number of energy groups

    Attributes
    ----------
    group_edges : NumPy array
        The energy group boundaries [MeV]
    num_groups : Integral
        The number of energy groups

    """

    def __init__(self):
        self._group_edges = None
        self._num_groups = None

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy this object, create a copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone.group_edges = copy.deepcopy(self._group_edges, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    @property
    def group_edges(self):
        return self._group_edges

    @property
    def num_groups(self):
        return self._num_groups

    @group_edges.setter
    def group_edges(self, edges):
        cv.check_type('group edges', edges, Iterable, Integral)
        cv.check_length('number of group edges', edges, 2)
        self._group_edges = np.array(edges)
        self._num_groups = len(edges)-1

    def __eq__(self, other):
        if not isinstance(other, EnergyGroups):
            return False
        elif self._group_edges != other._group_edges:
            return False

    def generate_bin_edges(self, start, stop, num_groups, type='linear'):
        """Generate equally or logarithmically-spaced energy group boundaries.

        Parameters
        ----------
        start : Real
            The lowest energy in MeV
        stop : Real
            The highest energy in MeV
        num_groups : Integral
            The number of energy groups
        type : str
            The spacing between groups ('linear' or 'logarithmic')

        """

        cv.check_type('first edge', start, Real)
        cv.check_type('last edge', stop, Real)
        cv.check_type('number of groups', num_groups, Integral)
        cv.check_type('type', type, basestring)
        cv.check_greater_than('first edge', start, 0, True)
        cv.check_greater_than('first edge', stop, start, False)
        cv.check_greater_than('number of groups', num_groups, 0)
        cv.check_value('type', type, ('linear', 'logarithmic'))

        if type == 'linear':
            self.group_edges = np.linspace(start, stop, num_groups+1)
        elif type == 'logarithmic':
            self.group_edges = \
                np.logspace(np.log10(start), np.log10(stop), num_groups+1)

        self._num_groups = num_groups

    def get_group(self, energy):
        """Returns the energy group in which the given energy resides.

        Parameters
        ----------
        energy : Real
            The energy of interest in MeV

        Returns
        -------
        Integral
            The energy group index, starting at 1 for the highest energies

        Raises
        ------
        ValueError
            If the group edges have not yet been set.

        """

        if self.group_edges is None:
            msg = 'Unable to get energy group for energy "{0}" eV since ' \
                  'the group edges have not yet been set'.format(energy)
            raise ValueError(msg)

        index = np.where(self.group_edges > energy)[0]
        group = self.num_groups - index
        return group

    def get_group_bounds(self, group):
        """Returns the energy boundaries for the energy group of interest.

        Parameters
        ----------
        group : Integral
            The energy group index, starting at 1 for the highest energies

        Returns
        -------
        2-tuple
            The low and high energy bounds for the group in MeV

        Raises
        ------
        ValueError
            If the group edges have not yet been set.

        """

        if self.group_edges is None:
            msg = 'Unable to get energy group bounds for group "{0}" since ' \
                  'the group edges have not yet been set'.format(group)
            raise ValueError(msg)

        lower = self.group_edges[self.num_groups-group]
        upper = self.group_edges[self.num_groups-group+1]
        return (lower, upper)

    def get_group_indices(self, groups='all'):
        """Returns the array indices for one or more energy groups.

        Parameters
        ----------
        groups : str, tuple
            The energy groups of interest - a tuple of the energy group indices,
            starting at 1 for the highest energies (default is 'all')

        Returns
        -------
        NumPy.ndarray
            The NumPy array indices for each energy group of interest

        Raises
        ------
        ValueError
            If the group edges have not yet been set, or if a group is requested
            that is outside the bounds of the number of energy groups.

        """

        if self.group_edges is None:
            msg = 'Unable to get energy group indices for groups "{0}" since ' \
                  'the group edges have not yet been set'.format(groups)
            raise ValueError(msg)

        if groups == 'all':
            indices = np.arange(self.num_groups)
        else:
            indices = np.zeros(len(groups), dtype=np.int64)

        for i, group in enumerate(groups):
            cv.check_greater_than('group', group, 0)
            cv.check_less_than('group', group, self.num_groups, True)
            indices[i] = group - 1

        return indices

    def get_condensed_groups(self, coarse_groups):
        """Return a coarsened version of this EnergyGroups object.

        This method merges together energy groups in this object into wider
        energy groups as defined by the list of groups specified by the user,
        and returns a new, coarse EnergyGroups object.

        Parameters
        ----------
        coarse_groups : list
            The energy groups of interest - a list of 2-tuples, each directly
             corresponding to one of the new coarse groups. The values in the
             2-tuples are upper/lower energy groups used to construct a new
             coarse group.

        Returns
        -------
        EnergyGroups
            A coarsened version of this EnergyGroups object.

        Raises
        ------
        ValueError
            If the group edges have not yet been set.
        """

        cv.check_type('group edges', coarse_groups, Iterable)
        for group in coarse_groups:
            cv.check_value('group edges', group, Iterable)
            cv.check_length('group edges', group, 2)
            cv.check_greater_than('lower group', group[0], 1, True)
            cv.check_less_than('lower group', group[0], self.num_groups, True)
            cv.check_greater_than('upper group', group[0], 1, True)
            cv.check_less_than('upper group', group[0], self.num_groups, True)
            cv.check_less_than('lower group', group[0], group[1], False)

        # Compute the group indices into the coarse group
        group_bounds = list()
        for group in coarse_groups:
            group_bounds.append(group[0])
        group_bounds.append(coarse_groups[-1][1])

        # Determine the indices mapping the fine-to-coarse energy groups
        group_bounds = np.asarray(group_bounds)
        group_indices = np.flipud(self.num_groups - group_bounds)
        group_indices[-1] += 1

        # Determine the edges between coarse energy groups and sort
        # in increasing order in case the user passed in unordered groups
        group_edges = self.group_edges[group_indices]
        group_edges = np.sort(group_edges)

        # Create a new condensed EnergyGroups object
        condensed_groups = EnergyGroups()
        condensed_groups.group_edges = group_edges
        return condensed_groups