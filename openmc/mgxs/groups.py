from collections import Iterable
from numbers import Real
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
    group_edges : Iterable of Real
        The energy group boundaries [MeV]

    Attributes
    ----------
    group_edges : Iterable of Real
        The energy group boundaries [MeV]
    num_groups : Integral
        The number of energy groups

    """

    def __init__(self, group_edges=None):
        self._group_edges = None

        if group_edges is not None:
            self.group_edges = group_edges

    def __deepcopy__(self, memo):
        existing = memo.get(id(self))

        # If this is the first time we have tried to copy object, create copy
        if existing is None:
            clone = type(self).__new__(type(self))
            clone._group_edges = copy.deepcopy(self.group_edges, memo)

            memo[id(self)] = clone

            return clone

        # If this object has been copied before, return the first copy made
        else:
            return existing

    def __eq__(self, other):
        if not isinstance(other, EnergyGroups):
            return False
        elif (self.group_edges != other.group_edges).all():
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(tuple(self.group_edges))

    @property
    def group_edges(self):
        return self._group_edges

    @property
    def num_groups(self):
        return len(self.group_edges) - 1

    @group_edges.setter
    def group_edges(self, edges):
        cv.check_type('group edges', edges, Iterable, Real)
        cv.check_greater_than('number of group edges', len(edges), 1)
        self._group_edges = np.array(edges)

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
            msg = 'Unable to get energy group for energy "{0}" MeV since ' \
                  'the group edges have not yet been set'.format(energy)
            raise ValueError(msg)

        index = np.where(self.group_edges > energy)[0][0]
        group = self.num_groups - index + 1
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

        cv.check_greater_than('group', group, 0)
        cv.check_less_than('group', group, self.num_groups, equality=True)

        lower = self.group_edges[self.num_groups-group]
        upper = self.group_edges[self.num_groups-group+1]
        return lower, upper

    def get_group_indices(self, groups='all'):
        """Returns the array indices for one or more energy groups.

        Parameters
        ----------
        groups : str, tuple
            The energy groups of interest - a tuple of the energy group indices,
            starting at 1 for the highest energies (default is 'all')

        Returns
        -------
        ndarray
            The ndarray array indices for each energy group of interest

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
            return np.arange(self.num_groups)
        else:
            indices = np.zeros(len(groups), dtype=np.int)

        for i, group in enumerate(groups):
            cv.check_greater_than('group', group, 0)
            cv.check_less_than('group', group, self.num_groups, equality=True)
            indices[i] = group - 1

        return indices

    def get_condensed_groups(self, coarse_groups):
        """Return a coarsened version of this EnergyGroups object.

        This method merges together energy groups in this object into wider
        energy groups as defined by the list of groups specified by the user,
        and returns a new, coarse EnergyGroups object.

        Parameters
        ----------
        coarse_groups : Iterable of 2-tuple
            The energy groups of interest - a list of 2-tuples, each directly
            corresponding to one of the new coarse groups. The values in the
            2-tuples are upper/lower energy groups used to construct a new
            coarse group. For example, if [(1,2), (3,4)] was used as the coarse
            groups, fine groups 1 and 2 would be merged into coarse group 1
            while fine groups 3 and 4 would be merged into coarse group 2.

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
            cv.check_type('group edges', group, Iterable)
            cv.check_length('group edges', group, 2)
            cv.check_greater_than('lower group', group[0], 1, True)
            cv.check_less_than('lower group', group[0], self.num_groups, True)
            cv.check_greater_than('upper group', group[0], 1, True)
            cv.check_less_than('upper group', group[0], self.num_groups, True)
            cv.check_less_than('lower group', group[0], group[1], False)

        # Compute the group indices into the coarse group
        group_bounds = [group[1] for group in coarse_groups]
        group_bounds.insert(0, coarse_groups[0][0])

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
