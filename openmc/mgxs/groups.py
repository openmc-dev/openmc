from collections.abc import Iterable
import copy
from numbers import Real

import numpy as np

import openmc.checkvalue as cv
import openmc.mgxs


class EnergyGroups:
    """An energy group structure used for multigroup cross-sections.

    Parameters
    ----------
    group_edges : Iterable of float or str
        The energy group boundaries in [eV] or the name of the group structure
        (Must be a valid key in the openmc.mgxs.GROUP_STRUCTURES dictionary).

        .. versionchanged:: 0.14.0
            Changed to allow a string specifying the group structure name.

    Attributes
    ----------
    group_edges : np.ndarray
        The energy group boundaries in [eV]
    num_groups : int
        The number of energy groups

    """

    def __init__(self, group_edges):
        if isinstance(group_edges, str):
            self._name = group_edges.upper()
            group_edges = openmc.mgxs.GROUP_STRUCTURES[self._name]

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
        elif self.num_groups != other.num_groups:
            return False
        elif np.allclose(self.group_edges, other.group_edges):
            return True
        else:
            return False

    def __hash__(self):
        return hash(tuple(self.group_edges))

    def __repr__(self):
        if hasattr(self, '_name'):
            return f"<EnergyGroups: {self.num_groups} groups ({self._name})>"
        else:
            return f"<EnergyGroups: {self.num_groups} groups>"

    @property
    def group_edges(self):
        return self._group_edges

    @group_edges.setter
    def group_edges(self, edges):
        cv.check_type('group edges', edges, Iterable, Real)
        cv.check_greater_than('number of group edges', len(edges), 1)
        self._group_edges = np.array(edges)

    @property
    def num_groups(self):
        return len(self.group_edges) - 1

    def get_group(self, energy):
        """Returns the energy group in which the given energy resides.

        Parameters
        ----------
        energy : float
            The energy of interest in eV

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

        index = np.where(self.group_edges > energy)[0][0]
        group = self.num_groups - index + 1
        return group

    def get_group_bounds(self, group):
        """Returns the energy boundaries for the energy group of interest.

        Parameters
        ----------
        group : int
            The energy group index, starting at 1 for the highest energies

        Returns
        -------
        2-tuple
            The low and high energy bounds for the group in eV

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
        numpy.ndarray
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
            indices = np.zeros(len(groups), dtype=int)

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
        openmc.mgxs.EnergyGroups
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
        return EnergyGroups(group_edges)

    def can_merge(self, other):
        """Determine if energy groups can be merged with another.

        Parameters
        ----------
        other : openmc.mgxs.EnergyGroups
            EnergyGroups to compare with

        Returns
        -------
        bool
            Whether the energy groups can be merged

        """

        if not isinstance(other, EnergyGroups):
            return False

        # If the energy group structures match then groups are mergeable
        if self == other:
            return True

        # This low energy edge coincides with other's high energy edge
        if self.group_edges[0] == other.group_edges[-1]:
            return True
        # This high energy edge coincides with other's low energy edge
        elif self.group_edges[-1] == other.group_edges[0]:
            return True
        else:
            return False

    def merge(self, other):
        """Merge this energy groups with another.

        Parameters
        ----------
        other : openmc.mgxs.EnergyGroups
            EnergyGroups to merge with

        Returns
        -------
        merged_groups : openmc.mgxs.EnergyGroups
            EnergyGroups resulting from the merge

        """

        if not self.can_merge(other):
            raise ValueError('Unable to merge energy groups')

        # Create deep copy to return as merged energy groups
        merged_groups = copy.deepcopy(self)

        # Merge unique filter bins
        merged_edges = np.concatenate((self.group_edges, other.group_edges))
        merged_edges = np.unique(merged_edges)
        merged_edges = sorted(merged_edges)

        # Assign merged edges to merged groups
        merged_groups.group_edges = list(merged_edges)
        return merged_groups


def convert_flux_groups(flux, source_groups, target_groups):
    """Convert flux spectrum between energy group structures.

    Uses flux-per-unit-lethargy conservation, which assumes constant flux per
    unit lethargy within each source group and distributes flux to target
    groups proportionally to their lethargy width.

    .. versionadded:: 0.15.4

    Parameters
    ----------
    flux : Iterable of float
        Flux values for source groups. Length must equal
        source_groups.num_groups.
    source_groups : EnergyGroups or str
        Energy group structure of the input flux with boundaries in [eV].
        Can be an EnergyGroups instance or the name of a group structure
        (e.g., 'CCFE-709').
    target_groups : EnergyGroups or str
        Target energy group structure with boundaries in [eV]. Can be an
        EnergyGroups instance or the name of a group structure
        (e.g., 'UKAEA-1102').

    Returns
    -------
    numpy.ndarray
        Flux values for target groups. Total flux is conserved for
        overlapping energy regions.

    Raises
    ------
    TypeError
        If source_groups or target_groups is not EnergyGroups or str
    ValueError
        If flux length doesn't match source_groups, or flux contains
        negative, NaN, or infinite values

    See Also
    --------
    EnergyGroups : Energy group structure class

    Notes
    -----
    The assumption of constant flux per unit lethargy within each source
    group is physically reasonable for most reactor spectra but is not
    exact. For best accuracy, use source spectra with sufficiently fine
    energy resolution.

    Examples
    --------
    Convert FNS 709-group flux to UKAEA-1102 structure:

    >>> import numpy as np
    >>> flux_709 = np.load('tests/fns_flux_709.npy')
    >>> flux_1102 = openmc.mgxs.convert_flux_groups(flux_709, 'CCFE-709', 'UKAEA-1102')

    Convert using EnergyGroups instances:

    >>> source = openmc.mgxs.EnergyGroups([1.0, 10.0, 100.0])
    >>> target = openmc.mgxs.EnergyGroups([1.0, 5.0, 10.0, 50.0, 100.0])
    >>> flux_target = openmc.mgxs.convert_flux_groups([1e8, 2e8], source, target)

    References
    ----------
    .. [1] J. J. Duderstadt and L. J. Hamilton, "Nuclear Reactor Analysis,"
       John Wiley & Sons, 1976.
    .. [2] M. Fleming and J.-Ch. Sublet, "FISPACT-II User Manual,"
       UKAEA-R(18)001, UK Atomic Energy Authority, 2018. See GRPCONVERT keyword.

    """
    # Handle string group structure names
    if isinstance(source_groups, str):
        source_groups = EnergyGroups(source_groups)
    if isinstance(target_groups, str):
        target_groups = EnergyGroups(target_groups)

    # Type validation
    cv.check_type('source_groups', source_groups, EnergyGroups)
    cv.check_type('target_groups', target_groups, EnergyGroups)

    # Convert flux to numpy array
    flux = np.asarray(flux, dtype=np.float64)
    if flux.ndim != 1:
        raise ValueError(f'flux must be 1-dimensional, got shape {flux.shape}')

    # Validate flux length matches source groups
    if len(flux) != source_groups.num_groups:
        raise ValueError(
            f'Length of flux ({len(flux)}) must equal number of source '
            f'groups ({source_groups.num_groups})'
        )

    # Check for invalid flux values
    if np.any(np.isnan(flux)):
        raise ValueError('flux contains NaN values')
    if np.any(np.isinf(flux)):
        raise ValueError('flux contains infinite values')
    if np.any(flux < 0):
        raise ValueError('flux values must be non-negative')

    # Get energy edges
    source_edges = source_groups.group_edges
    target_edges = target_groups.group_edges
    num_target = target_groups.num_groups

    # Initialize output array
    flux_target = np.zeros(num_target)

    # Main conversion loop: distribute flux using lethargy weighting
    for idx_src, flux_src in enumerate(flux):
        if flux_src == 0:
            continue

        e_low_src = source_edges[idx_src]
        e_high_src = source_edges[idx_src + 1]
        lethargy_src = np.log(e_high_src / e_low_src)

        for idx_tgt in range(num_target):
            e_low_tgt = target_edges[idx_tgt]
            e_high_tgt = target_edges[idx_tgt + 1]

            # Skip non-overlapping groups
            if e_high_tgt <= e_low_src or e_low_tgt >= e_high_src:
                continue

            # Calculate overlap region
            e_low_overlap = max(e_low_src, e_low_tgt)
            e_high_overlap = min(e_high_src, e_high_tgt)
            lethargy_overlap = np.log(e_high_overlap / e_low_overlap)

            # Distribute flux proportionally to lethargy fraction
            flux_target[idx_tgt] += flux_src * (lethargy_overlap / lethargy_src)

    return flux_target
