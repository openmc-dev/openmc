from abc import ABCMeta, abstractmethod
from collections import Iterable

from openmc.checkvalue import check_type


class Region(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __str__(self):
        return ''


class Intersection(Region):
    """Intersection of two or more regions.

    Parameters
    ----------
    *nodes
        Regions to take the intersection of

    Attributes
    ----------
    nodes : tuple of Region
        Regions to take the intersection of

    """

    def __init__(self, *nodes):
        self.nodes = nodes

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, nodes):
        check_type('nodes', nodes, Iterable, Region)
        self._nodes = nodes

    def __str__(self):
        return '(' + ' '.join(map(str, self.nodes)) + ')'


class Union(Region):
    """Union of two or more regions.

    Parameters
    ----------
    *nodes
        Regions to take the union of

    Attributes
    ----------
    nodes : tuple of Region
        Regions to take the union of

    """

    def __init__(self, *nodes):
        self.nodes = nodes

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, nodes):
        check_type('nodes', nodes, Iterable, Region)
        self._nodes = nodes

    def __str__(self):
        return '(' + ' ^ '.join(map(str, self.nodes)) + ')'


class Complement(Region):
    """Complement of a region.

    Parameters
    ----------
    node : Region
        Region to take the complement of

    Attributes
    ----------
    node : Region
        Regions to take the complement of

    """

    def __init__(self, node):
        self.node = node

    @property
    def node(self):
        return self._node

    @node.setter
    def node(self, node):
        check_type('node', node, Region)
        self._node = node

    def __str__(self):
        return '~' + str(self.node)
