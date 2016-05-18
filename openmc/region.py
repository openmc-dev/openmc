from abc import ABCMeta, abstractmethod
from collections import Iterable

import numpy as np

from openmc.checkvalue import check_type


class Region(object):
    """Region of space that can be assigned to a cell.

    Region is an abstract base class that is inherited by
    :class:`openmc.Halfspace`, :class:`openmc.Intersection`,
    :class:`openmc.Union`, and :class:`openmc.Complement`. Each of those
    respective classes are typically not instantiated directly but rather are
    created through operators of the Surface and Region classes.

    """

    __metaclass__ = ABCMeta

    def __and__(self, other):
        return Intersection(self, other)

    def __or__(self, other):
        return Union(self, other)

    def __invert__(self):
        return Complement(self)

    @abstractmethod
    def __contains__(self, point):
        return False

    @abstractmethod
    def __str__(self):
        return ''

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        elif str(self) != str(other):
            return False
        else:
            return True

    def __ne__(self, other):
        return not self == other

    @staticmethod
    def from_expression(expression, surfaces):
        """Generate a region given an infix expression.

        Parameters
        ----------
        expression : str
            Boolean expression relating surface half-spaces. The possible
            operators are union '|', intersection ' ', and complement '~'. For
            example, '(1 -2) | 3 ~(4 -5)'.
        surfaces : dict
            Dictionary whose keys are suface IDs that appear in the Boolean
            expression and whose values are Surface objects.

        """

        # Strip leading and trailing whitespace
        expression = expression.strip()

        # Convert the string expression into a list of tokens, i.e., operators
        # and surface half-spaces, representing the expression in infix
        # notation.
        i = 0
        i_start = -1
        tokens = []
        while i < len(expression):
            if expression[i] in '()|~ ':
                # If special character appears immediately after a non-operator,
                # create a token with the apporpriate half-space
                if i_start >= 0:
                    j = int(expression[i_start:i])
                    if j < 0:
                        tokens.append(-surfaces[abs(j)])
                    else:
                        tokens.append(+surfaces[abs(j)])

                if expression[i] in '()|~':
                    # For everything other than intersection, add the operator
                    # to the list of tokens
                    tokens.append(expression[i])
                else:
                    # Find next non-space character
                    while expression[i+1] == ' ':
                        i += 1

                    # If previous token is a halfspace or right parenthesis and next token
                    # is not a left parenthese or union operator, that implies that the
                    # whitespace is to be interpreted as an intersection operator
                    if (i_start >= 0 or tokens[-1] == ')') and \
                       expression[i+1] not in ')|':
                        tokens.append(' ')

                i_start = -1
            else:
                # Check for invalid characters
                if expression[i] not in '-+0123456789':
                    raise SyntaxError("Invalid character '{}' in expression"
                                      .format(expression[i]))

                # If we haven't yet reached the start of a word, start one
                if i_start < 0:
                    i_start = i
            i += 1

        # If we've reached the end and we're still in a word, create a
        # half-space token and add it to the list
        if i_start >= 0:
            j = int(expression[i_start:])
            if j < 0:
                tokens.append(-surfaces[abs(j)])
            else:
                tokens.append(+surfaces[abs(j)])

        # The functions below are used to apply an operator to operands on the
        # output queue during the shunting yard algorithm.
        def can_be_combined(region):
            return isinstance(region, Complement) or hasattr(region, 'surface')

        def apply_operator(output, operator):
            r2 = output.pop()
            if operator == ' ':
                r1 = output.pop()
                if isinstance(r1, Intersection) and can_be_combined(r2):
                    r1.nodes.append(r2)
                    output.append(r1)
                elif isinstance(r2, Intersection) and can_be_combined(r1):
                    r2.nodes.insert(0, r1)
                    output.append(r2)
                elif isinstance(r1, Intersection) and isinstance(r2, Intersection):
                    r1.nodes += r2.nodes
                    output.append(r1)
                else:
                    output.append(Intersection(r1, r2))
            elif operator == '|':
                r1 = output.pop()
                if isinstance(r1, Union) and can_be_combined(r2):
                    r1.nodes.append(r2)
                    output.append(r1)
                elif isinstance(r2, Union) and can_be_combined(r1):
                    r2.nodes.insert(0, r1)
                    output.append(r2)
                elif isinstance(r1, Union) and isinstance(r2, Union):
                    r1.nodes += r2.nodes
                    output.append(r1)
                else:
                    output.append(Union(r1, r2))
            elif operator == '~':
                output.append(Complement(r2))

        # The following is an implementation of the shunting yard algorithm to
        # generate an abstract syntax tree for the region expression.
        output = []
        stack = []
        precedence = {'|': 1, ' ': 2, '~': 3}
        associativity = {'|': 'left', ' ': 'left', '~': 'right'}
        for token in tokens:
            if token in (' ', '|', '~'):
                # Normal operators
                while stack:
                    op = stack[-1]
                    if (op not in ('(', ')') and
                        ((associativity[token] == 'right' and
                          precedence[token] < precedence[op]) or
                         (associativity[token] == 'left' and
                          precedence[token] <= precedence[op]))):
                        apply_operator(output, stack.pop())
                    else:
                        break
                stack.append(token)
            elif token == '(':
                # Left parentheses
                stack.append(token)
            elif token == ')':
                # Right parentheses
                while stack[-1] != '(':
                    apply_operator(output, stack.pop())
                    if len(stack) == 0:
                        raise SyntaxError('Mismatched parentheses in '
                                          'region specification.')
                stack.pop()
            else:
                # Surface halfspaces
                output.append(token)
        while stack:
            if stack[-1] in '()':
                raise SyntaxError('Mismatched parentheses in region '
                                  'specification.')
            apply_operator(output, stack.pop())

        # Since we are generating an abstract syntax tree rather than a reverse
        # Polish notation expression, the output queue should have a single item
        # at the end
        return output[0]


class Intersection(Region):
    """Intersection of two or more regions.

    Instances of Intersection are generally created via the __and__ operator
    applied to two instances of :class:`openmc.Region`. This is illustrated in
    the following example:

    >>> equator = openmc.ZPlane(z0=0.0)
    >>> earth = openmc.Sphere(R=637.1e6)
    >>> northern_hemisphere = -earth & +equator
    >>> southern_hemisphere = -earth & -equator
    >>> type(northern_hemisphere)
    <class 'openmc.region.Intersection'>

    Parameters
    ----------
    \*nodes
        Regions to take the intersection of

    Attributes
    ----------
    nodes : list of openmc.Region
        Regions to take the intersection of
    bounding_box : tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box

    """

    def __init__(self, *nodes):
        self.nodes = list(nodes)

    def __iter__(self):
        for n in self.nodes:
            yield n

    def __contains__(self, point):
        """Check whether a point is contained in the region.

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates, :math:`(x',y',z')`, of the point

        Returns
        -------
        bool
            Whether the point is in the region

        """
        return all(point in n for n in self.nodes)

    def __str__(self):
        return '(' + ' '.join(map(str, self.nodes)) + ')'

    @property
    def nodes(self):
        return self._nodes

    @property
    def bounding_box(self):
        lower_left = np.array([-np.inf, -np.inf, -np.inf])
        upper_right = np.array([np.inf, np.inf, np.inf])
        for n in self.nodes:
            lower_left_n, upper_right_n = n.bounding_box
            lower_left[:] = np.maximum(lower_left, lower_left_n)
            upper_right[:] = np.minimum(upper_right, upper_right_n)
        return lower_left, upper_right

    @nodes.setter
    def nodes(self, nodes):
        check_type('nodes', nodes, Iterable, Region)
        self._nodes = nodes


class Union(Region):
    """Union of two or more regions.

    Instances of Union are generally created via the __or__ operator applied to
    two instances of :class:`openmc.Region`. This is illustrated in the
    following example:

    >>> s1 = openmc.ZPlane(z0=0.0)
    >>> s2 = openmc.Sphere(R=637.1e6)
    >>> type(-s2 | +s1)
    <class 'openmc.region.Union'>

    Parameters
    ----------
    \*nodes
        Regions to take the union of

    Attributes
    ----------
    nodes : tuple of openmc.Region
        Regions to take the union of
    bounding_box : tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box

    """

    def __init__(self, *nodes):
        self.nodes = list(nodes)

    def __iter__(self):
        for n in self.nodes:
            yield n

    def __contains__(self, point):
        """Check whether a point is contained in the region.

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates, :math:`(x',y',z')`, of the point

        Returns
        -------
        bool
            Whether the point is in the region

        """
        return any(point in n for n in self.nodes)

    def __str__(self):
        return '(' + ' | '.join(map(str, self.nodes)) + ')'

    @property
    def nodes(self):
        return self._nodes

    @property
    def bounding_box(self):
        lower_left = np.array([np.inf, np.inf, np.inf])
        upper_right = np.array([-np.inf, -np.inf, -np.inf])
        for n in self.nodes:
            lower_left_n, upper_right_n = n.bounding_box
            lower_left[:] = np.minimum(lower_left, lower_left_n)
            upper_right[:] = np.maximum(upper_right, upper_right_n)
        return lower_left, upper_right

    @nodes.setter
    def nodes(self, nodes):
        check_type('nodes', nodes, Iterable, Region)
        self._nodes = nodes


class Complement(Region):
    """Complement of a region.

    The Complement of an existing :class:`openmc.Region` can be created by using
    the __invert__ operator as the following example demonstrates:

    >>> xl = openmc.XPlane(x0=-10.0)
    >>> xr = openmc.XPlane(x0=10.0)
    >>> yl = openmc.YPlane(y0=-10.0)
    >>> yr = openmc.YPlane(y0=10.0)
    >>> inside_box = +xl & -xr & +yl & -yl
    >>> outside_box = ~inside_box
    >>> type(outside_box)
    <class 'openmc.region.Complement'>

    Parameters
    ----------
    node : openmc.Region
        Region to take the complement of

    Attributes
    ----------
    node : openmc.Region
        Regions to take the complement of
    bounding_box : tuple of numpy.array
        Lower-left and upper-right coordinates of an axis-aligned bounding box

    """

    def __init__(self, node):
        self.node = node

    def __contains__(self, point):
        """Check whether a point is contained in the region.

        Parameters
        ----------
        point : 3-tuple of float
            Cartesian coordinates, :math:`(x',y',z')`, of the point

        Returns
        -------
        bool
            Whether the point is in the region

        """
        return point not in self.node

    def __str__(self):
        return '~' + str(self.node)

    @property
    def node(self):
        return self._node

    @node.setter
    def node(self, node):
        check_type('node', node, Region)
        self._node = node

    @property
    def bounding_box(self):
        # Use De Morgan's laws to distribute the complement operator so that it
        # only applies to surface half-spaces, thus allowing us to calculate the
        # bounding box in the usual recursive manner.
        if isinstance(self.node, Union):
            temp_region = Intersection(*[~n for n in self.node.nodes])
        elif isinstance(self.node, Intersection):
            temp_region = Union(*[~n for n in self.node.nodes])
        elif isinstance(self.node, Complement):
            temp_region = self.node.node
        else:
            temp_region = ~self.node
        return temp_region.bounding_box
