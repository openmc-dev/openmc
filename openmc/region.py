from abc import ABCMeta, abstractmethod
from collections import Iterable

from openmc.checkvalue import check_type


class Region(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __str__(self):
        return ''

    @classmethod
    def from_expression(cls, expression, surfaces):
        """Generate a region given an infix expression.

        Parameters
        ----------
        expression : str
            Boolean expression relating surface half-spaces. The possible
            operators are union '^', intersection ' ', and complement '~'. For
            example, '(1 -2) ^ 3 ~(4 -5)'.
        surfaces : dict
            Dictionary whose keys are suface IDs that appear in the Boolean
            expression and whose values are Surface objects.

        """

        # Convert the string expression into a list of tokens, i.e., operators
        # and surface half-spaces, representing the expression in infix
        # notation.
        i = 0
        i_start = -1
        tokens = []
        while i < len(expression):
            if expression[i] in '()^~ ':
                # If special character appears immediately after a non-operator,
                # create a token with the apporpriate half-space
                if i_start >= 0:
                    j = int(expression[i_start:i])
                    if j < 0:
                        tokens.append(surfaces[abs(j)].negative)
                    else:
                        tokens.append(surfaces[abs(j)].positive)

                if expression[i] in '()^~':
                    # For everything other than intersection, add the operator
                    # to the list of tokens
                    tokens.append(expression[i])
                else:
                    # For spaces, we need to check the context further. If it
                    # doesn't appear before a right parentheses or union, it is
                    # interpreted to be as an intersection operator
                    while expression[i+1] == ' ':
                        i += 1

                    if i_start >= 0 and expression[i+1] not in ')^':
                        tokens.append(' ')

                i_start = -1
            else:
                # Check for invalid characters
                if expression[i] not in '-0123456789':
                    raise SyntaxError('Invalid character in expression')

                # If we haven't yet reached the start of a word, start one
                if i_start < 0:
                    i_start = i
            i += 1

        # If we've reached the end and we're still in a word, create a
        # half-space token and add it to the list
        if i_start >= 0:
            j = int(expression[i_start:])
            if j < 0:
                tokens.append(surfaces[abs(j)].negative)
            else:
                tokens.append(surfaces[abs(j)].positive)

        # This function is used below to apply an operator to operands on the
        # output queue during the shunting yard algorithm.
        def apply_operator(output, operator):
            r2 = output.pop()
            if operator == ' ':
                r1 = output.pop()
                output.append(Intersection(r1, r2))
            elif operator == '^':
                r1 = output.pop()
                output.append(Union(r1, r2))
            elif operator == '~':
                output.append(Complement(r2))


        # The following is an implementation of the shunting yard algorithm to
        # generate an abstract syntax tree for the region expression.
        output = []
        stack = []
        precedence = {'^': 1, ' ': 2, '~': 3}
        associativity = {'^': 'left', ' ': 'left', '~': 'right'}
        for token in tokens:
            if token in (' ', '^', '~'):
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
        self.nodes = list(nodes)

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
        self.nodes = list(nodes)

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
