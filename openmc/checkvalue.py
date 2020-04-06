import copy
from collections.abc import Iterable

import numpy as np


def check_type(name, value, expected_type, expected_iter_type=None, *, none_ok=False):
    """Ensure that an object is of an expected type. Optionally, if the object is
    iterable, check that each element is of a particular type.

    Parameters
    ----------
    name : str
        Description of value being checked
    value : object
        Object to check type of
    expected_type : type or Iterable of type
        type to check object against
    expected_iter_type : type or Iterable of type or None, optional
        Expected type of each element in value, assuming it is iterable. If
        None, no check will be performed.
    none_ok : bool, optional
        Whether None is allowed as a value

    """
    if none_ok and value is None:
        return

    if not isinstance(value, expected_type):
        if isinstance(expected_type, Iterable):
            msg = 'Unable to set "{}" to "{}" which is not one of the ' \
                  'following types: "{}"'.format(name, value, ', '.join(
                      [t.__name__ for t in expected_type]))
        else:
            msg = 'Unable to set "{}" to "{}" which is not of type "{}"'.format(
                name, value, expected_type.__name__)
        raise TypeError(msg)

    if expected_iter_type:
        if isinstance(value, np.ndarray):
            if not issubclass(value.dtype.type, expected_iter_type):
                msg = 'Unable to set "{}" to "{}" since each item must be ' \
                      'of type "{}"'.format(name, value,
                                            expected_iter_type.__name__)
                raise TypeError(msg)
            else:
                return

        for item in value:
            if not isinstance(item, expected_iter_type):
                if isinstance(expected_iter_type, Iterable):
                    msg = 'Unable to set "{}" to "{}" since each item must be ' \
                          'one of the following types: "{}"'.format(
                              name, value, ', '.join([t.__name__ for t in
                                                      expected_iter_type]))
                else:
                    msg = 'Unable to set "{}" to "{}" since each item must be ' \
                          'of type "{}"'.format(name, value,
                                                expected_iter_type.__name__)
                raise TypeError(msg)


def check_iterable_type(name, value, expected_type, min_depth=1, max_depth=1):
    """Ensure that an object is an iterable containing an expected type.

    Parameters
    ----------
    name : str
        Description of value being checked
    value : Iterable
        Iterable, possibly of other iterables, that should ultimately contain
        the expected type
    expected_type : type
        type that the iterable should contain
    min_depth : int
        The minimum number of layers of nested iterables there should be before
        reaching the ultimately contained items
    max_depth : int
        The maximum number of layers of nested iterables there should be before
        reaching the ultimately contained items
    """
    # Initialize the tree at the very first item.
    tree = [value]
    index = [0]

    # Traverse the tree.
    while index[0] != len(tree[0]):
        # If we are done with this level of the tree, go to the next branch on
        # the level above this one.
        if index[-1] == len(tree[-1]):
            del index[-1]
            del tree[-1]
            index[-1] += 1
            continue

        # Get a string representation of the current index in case we raise an
        # exception.
        form = '[' + '{:d}, ' * (len(index)-1) + '{:d}]'
        ind_str = form.format(*index)

        # What is the current item we are looking at?
        current_item = tree[-1][index[-1]]

        # If this item is of the expected type, then we've reached the bottom
        # level of this branch.
        if isinstance(current_item, expected_type):
            # Is this deep enough?
            if len(tree) < min_depth:
                msg = 'Error setting "{0}": The item at {1} does not meet the '\
                      'minimum depth of {2}'.format(name, ind_str, min_depth)
                raise TypeError(msg)

            # This item is okay.  Move on to the next item.
            index[-1] += 1

        # If this item is not of the expected type, then it's either an error or
        # on a deeper level of the tree.
        else:
            if isinstance(current_item, Iterable):
                # The tree goes deeper here, let's explore it.
                tree.append(current_item)
                index.append(0)

                # But first, have we exceeded the max depth?
                if len(tree) > max_depth:
                    msg = 'Error setting {0}: Found an iterable at {1}, items '\
                          'in that iterable exceed the maximum depth of {2}' \
                          .format(name, ind_str, max_depth)
                    raise TypeError(msg)

            else:
                # This item is completely unexpected.
                msg = "Error setting {0}: Items must be of type '{1}', but " \
                      "item at {2} is of type '{3}'"\
                      .format(name, expected_type.__name__, ind_str,
                              type(current_item).__name__)
                raise TypeError(msg)


def check_length(name, value, length_min, length_max=None):
    """Ensure that a sized object has length within a given range.

    Parameters
    ----------
    name : str
        Description of value being checked
    value : collections.Sized
        Object to check length of
    length_min : int
        Minimum length of object
    length_max : int or None, optional
        Maximum length of object. If None, it is assumed object must be of
        length length_min.

    """

    if length_max is None:
        if len(value) < length_min:
            msg = 'Unable to set "{}" to "{}" since it must be at least of ' \
                  'length "{}"'.format(name, value, length_min)
            raise ValueError(msg)
    elif not length_min <= len(value) <= length_max:
        if length_min == length_max:
            msg = 'Unable to set "{}" to "{}" since it must be of ' \
                  'length "{}"'.format(name, value, length_min)
        else:
            msg = 'Unable to set "{}" to "{}" since it must have length ' \
                  'between "{}" and "{}"'.format(name, value, length_min,
                                                 length_max)
        raise ValueError(msg)


def check_value(name, value, accepted_values):
    """Ensure that an object's value is contained in a set of acceptable values.

    Parameters
    ----------
    name : str
        Description of value being checked
    value : collections.Iterable
        Object to check
    accepted_values : collections.Container
        Container of acceptable values

    """

    if value not in accepted_values:
        msg = 'Unable to set "{0}" to "{1}" since it is not in "{2}"'.format(
            name, value, accepted_values)
        raise ValueError(msg)


def check_less_than(name, value, maximum, equality=False):
    """Ensure that an object's value is less than a given value.

    Parameters
    ----------
    name : str
        Description of the value being checked
    value : object
        Object to check
    maximum : object
        Maximum value to check against
    equality : bool, optional
        Whether equality is allowed. Defaults to False.

    """

    if equality:
        if value > maximum:
            msg = 'Unable to set "{0}" to "{1}" since it is greater than ' \
                  '"{2}"'.format(name, value, maximum)
            raise ValueError(msg)
    else:
        if value >= maximum:
            msg = 'Unable to set "{0}" to "{1}" since it is greater than ' \
                  'or equal to "{2}"'.format(name, value, maximum)
            raise ValueError(msg)


def check_greater_than(name, value, minimum, equality=False):
    """Ensure that an object's value is greater than a given value.

    Parameters
    ----------
    name : str
        Description of the value being checked
    value : object
        Object to check
    minimum : object
        Minimum value to check against
    equality : bool, optional
        Whether equality is allowed. Defaults to False.

    """

    if equality:
        if value < minimum:
            msg = 'Unable to set "{0}" to "{1}" since it is less than ' \
                  '"{2}"'.format(name, value, minimum)
            raise ValueError(msg)
    else:
        if value <= minimum:
            msg = 'Unable to set "{0}" to "{1}" since it is less than ' \
                  'or equal to "{2}"'.format(name, value, minimum)
            raise ValueError(msg)


def check_filetype_version(obj, expected_type, expected_version):
    """Check filetype and version of an HDF5 file.

    Parameters
    ----------
    obj : h5py.File
        HDF5 file to check
    expected_type : str
        Expected file type, e.g. 'statepoint'
    expected_version : int
        Expected major version number.

    """
    try:
        this_filetype = obj.attrs['filetype'].decode()
        this_version = obj.attrs['version']

        # Check filetype
        if this_filetype != expected_type:
            raise IOError('{} is not a {} file.'.format(
                obj.filename, expected_type))

        # Check version
        if this_version[0] != expected_version:
            raise IOError('{} file has a version of {} which is not '
                          'consistent with the version expected by OpenMC, {}'
                          .format(this_filetype,
                                  '.'.join(str(v) for v in this_version),
                                  expected_version))
    except AttributeError:
        raise IOError('Could not read {} file. This most likely means the '
                      'file was produced by a different version of OpenMC than '
                      'the one you are using.'.format(obj.filename))


class CheckedList(list):
    """A list for which each element is type-checked as it's added

    Parameters
    ----------
    expected_type : type or Iterable of type
        Type(s) which each element should be
    name : str
        Name of data being checked
    items : Iterable, optional
        Items to initialize the list with

    """

    def __init__(self, expected_type, name, items=None):
        super().__init__()
        self.expected_type = expected_type
        self.name = name
        if items is not None:
            for item in items:
                self.append(item)

    def __add__(self, other):
        new_instance = copy.copy(self)
        new_instance += other
        return new_instance

    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        check_type('CheckedList add operand', other, Iterable,
                   self.expected_type)
        for item in other:
            self.append(item)
        return self

    def append(self, item):
        """Append item to list

        Parameters
        ----------
        item : object
            Item to append

        """
        check_type(self.name, item, self.expected_type)
        super().append(item)

    def insert(self, index, item):
        """Insert item before index

        Parameters
        ----------
        index : int
            Index in list
        item : object
            Item to insert

        """
        check_type(self.name, item, self.expected_type)
        super().insert(index, item)
