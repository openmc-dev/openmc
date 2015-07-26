def check_type(name, value, expected_type, expected_iter_type=None):
    """Ensure that an object is of an expected type. Optionally, if the object is
    iterable, check that each element is of a particular type.

    Parameters
    ----------
    name : str
        Description of value being checked
    value : object
        Object to check type of
    expected_type : type
        type to check object against
    expected_iter_type : type or None, optional
        Expected type of each element in value, assuming it is iterable. If
        None, no check will be performed.

    """

    if not isinstance(value, expected_type):
        msg = 'Unable to set {0} to {1} which is not of type {2}'.format(
            name, value, expected_type.__name__)
        raise ValueError(msg)

    if expected_iter_type:
        for item in value:
            if not isinstance(item, expected_iter_type):
                msg = 'Unable to set {0} to {1} since each item must be ' \
                      'of type {2}'.format(name, value,
                                           expected_iter_type.__name__)
                raise ValueError(msg)


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
        if len(value) != length_min:
            msg = 'Unable to set {0} to {1} since it must be of ' \
                  'length {2}'.format(name, value, length_min)
            raise ValueError(msg)
    elif not length_min <= len(value) <= length_max:
        if length_min == length_max:
            msg = 'Unable to set {0} to {1} since it must be of ' \
                  'length {2}'.format(name, value, length_min)
        else:
            msg = 'Unable to set {0} to {1} since it must have length ' \
                  'between {2} and {3}'.format(name, value, length_min,
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
        msg = 'Unable to set {0} to {1} since it is not in {2}'.format(
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
        Whether equality is allowed. Defaluts to False.

    """

    if equality:
        if value > maximum:
            msg = 'Unable to set {0} to {1} since it is greater than ' \
                  '{2}'.format(name, value, maximum)
            raise ValueError(msg)
    else:
        if value >= maximum:
            msg = 'Unable to set {0} to {1} since it is greater than ' \
                  'or equal to {2}'.format(name, value, maximum)
            raise ValueError(msg)

def check_greater_than(name, value, minimum, equality=False):
    """Ensure that an object's value is less than a given value.

    Parameters
    ----------
    name : str
        Description of the value being checked
    value : object
        Object to check
    minimum : object
        Minimum value to check against
    equality : bool, optional
        Whether equality is allowed. Defaluts to False.

    """

    if equality:
        if value < minimum:
            msg = 'Unable to set {0} to {1} since it is less than ' \
                  '{2}'.format(name, value, minimum)
            raise ValueError(msg)
    else:
        if value <= minimum:
            msg = 'Unable to set {0} to {1} since it is less than ' \
                  'or equal to {2}'.format(name, value, minimum)
            raise ValueError(msg)
