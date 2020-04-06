from numbers import Integral
from warnings import warn

import numpy as np

import openmc.checkvalue as cv


class EqualityMixin:
    """A Class which provides generic __eq__ and __ne__ functionality which
    can easily be inherited by downstream classes.
    """

    def __eq__(self, other):
        if isinstance(other, type(self)):
            for key, value in self.__dict__.items():
                if not np.array_equal(value, other.__dict__.get(key)):
                    return False
        else:
            return False

        return True


class IDWarning(UserWarning):
    pass


class IDManagerMixin:
    """A Class which automatically manages unique IDs.

    This mixin gives any subclass the ability to assign unique IDs through an
    'id' property and keeps track of which ones have already been
    assigned. Crucially, each subclass must define class variables 'next_id' and
    'used_ids' as they are used in the 'id' property that is supplied here.

    """

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, uid):
        # The first time this is called for a class, we search through the MRO
        # to determine which class actually holds next_id and used_ids. Since
        # next_id is an integer (immutable), we can't modify it directly through
        # the instance without just creating a new attribute
        try:
            cls = self._id_class
        except AttributeError:
            for cls in self.__class__.__mro__:
                if 'next_id' in cls.__dict__:
                    break

        if uid is None:
            while cls.next_id in cls.used_ids:
                cls.next_id += 1
            self._id = cls.next_id
            cls.used_ids.add(cls.next_id)
        else:
            name = cls.__name__
            cv.check_type('{} ID'.format(name), uid, Integral)
            cv.check_greater_than('{} ID'.format(name), uid, 0, equality=True)
            if uid in cls.used_ids:
                msg = 'Another {} instance already exists with id={}.'.format(
                    name, uid)
                warn(msg, IDWarning)
            else:
                cls.used_ids.add(uid)
            self._id = uid


def reset_auto_ids():
    """Reset counters for all auto-generated IDs"""
    for cls in IDManagerMixin.__subclasses__():
        cls.used_ids.clear()
        cls.next_id = 1


def reserve_ids(ids, cls=None):
    """Reserve a set of IDs that won't be used for auto-generated IDs.

    Parameters
    ----------
    ids : iterable of int
        IDs to reserve
    cls : type or None
        Class for which IDs should be reserved (e.g., :class:`openmc.Cell`). If
        None, all classes that have auto-generated IDs will be used.

    """
    if cls is None:
        for cls in IDManagerMixin.__subclasses__():
            cls.used_ids |= set(ids)
    else:
        cls.used_ids |= set(ids)


def set_auto_id(next_id):
    """Set the next ID for auto-generated IDs.

    Parameters
    ----------
    next_id : int
        The next ID to assign to objects with auto-generated IDs.

    """
    for cls in IDManagerMixin.__subclasses__():
        cls.next_id = next_id
