import numpy as np


class EqualityMixin(object):
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

    def __ne__(self, other):
        return not self.__eq__(other)
