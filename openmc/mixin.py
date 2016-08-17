import numpy as np


class Equality(object):
    """A Class which provides generic __eq__ and __ne__ functionality which
    can easily be inherited by downstream classes.
    """

    def __eq__(self, other):
        eqval = True
        if isinstance(other, type(self)):
            for key, value in self.__dict__.items():
                if key in other.__dict__:
                    if not np.array_equal(value, other.__dict__.get(key)):
                        eqval = False
                else:
                    eqval = False
        else:
            eqval = False

        return eqval

    def __ne__(self, other):
        return not self.__eq__(other)
