from checkvalue import *
from checkvalue import _isinstance

import numpy as np

zs = np.zeros((2,))

print _isinstance(zs[0], Integral)
print _isinstance(zs[0], Real)
print _isinstance(zs[0], (Integral, Real))

print check_iterable_type('thing', zs, (Real, Integral))
