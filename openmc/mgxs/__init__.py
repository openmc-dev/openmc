import numpy as np

from openmc.mgxs.groups import EnergyGroups
from openmc.mgxs.library import Library
from openmc.mgxs.mgxs import *
from openmc.mgxs.mdgxs import *

GROUP_STRUCTURES = {}
"""Dictionary of commonly used energy group structures, including "CASMO-X" (where X
is 2, 4, 8, 16, 25, 40 or 70) from the CASMO_ lattice physics code.

.. _CASMO: https://www.studsvik.com/SharepointFiles/CASMO-5%20Development%20and%20Applications.pdf

"""

GROUP_STRUCTURES['CASMO-2'] = np.array([
  0., 6.25e-1, 2.e7])
GROUP_STRUCTURES['CASMO-4'] = np.array([
  0., 6.25e-1, 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-8'] = np.array([
  0., 5.8e-2, 1.4e-1, 2.8e-1, 6.25e-1, 4., 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-16'] = np.array([
  0., 3.e-2, 5.8e-2, 1.4e-1, 2.8e-1, 3.5e-1, 6.25e-1, 8.5e-1,
  9.72e-1, 1.02, 1.097, 1.15, 1.3,  4., 5.53e3, 8.21e5, 2.e7])
GROUP_STRUCTURES['CASMO-25']  = np.array([
    0., 3.e-2, 5.8e-2, 1.4e-1, 2.8e-1, 3.5e-1, 6.25e-1, 9.72e-1, 1.02, 1.097,
    1.15, 1.855, 4., 9.877, 1.5968e1, 1.4873e2, 5.53e3, 9.118e3, 1.11e5, 5.e5,
    8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
GROUP_STRUCTURES['CASMO-40'] = np.array([
    0., 1.5e-2, 3.e-2, 4.2e-2, 5.8e-2, 8.e-2, 1.e-1, 1.4e-1,
    1.8e-1, 2.2e-1, 2.8e-1, 3.5e-1, 6.25e-1, 8.5e-1, 9.5e-1,
    9.72e-1, 1.02, 1.097, 1.15, 1.3, 1.5, 1.855, 2.1, 2.6, 3.3, 4.,
    9.877, 1.5968e1, 2.77e1, 4.8052e1, 1.4873e2, 5.53e3, 9.118e3,
    1.11e5, 5.e5, 8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
GROUP_STRUCTURES['CASMO-70'] = np.array([
    0., 5.e-3, 1.e-2, 1.5e-2, 2.e-2, 2.5e-2, 3.e-2, 3.5e-2, 4.2e-2,
    5.e-2, 5.8e-2, 6.7e-2, 8.e-2, 1.e-1, 1.4e-1, 1.8e-1, 2.2e-1,
    2.5e-1, 2.8e-1, 3.e-1, 3.2e-1, 3.5e-1, 4.e-1, 5.e-1, 6.25e-1,
    7.8e-1, 8.5e-1, 9.1e-1, 9.5e-1, 9.72e-1, 9.96e-1, 1.02, 1.045,
    1.071, 1.097, 1.123, 1.15, 1.3, 1.5, 1.855, 2.1, 2.6, 3.3, 4.,
    9.877, 1.5968e1, 2.77e1, 4.8052e1, 7.5501e1, 1.4873e2,
    3.6726e2, 9.069e2, 1.4251e3, 2.2395e3, 3.5191e3, 5.53e3,
    9.118e3, 1.503e4, 2.478e4, 4.085e4, 6.734e4, 1.11e5, 1.83e5,
    3.025e5, 5.e5, 8.21e5, 1.353e6, 2.231e6, 3.679e6, 6.0655e6, 2.e7])
