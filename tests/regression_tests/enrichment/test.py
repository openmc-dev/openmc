import os
import sys

import numpy as np

from openmc import Material
from openmc.data import NATURAL_ABUNDANCE, atomic_mass


def test_enrichment():
    # This test doesn't require an OpenMC run.  We just need to make sure the
    # element.expand() method expands Uranium to the proper enrichment.

    uranium = Material()
    uranium.add_element('U', 1.0, 'wo', 4.95)
    densities = uranium.get_nuclide_densities()

    sum_densities = 0.
    for nuc in densities.keys():
        assert nuc in ('U234', 'U235', 'U236', 'U238')
        sum_densities += densities[nuc][1]

    # Compute the weight percent U235
    enrichment = densities['U235'][1] / sum_densities
    assert np.isclose(enrichment, 0.0495, rtol=1.e-8)

    # Compute the ratio of U234/U235
    u234_to_u235 = densities['U234'][1] / densities['U235'][1]
    assert np.isclose(u234_to_u235, 0.0089, rtol=1.e-8)

    # Compute the ratio of U236/U235
    u236_to_u235 = densities['U236'][1] / densities['U235'][1]
    assert np.isclose(u236_to_u235, 0.0046, rtol=1.e-8)
