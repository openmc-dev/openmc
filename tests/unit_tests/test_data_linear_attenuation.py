import numpy as np

import openmc.data.photon_attenuation as linear_attenuation
import pytest
from openmc.data.photon_attenuation import linear_attenuation_xs

import openmc.data

PHOTON_MTS = (502, 504, 515, 517, 522)


@pytest.mark.parametrize("symbol", ["Cu", "Pu"])
def test_linear_attenuation_xs_matches_sum(elements_endf, symbol, monkeypatch):
    """linear_attenuation_xs should reproduce the sum of the relevant
    reaction channels from IncidentPhoton.reactions.
    """
    element = elements_endf[symbol]
    assert isinstance(element, openmc.data.IncidentPhoton)

    # Stub out the data lookup so we don't depend on a DataLibrary/cross_sections.xml
    monkeypatch.setattr(linear_attenuation, "_get_photon_data", lambda name: element)

    # Call the helper
    xs_sum = linear_attenuation_xs(symbol, temperature=293.6)

    # If the element has no relevant reactions, helper should return None
    has_relevant = any(mt in element.reactions for mt in PHOTON_MTS)
    if not has_relevant:
        assert xs_sum is None
        return

    assert isinstance(xs_sum, openmc.data.Sum)

    # Compare against explicit sum of reaction cross sections
    energy = np.logspace(2, 4, 50)  
    expected = np.zeros_like(energy)
    for mt in PHOTON_MTS:
        if mt in element.reactions:
            expected += element.reactions[mt].xs(energy)

    actual = xs_sum(energy)
    assert np.allclose(actual, expected)


def test_linear_attenuation_xs_returns_none_when_no_photon_data(monkeypatch):
    """If _get_photon_data returns None, the helper should return None."""
    # Force _get_photon_data to return None regardless of nuclide
    monkeypatch.setattr(linear_attenuation, "_get_photon_data", lambda name: None)

    xs_sum = linear_attenuation_xs("NonExistent", temperature=300.0)
    assert xs_sum is None


# def test_linear_attenuation_xs_temperature_fallback(monkeypatch):
#     """When exact temperature is not present, the closest available
#     temperature should be selected from the xs dict.
#     """
#
#     class DummyXS:
#         def __init__(self, value: float):
#             self._value = value
#
#         def __call__(self, E):
#             E = np.asanyarray(E)
#             return np.full_like(E, self._value, dtype=float)
#
#     class DummyReaction:
#         def __init__(self, mt: int, xs):
#             self.mt = mt
#             self.xs = xs
#
#     class DummyPhotonData:
#         def __init__(self):
#             # xs for two temperatures, keyed as "<T>K"
#             self.reactions = {
#                 502: DummyReaction(502, {"290K": DummyXS(1.0), "600K": DummyXS(2.0)}),
#                 504: DummyReaction(504, {"290K": DummyXS(10.0), "600K": DummyXS(20.0)}),
#             }
#
#     dummy_data = DummyPhotonData()
#
#     # Use dummy photon data instead of reading from files/DataLibrary
#     monkeypatch.setattr(photon_xs, "_get_photon_data", lambda name: dummy_data)
#
#     energy = np.array([1.0, 2.0, 5.0])
#
#     # 295 K is closer to 290 K -> expect use of 290K datasets
#     xs_295 = linear_attenuation_xs("dummy", temperature=295.0)
#     assert isinstance(xs_295, photon_xs.Sum)
#     vals_295 = xs_295(energy)
#     # 502: 1.0, 504: 10.0  -> total 11.0
#     assert np.allclose(vals_295, 11.0)
#
#     # 500 K is closer to 600 K -> expect use of 600K datasets
#     xs_500 = linear_attenuation_xs("dummy", temperature=500.0)
#     assert isinstance(xs_500, photon_xs.Sum)
#     vals_500 = xs_500(energy)
#     # 502: 2.0, 504: 20.0  -> total 22.0
#     assert np.allclose(vals_500, 22.0)
