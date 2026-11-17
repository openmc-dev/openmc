"""The configured QoI is independent of transport/perturbation constituents."""
from types import SimpleNamespace
import pytest
from sandbox.analysis_codes import production


def configure(monkeypatch):
    monkeypatch.setattr(production,'settings',lambda:SimpleNamespace(
        PRODUCTION_CHANNELS={'Li6':[105],'Li7':list(range(52,83))}))


def test_lithium_qoi_with_lead_without_exclusion_setting(monkeypatch):
    configure(monkeypatch)
    selected=production.response_mapping({'Li6','Li7','Pb204','Pb206','Pb207','Pb208'})
    assert selected=={'Li6':[105],'Li7':list(range(52,83))}


def test_absent_configured_constituents_are_not_scored(monkeypatch):
    configure(monkeypatch)
    assert production.response_mapping({'Li6','Pb208'})=={'Li6':[105]}


def test_unknown_requested_production_channel_still_needs_verification(monkeypatch):
    monkeypatch.setattr(production,'settings',lambda:SimpleNamespace(PRODUCTION_CHANNELS={'Pb208':[105]}))
    with pytest.raises(ValueError,match='Unverified production channels'):
        production.response_mapping({'Pb208'})


def test_empty_qoi_rejected(monkeypatch):
    configure(monkeypatch)
    with pytest.raises(ValueError,match='No configured production'):
        production.response_mapping({'Pb208'})
