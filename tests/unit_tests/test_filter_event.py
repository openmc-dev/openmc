import numpy as np
import pytest
import openmc
from openmc.data.reaction import REACTION_NAME


def test_event_filter_construction_with_strings():
    f = openmc.EventFilter(['(n,elastic)', '(n,gamma)'])
    assert len(f.bins) == 2
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_event_filter_construction_with_mt():
    f = openmc.EventFilter([2, 102])
    assert len(f.bins) == 2
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_event_filter_mixed():
    f = openmc.EventFilter([2, '(n,gamma)'])
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_event_filter_single_bin_string():
    f = openmc.EventFilter('(n,elastic)')
    assert len(f.bins) == 1
    assert f.bins[0] == '(n,elastic)'


def test_event_filter_single_bin_mt():
    f = openmc.EventFilter(2)
    assert len(f.bins) == 1
    assert f.bins[0] == '(n,elastic)'


def test_event_filter_invalid_mt():
    with pytest.raises(ValueError, match="No known reaction"):
        openmc.EventFilter([999999])


def test_event_filter_invalid_string():
    with pytest.raises(ValueError, match="Unknown reaction name"):
        openmc.EventFilter(['not-a-reaction'])


def test_event_filter_invalid_type():
    with pytest.raises(TypeError, match="Expected str or int"):
        openmc.EventFilter([3.14])


def test_event_filter_xml_roundtrip():
    f = openmc.EventFilter([2, 102], filter_id=42)
    elem = f.to_xml_element()
    f2 = openmc.EventFilter.from_xml_element(elem)
    assert f2.id == 42
    assert len(f2.bins) == 2
    assert f2.bins[0] == '(n,elastic)'
    assert f2.bins[1] == '(n,gamma)'


def test_event_filter_num_bins():
    f = openmc.EventFilter(['(n,elastic)', '(n,fission)', '(n,gamma)'])
    assert f.num_bins == 3


def test_event_filter_repr():
    f = openmc.EventFilter([2, 102])
    r = repr(f)
    assert 'EventFilter' in r


def test_event_filter_short_name():
    assert openmc.EventFilter.short_name == 'Event'
