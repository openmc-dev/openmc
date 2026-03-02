import pytest
import openmc


def test_reaction_filter_construction_with_strings():
    f = openmc.ReactionFilter(['(n,elastic)', '(n,gamma)'])
    assert len(f.bins) == 2
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_reaction_filter_construction_with_mt():
    f = openmc.ReactionFilter([2, 102])
    assert len(f.bins) == 2
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_reaction_filter_mixed():
    f = openmc.ReactionFilter([2, '(n,gamma)'])
    assert f.bins[0] == '(n,elastic)'
    assert f.bins[1] == '(n,gamma)'


def test_reaction_filter_single_bin_string():
    f = openmc.ReactionFilter('(n,elastic)')
    assert len(f.bins) == 1
    assert f.bins[0] == '(n,elastic)'


def test_reaction_filter_single_bin_mt():
    f = openmc.ReactionFilter(2)
    assert len(f.bins) == 1
    assert f.bins[0] == '(n,elastic)'


def test_reaction_filter_single_bin_naming():
    f = openmc.ReactionFilter('total')
    assert len(f.bins) == 1
    assert f.bins[0] == '(n,total)'


def test_reaction_filter_invalid_mt():
    with pytest.raises(ValueError, match="No known reaction"):
        openmc.ReactionFilter([999999])


def test_reaction_filter_invalid_string():
    with pytest.raises(ValueError, match="Unknown reaction name"):
        openmc.ReactionFilter(['not-a-reaction'])


def test_reaction_filter_invalid_type():
    with pytest.raises(TypeError, match="Expected str or int"):
        openmc.ReactionFilter([3.14])


def test_reaction_filter_xml_roundtrip():
    f = openmc.ReactionFilter([2, 102], filter_id=42)
    elem = f.to_xml_element()
    f2 = openmc.ReactionFilter.from_xml_element(elem)
    assert f2.id == 42
    assert len(f2.bins) == 2
    assert f2.bins[0] == '(n,elastic)'
    assert f2.bins[1] == '(n,gamma)'


def test_reaction_filter_num_bins():
    f = openmc.ReactionFilter(['(n,elastic)', '(n,fission)', '(n,gamma)'])
    assert f.num_bins == 3


def test_reaction_filter_repr():
    f = openmc.ReactionFilter([2, 102])
    r = repr(f)
    assert 'ReactionFilter' in r


def test_reaction_filter_short_name():
    assert openmc.ReactionFilter.short_name == 'Reaction'


def test_reaction_filter_total_warning():
    """Test that using 'total' emits a warning about ambiguity."""
    with pytest.warns(UserWarning, match="ambiguous"):
        f = openmc.ReactionFilter(['total'])
    assert f.bins[0] == '(n,total)'
