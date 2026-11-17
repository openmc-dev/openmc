"""Reaction metadata round trips and backward compatibility."""

import h5py
import pytest

import openmc.data


@pytest.fixture
def data():
    data = openmc.data.IncidentNeutron('Am242_m1', 95, 242, 1, 240.0, [])
    data.reactions[51] = openmc.data.Reaction(51)
    data[51].q_reaction = 48600.0
    return data


@pytest.mark.parametrize('values', [(None, None, None), (0.0, 0.0, 0),
                                   (48600.0, -1.2e6, 33)])
def test_metadata_hdf5(tmp_path, data, values):
    data.excitation_energy, data[51].q_mass_difference, data[51].breakup_flag = values
    path = tmp_path / 'data.h5'
    data.export_to_hdf5(path)
    restored = openmc.data.IncidentNeutron.from_hdf5(path)
    assert (restored.excitation_energy, restored[51].q_mass_difference,
            restored[51].breakup_flag) == values
    assert restored[51].q_reaction == 48600.0
    with h5py.File(path) as handle:
        assert 'Q_value' not in handle['Am242_m1/reactions/reaction_051'].attrs


def test_legacy_q_value(tmp_path, data):
    path = tmp_path / 'data.h5'
    data.export_to_hdf5(path)
    with h5py.File(path, 'r+') as handle:
        attrs = handle['Am242_m1/reactions/reaction_051'].attrs
        attrs['Q_value'] = -1.0e6

    # Prefer q_reaction when both names are present.
    assert openmc.data.IncidentNeutron.from_hdf5(path)[51].q_reaction == 48600.0
    with h5py.File(path, 'r+') as handle:
        del handle['Am242_m1/reactions/reaction_051'].attrs['q_reaction']
        handle.attrs['version'] = [3, 0]

    restored = openmc.data.IncidentNeutron.from_hdf5(path)
    assert restored[51].q_reaction == -1.0e6
    assert restored[51].q_mass_difference is None
    assert restored[51].breakup_flag is None
    assert restored.excitation_energy is None


def test_q_value_deprecated():
    reaction = openmc.data.Reaction(51)
    with pytest.warns(FutureWarning, match='use q_reaction'):
        reaction.q_value = 48600.0
    assert reaction.q_reaction == 48600.0
    reaction.q_reaction = -1.0e6
    with pytest.warns(FutureWarning, match='use q_reaction'):
        assert reaction.q_value == -1.0e6
