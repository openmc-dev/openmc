from pathlib import Path
from unittest import mock

import numpy as np
import pytest

import openmc
import openmc.deplete
import openmc.lib
from openmc.deplete import Chain


def test_materials_deplete():
    pristine_material_1 = openmc.Material()
    pristine_material_1.add_nuclide("Ni58", 1.)
    pristine_material_1.set_density("g/cm3", 7.87)
    pristine_material_1.depletable = True
    pristine_material_1.temperature = 293.6
    pristine_material_1.volume = 1.

    pristine_material_2 = openmc.Material()
    pristine_material_2.add_nuclide("Ni60", 1.)
    pristine_material_2.set_density("g/cm3", 7.87)
    pristine_material_2.depletable = True
    pristine_material_2.temperature = 293.6
    pristine_material_2.volume = 1.

    pristine_materials = openmc.Materials([pristine_material_1, pristine_material_2])

    mg_flux = [0.5e11] * 42

    chain = Chain.from_xml(
        Path(__file__).parents[1] / "chain_ni.xml"
    )

    depleted_material = pristine_materials.deplete(
        multigroup_fluxes=[mg_flux, mg_flux],
        energy_group_structures=["VITAMIN-J-42", "VITAMIN-J-42"],
        timesteps=[100, 100],
        source_rates=[1e19, 0.0],
        timestep_units="d",
        chain_file=chain,
    )

    assert list(depleted_material.keys()) == [pristine_material_1.id, pristine_material_2.id]
    for mat_id, materials in depleted_material.items():
        for i_step, material in enumerate(materials):
            assert isinstance(material, openmc.Material)
            if i_step > 0:
                assert len(material.get_nuclides()) > 1
            assert mat_id == material.id

    mats = depleted_material[pristine_material_1.id]
    Co58_mat_1_step_0 = mats[0].get_nuclide_atom_densities("Co58").get("Co58", 0.0)
    Co58_mat_1_step_1 = mats[1].get_nuclide_atom_densities("Co58")["Co58"]
    Co58_mat_1_step_2 = mats[2].get_nuclide_atom_densities("Co58")["Co58"]

    assert Co58_mat_1_step_0 == 0.0
    # Co58 is the main activation product of Ni58 in the first irradiation step.
    # It then decays in the second cooling step (flux = 0)
    assert Co58_mat_1_step_1 > 0.0 and Co58_mat_1_step_1 > Co58_mat_1_step_2

    Ni59_mat_1_step_0 = mats[0].get_nuclide_atom_densities("Ni59").get("Ni59", 0.0)
    Ni59_mat_1_step_1 = mats[1].get_nuclide_atom_densities("Ni59")["Ni59"]
    Ni59_mat_1_step_2 = mats[2].get_nuclide_atom_densities("Ni59")["Ni59"]

    assert Ni59_mat_1_step_0 == 0.0
    # Ni59 is one of the main activation product of Ni60 in the first irradiation
    # step. It then decays in the second cooling step (flux = 0)
    assert Ni59_mat_1_step_1 > 0.0 and Ni59_mat_1_step_1 > Ni59_mat_1_step_2


def test_export_duplicate_materials_to_xml(run_in_tmpdir):
    """
    Test exporting Materials to xml with a duplicate and checking that only
    unique entities are exported.
    """
    my_mat = openmc.Material(name="my_mat")
    my_mat2 = openmc.Material(name="my_mat2")

    materials = openmc.Materials([my_mat, my_mat2, my_mat])

    materials.export_to_xml("materials.xml")

    materials_in = openmc.Materials.from_xml("materials.xml")
    assert len(materials_in) == 2


def test_materials_deplete_length_mismatch():
    mats = openmc.Materials([openmc.Material()])

    with pytest.raises(ValueError, match="multigroup_fluxes length"):
        mats.deplete(
            multigroup_fluxes=[],
            energy_group_structures=["VITAMIN-J-42"],
            timesteps=[1.0],
            source_rates=1.0,
        )

    with pytest.raises(ValueError, match="energy_group_structures length"):
        mats.deplete(
            multigroup_fluxes=[[1.0]],
            energy_group_structures=[],
            timesteps=[1.0],
            source_rates=1.0,
        )


def test_materials_deplete_missing_volume(monkeypatch):
    mat = openmc.Material()
    mat.add_nuclide("Ni58", 1.0)
    mat.set_density("g/cm3", 7.87)

    mats = openmc.Materials([mat])

    class DummySession:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(openmc.lib, "TemporarySession", DummySession)

    chain = Path(__file__).parents[1] / "chain_ni.xml"
    with pytest.raises(ValueError, match="has no volume"):
        mats.deplete(
            multigroup_fluxes=[[1.0]],
            energy_group_structures=["VITAMIN-J-42"],
            timesteps=[1.0],
            source_rates=1.0,
            chain_file=chain,
        )


def test_materials_deplete_groups_by_energy_temperature(monkeypatch):
    """Materials.deplete groups materials by (energy, temperature): the group
    cross section table is built once per distinct group, the micros keep
    material order, and each equals the ungrouped from_multigroup_flux result.
    """
    import openmc.deplete.microxs as microxs_mod

    chain = Path(__file__).parents[1] / "chain_ni.xml"
    energy = "VITAMIN-J-42"
    n_groups = 42

    rng = np.random.default_rng(0)
    fluxes = [rng.random(n_groups) for _ in range(3)]
    temps = [293.6, 293.6, 600.0]

    def make_mat(nuclide, temperature):
        m = openmc.Material()
        m.add_nuclide(nuclide, 1.0)
        m.set_density("g/cm3", 7.87)
        m.depletable = True
        m.temperature = temperature
        m.volume = 1.0
        return m

    # mat0 and mat1 share (energy, T=293.6); mat2 differs only in temperature
    mats = openmc.Materials([
        make_mat("Ni58", temps[0]),
        make_mat("Ni60", temps[1]),
        make_mat("Ni62", temps[2]),
    ])

    # Spy on the table build and capture the micros the operator receives
    build_spy = mock.MagicMock(wraps=microxs_mod._build_xs_table_ce)
    monkeypatch.setattr(microxs_mod, "_build_xs_table_ce", build_spy)

    captured = {}
    def fake_operator(*args, **kwargs):
        captured["micros"] = kwargs["micros"]
        raise StopIteration
    monkeypatch.setattr(openmc.deplete, "IndependentOperator", fake_operator)

    with pytest.raises(StopIteration):
        mats.deplete(
            multigroup_fluxes=fluxes,
            energy_group_structures=[energy, energy, energy],
            timesteps=[1.0],
            source_rates=1.0,
            chain_file=chain,
        )

    # Two distinct (energy, temperature) groups -> two builds, not three
    assert build_spy.call_count == 2

    micros = captured["micros"]
    assert len(micros) == 3
    assert all(m is not None for m in micros)

    # Each per-material micro equals the ungrouped single-flux result
    with openmc.lib.TemporarySession():
        for flux, temperature, micro in zip(fluxes, temps, micros):
            ref = openmc.deplete.MicroXS.from_multigroup_flux(
                energies=energy, multigroup_flux=flux, chain_file=chain,
                temperature=temperature)
            assert micro.data == pytest.approx(ref.data)
