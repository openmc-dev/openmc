"""D1S ownership and snapshot lifetime, using the repository's evaluated chain."""
from pathlib import Path

import pytest

import openmc
import openmc.lib
from openmc import exceptions


CHAIN = Path(__file__).parents[1] / 'chain_ni.xml'


def test_d1s_discovery_and_snapshot_lifetime(tmp_path, monkeypatch, mpi_intracomm):
    monkeypatch.setenv('OPENMC_CHAIN_FILE', str(CHAIN))
    material = openmc.Material()
    material.add_nuclide('Ni60', 1.0)
    material.set_density('g/cm3', 1.0)
    sphere = openmc.Sphere(r=1, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-sphere)
    model = openmc.Model(geometry=openmc.Geometry([cell]))
    model.settings.run_mode = 'fixed source'
    model.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    model.settings.particles = 10
    model.settings.batches = 2
    model.settings.use_decay_photons = True
    model.settings.photon_transport = True
    model.settings.output = {'path': str(tmp_path), 'summary': False}
    try:
        model.init_lib(directory=str(tmp_path), output=False, intracomm=mpi_intracomm)
        products = openmc.lib.get_reaction_products('Ni60', include_data=True)
        decay = [(p, d) for p in products for d in p['distributions']
                 if d['angle_energy_type'] == 'decay-photon']
        assert decay
        for p, d in decay:
            assert p['particle'] == 'photon'
            assert p['emission_mode'] == 'delayed'
            assert p['parent_nuclide_index'] >= 0
            assert d['read_supported']
            assert not d['mutation_supported']
            assert not d['energy_mutation_supported']
            assert 'borrowed' in d['mutation_reason']
            assert d['data']['ownership'] == 'borrowed-from-chain-nuclide'
            assert d['data']['spectrum']['type'] in ('discrete', 'mixture', 'tabular')
        # Inspection also works after simulation initialization and finalization.
        openmc.lib.simulation_init()
        during = openmc.lib.get_reaction_products('Ni60')
        assert len(during) == len(products)
        openmc.lib.simulation_finalize()
        assert len(openmc.lib.get_reaction_products('Ni60')) == len(products)
        p, d = decay[0]
        snapshot = openmc.lib.get_emission_distribution(
            'Ni60', p['mt'], p['product_index'], d['distribution_index'])
        index = openmc.lib.nuclides['Ni60']._index
    finally:
        model.finalize_lib()
    # Copies remain usable after all native transport/chain objects are gone.
    assert snapshot['type'] == 'decay-photon'
    assert snapshot['spectrum']['integral'] > 0
    with pytest.raises(exceptions.OutOfBoundsError):
        openmc.lib.get_reaction_products(index)
