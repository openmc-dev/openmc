from dataclasses import replace
import numpy as np
import pytest
from sandbox.tokomak_analysis.config import Config
from sandbox.tokomak_analysis.materials import make_materials


def test_pebble_volume_density_and_enrichment():
    c=Config();m=make_materials(c)
    bed=m['breeder0'];solid_density=c.li2o_solid_density
    he_density=c.purge_pressure_mpa*1e6*.004002602/(8.314462618*c.breeder_temperature)/1000
    assert bed.get_mass_density()==pytest.approx(c.breeder_packing*solid_density+(1-c.breeder_packing)*he_density)
    atoms=bed.get_nuclide_atom_densities()
    assert atoms['Li6']/(atoms['Li6']+atoms['Li7'])==pytest.approx(c.enrichment[0])
    assert (atoms['Li6']+atoms['Li7'])/sum(v for k,v in atoms.items() if k.startswith('O'))==pytest.approx(2)
    assert all(x.get_mass_density()>0 for x in m.values())


def test_alternative_multiplier():
    m=make_materials(replace(Config(),multiplier='Be12Ti'))
    names=m['multiplier'].get_nuclides()
    assert 'Be9' in names and any(n.startswith('Ti') for n in names)
