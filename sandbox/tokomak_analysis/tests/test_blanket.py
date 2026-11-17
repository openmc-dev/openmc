from dataclasses import replace
import numpy as np
from sandbox.tokomak_analysis import Config,build_model
from sandbox.tokomak_analysis.geometry.primitives import polygon_volume


def test_parent_module_instances_and_layers(assembly):
    a=assembly;a.model.geometry.determine_paths()
    modules=[i for i,r in a.registry.records.items() if r['component']=='blanket_module']
    assert len(modules)==a.config.poloidal_modules
    assert sum(len(a.registry.cells[i].paths) for i in modules)==a.config.tf_coils*a.config.toroidal_modules_per_sector*a.config.poloidal_modules
    breeder=[r for r in a.registry.records.values() if r.get('kind')=='breeder']
    assert len(breeder)==a.config.poloidal_modules*a.config.breeder_cycles
    for x,y in zip(a.layers,a.layers[1:]):np.testing.assert_allclose(x['outer'].points,y['inner'].points,atol=1e-10)


def test_fidelity_envelopes_match():
    geometries=[build_model(replace(Config(),fidelity=i),tallies=False) for i in (1,2,3)]
    reference=geometries[0].contours['blanket_outer'].points
    for a in geometries[1:]:np.testing.assert_array_equal(a.contours['blanket_outer'].points,reference)
    assert len(geometries[0].layers)<len(geometries[1].layers)
    assert any('channels/ribs' in c.name for c in geometries[2].registry.cells.values())
