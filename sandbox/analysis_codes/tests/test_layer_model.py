"""Alternating-cube topology, material settings and beam placement."""
from pathlib import Path
import numpy as np
import openmc
import pytest
from sandbox.analysis_codes.context import load_study

ROOT=Path(__file__).resolve().parents[2]/'studies'


def test_cube_has_twenty_nonoverlapping_layers():
    study=load_study(ROOT/'Be_Li_layers');m=study.model
    layout=m.layer_layout()
    assert len(layout)==20
    assert sum(name=='Li' for _,_,name in layout)==10
    assert sum(name=='Be' for _,_,name in layout)==10
    assert sum(right-left for left,right,_ in layout)==pytest.approx(20)
    openmc.reset_auto_ids();model=m.build_model(25,4)
    for j,(left,right,name) in enumerate(layout):
        point=np.zeros(3);point[0]=(left+right)/2
        cell=model.geometry.find(point)[-1]
        assert cell.name==f'Layer {j:02d} {name}'
        assert cell.fill.name==('Lithium' if name=='Li' else 'Beryllium')
    assert model.geometry.find(m.source_point())[-1].fill is None
    assert model.geometry.find((.5,0,-9.99))[-1].fill.name=='Lithium'
    assert not model.geometry.find((12,0,0))


@pytest.mark.parametrize('layer_axis',['x','y','z'])
@pytest.mark.parametrize('beam_axis',[0,1,2])
def test_automatic_beam_matches_configured_axes(layer_axis,beam_axis):
    study=load_study(ROOT/'Be_Li_layers');c=study.config
    c.LAYER_AXIS=layer_axis;c.SOURCE_DIRECTION=tuple(1 if j==beam_axis else 0 for j in range(3))
    point=np.asarray(study.model.source_point())
    assert point[beam_axis]==pytest.approx(-10.01)
    openmc.reset_auto_ids();model=study.model.build_model()
    entered=point.copy();entered[beam_axis]= -9.99
    assert model.geometry.find(entered)[-1].fill is not None


def test_layer_widths_and_explicit_source():
    study=load_study(ROOT/'Be_Li_layers');c=study.config
    c.CUBE_SIDE_CM=5;c.LI_LAYER_THICKNESS_CM=2;c.BE_LAYER_THICKNESS_CM=1
    assert study.model.layer_layout()==[(-2.5,-.5,'Li'),(-.5,.5,'Be'),(.5,2.5,'Li')]
    c.SOURCE_POSITION_CM=(0,0,-3)
    assert study.model.source_point()==(0,0,-3)
