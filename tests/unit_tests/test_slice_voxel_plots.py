"""Tests for SlicePlot and VoxelPlot classes

This module tests the functionality of the new SlicePlot and VoxelPlot
classes that replace the legacy Plot class.
"""
import warnings
from pathlib import Path

import numpy as np
import pytest

import openmc
from openmc.plots import _SVG_COLORS


def test_slice_plot_initialization():
    """Test SlicePlot initialization with defaults"""
    plot = openmc.SlicePlot()
    assert plot.width == [4.0, 4.0]
    assert plot.pixels == [400, 400]
    assert plot.basis == 'xy'
    assert plot.origin == [0., 0., 0.]


def test_slice_plot_width_validation():
    """Test that SlicePlot only accepts 2 values for width"""
    plot = openmc.SlicePlot()
    
    # Should accept 2 values
    plot.width = [10.0, 20.0]
    assert plot.width == [10.0, 20.0]
    
    # Should reject 1 value
    with pytest.raises(ValueError, match='must be of length "2"'):
        plot.width = [10.0]
    
    # Should reject 3 values
    with pytest.raises(ValueError, match='must be of length "2"'):
        plot.width = [10.0, 20.0, 30.0]


def test_slice_plot_pixels_validation():
    """Test that SlicePlot only accepts 2 values for pixels"""
    plot = openmc.SlicePlot()
    
    # Should accept 2 values
    plot.pixels = [100, 200]
    assert plot.pixels == [100, 200]
    
    # Should reject 1 value
    with pytest.raises(ValueError, match='must be of length "2"'):
        plot.pixels = [100]
    
    # Should reject 3 values
    with pytest.raises(ValueError, match='must be of length "2"'):
        plot.pixels = [100, 200, 300]


def test_slice_plot_basis():
    """Test that SlicePlot has basis attribute"""
    plot = openmc.SlicePlot()
    
    # Test all valid basis values
    for basis in ['xy', 'xz', 'yz']:
        plot.basis = basis
        assert plot.basis == basis
    
    # Test invalid basis
    with pytest.raises(ValueError):
        plot.basis = 'invalid'


def test_slice_plot_meshlines():
    """Test that SlicePlot has meshlines attribute"""
    plot = openmc.SlicePlot()
    
    meshlines = {
        'type': 'tally',
        'id': 1,
        'linewidth': 2,
        'color': (255, 0, 0)
    }
    plot.meshlines = meshlines
    assert plot.meshlines == meshlines


def test_slice_plot_xml_roundtrip():
    """Test SlicePlot XML serialization and deserialization"""
    plot = openmc.SlicePlot(name='test_slice')
    plot.width = [15.0, 25.0]
    plot.pixels = [150, 250]
    plot.basis = 'xz'
    plot.origin = [1.0, 2.0, 3.0]
    plot.color_by = 'material'
    plot.filename = 'test_plot'
    
    # Convert to XML and back
    elem = plot.to_xml_element()
    new_plot = openmc.SlicePlot.from_xml_element(elem)
    
    # Check all attributes preserved
    assert new_plot.name == plot.name
    assert new_plot.width == pytest.approx(plot.width)
    assert new_plot.pixels == tuple(plot.pixels)
    assert new_plot.basis == plot.basis
    assert new_plot.origin == pytest.approx(plot.origin)
    assert new_plot.color_by == plot.color_by
    assert new_plot.filename == plot.filename


def test_slice_plot_from_geometry():
    """Test creating SlicePlot from geometry"""
    # Create simple geometry
    s = openmc.Sphere(r=10.0, boundary_type='vacuum')
    c = openmc.Cell(region=-s)
    univ = openmc.Universe(cells=[c])
    geom = openmc.Geometry(univ)
    
    # Test all basis options
    for basis in ['xy', 'xz', 'yz']:
        plot = openmc.SlicePlot.from_geometry(geom, basis=basis)
        assert plot.basis == basis
        assert plot.width == pytest.approx([20.0, 20.0])
        assert plot.origin == pytest.approx([0.0, 0.0, 0.0])


def test_voxel_plot_initialization():
    """Test VoxelPlot initialization with defaults"""
    plot = openmc.VoxelPlot()
    assert plot.width == [4.0, 4.0, 4.0]
    assert plot.pixels == [400, 400, 400]
    assert plot.origin == [0., 0., 0.]


def test_voxel_plot_width_validation():
    """Test that VoxelPlot only accepts 3 values for width"""
    plot = openmc.VoxelPlot()
    
    # Should accept 3 values
    plot.width = [10.0, 20.0, 30.0]
    assert plot.width == [10.0, 20.0, 30.0]
    
    # Should reject 2 values
    with pytest.raises(ValueError, match='must be of length "3"'):
        plot.width = [10.0, 20.0]
    
    # Should reject 1 value
    with pytest.raises(ValueError, match='must be of length "3"'):
        plot.width = [10.0]


def test_voxel_plot_pixels_validation():
    """Test that VoxelPlot only accepts 3 values for pixels"""
    plot = openmc.VoxelPlot()
    
    # Should accept 3 values
    plot.pixels = [100, 200, 300]
    assert plot.pixels == [100, 200, 300]
    
    # Should reject 2 values
    with pytest.raises(ValueError, match='must be of length "3"'):
        plot.pixels = [100, 200]
    
    # Should reject 1 value
    with pytest.raises(ValueError, match='must be of length "3"'):
        plot.pixels = [100]


def test_voxel_plot_no_basis():
    """Test that VoxelPlot does not have basis attribute"""
    plot = openmc.VoxelPlot()
    assert not hasattr(plot, 'basis')
    assert not hasattr(plot, '_basis')


def test_voxel_plot_no_meshlines():
    """Test that VoxelPlot does not have meshlines attribute"""
    plot = openmc.VoxelPlot()
    assert not hasattr(plot, 'meshlines')
    assert not hasattr(plot, '_meshlines')


def test_voxel_plot_has_to_vtk():
    """Test that VoxelPlot has to_vtk method"""
    plot = openmc.VoxelPlot()
    assert hasattr(plot, 'to_vtk')
    assert callable(plot.to_vtk)


def test_voxel_plot_xml_roundtrip():
    """Test VoxelPlot XML serialization and deserialization"""
    plot = openmc.VoxelPlot(name='test_voxel')
    plot.width = [10.0, 20.0, 30.0]
    plot.pixels = [100, 200, 300]
    plot.origin = [1.0, 2.0, 3.0]
    plot.color_by = 'cell'
    plot.filename = 'voxel_plot'
    
    # Convert to XML and back
    elem = plot.to_xml_element()
    new_plot = openmc.VoxelPlot.from_xml_element(elem)
    
    # Check all attributes preserved
    assert new_plot.name == plot.name
    assert new_plot.width == pytest.approx(plot.width)
    assert new_plot.pixels == tuple(plot.pixels)
    assert new_plot.origin == pytest.approx(plot.origin)
    assert new_plot.color_by == plot.color_by
    assert new_plot.filename == plot.filename


def test_plot_deprecation_warning():
    """Test that Plot class raises deprecation warning"""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        plot = openmc.Plot()
        
        assert len(w) == 1
        assert issubclass(w[0].category, DeprecationWarning)
        assert "deprecated" in str(w[0].message).lower()


def test_plot_slice_compatibility():
    """Test Plot class with slice type"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot = openmc.Plot()
        
        plot.type = 'slice'
        plot.width = [10.0, 20.0]
        plot.pixels = [100, 200]
        plot.basis = 'yz'
        
        assert plot.type == 'slice'
        assert plot.width == [10.0, 20.0]
        assert plot.pixels == [100, 200]
        assert plot.basis == 'yz'


def test_plot_voxel_compatibility():
    """Test Plot class with voxel type"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot = openmc.Plot()
        
        plot.type = 'voxel'
        plot.width = [10.0, 20.0, 30.0]
        plot.pixels = [100, 200, 300]
        
        assert plot.type == 'voxel'
        assert plot.width == [10.0, 20.0, 30.0]
        assert plot.pixels == [100, 200, 300]


def test_plot_xml_roundtrip_slice():
    """Test XML roundtrip for Plot with slice type"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot = openmc.Plot(name='legacy_slice')
        plot.type = 'slice'
        plot.width = [5.0, 10.0]
        plot.pixels = [50, 100]
        plot.basis = 'xz'
        
        elem = plot.to_xml_element()
        new_plot = openmc.Plot.from_xml_element(elem)
        
        assert new_plot.type == plot.type
        assert new_plot.width == pytest.approx(plot.width)
        assert new_plot.pixels == tuple(plot.pixels)
        assert new_plot.basis == plot.basis


def test_plot_xml_roundtrip_voxel():
    """Test XML roundtrip for Plot with voxel type"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plot = openmc.Plot(name='legacy_voxel')
        plot.type = 'voxel'
        plot.width = [5.0, 10.0, 15.0]
        plot.pixels = [50, 100, 150]
        
        elem = plot.to_xml_element()
        new_plot = openmc.Plot.from_xml_element(elem)
        
        assert new_plot.type == plot.type
        assert new_plot.width == pytest.approx(plot.width)
        assert new_plot.pixels == tuple(plot.pixels)


def test_plots_collection_mixed_types():
    """Test Plots collection with different plot types"""
    slice_plot = openmc.SlicePlot(name='slice')
    voxel_plot = openmc.VoxelPlot(name='voxel')
    wireframe_plot = openmc.WireframeRayTracePlot(name='wireframe')
    
    plots = openmc.Plots([slice_plot, voxel_plot, wireframe_plot])
    
    assert len(plots) == 3
    assert isinstance(plots[0], openmc.SlicePlot)
    assert isinstance(plots[1], openmc.VoxelPlot)
    assert isinstance(plots[2], openmc.WireframeRayTracePlot)


def test_plots_collection_xml_roundtrip(run_in_tmpdir):
    """Test XML export and import with new plot types"""
    s1 = openmc.SlicePlot(name='slice1')
    s1.width = [10.0, 20.0]
    s1.basis = 'xz'
    
    v1 = openmc.VoxelPlot(name='voxel1')
    v1.width = [10.0, 20.0, 30.0]
    
    plots = openmc.Plots([s1, v1])
    plots.export_to_xml()
    
    # Read back
    new_plots = openmc.Plots.from_xml()
    
    assert len(new_plots) == 2
    assert isinstance(new_plots[0], openmc.SlicePlot)
    assert isinstance(new_plots[1], openmc.VoxelPlot)
    assert new_plots[0].name == 'slice1'
    assert new_plots[1].name == 'voxel1'
    assert new_plots[0].basis == 'xz'
    assert new_plots[0].width == pytest.approx([10.0, 20.0])
    assert new_plots[1].width == pytest.approx([10.0, 20.0, 30.0])