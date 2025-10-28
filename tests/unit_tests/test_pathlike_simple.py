"""Simple test for PathLike filename support"""

from pathlib import Path
import tempfile
import os
import pytest
import openmc
from openmc.checkvalue import check_type, PathLike

# Test the type checking directly
def test_pathlike_type_checking():
    """Test that PathLike type checking works correctly"""
    
    # Test with string (should work)
    check_type('filename', 'test.txt', PathLike)
    
    # Test with Path object (should work)
    path_obj = Path('test.txt')
    check_type('filename', path_obj, PathLike)
    
    # Test with Path object containing subdirectories (should work)
    path_with_subdir = Path('subdir') / 'test.txt'
    check_type('filename', path_with_subdir, PathLike)
    
    # Test with invalid type (should raise TypeError)
    with pytest.raises(TypeError):
        check_type('filename', 123, PathLike)

def test_plot_filename_pathlike():
    """Test that plot filename accepts Path objects"""
    
    plot = openmc.Plot()
    
    # Test with string (should still work)
    plot.filename = "test_plot"
    assert plot.filename == "test_plot"
    
    # Test with Path object
    path_obj = Path("test_plot_path")
    plot.filename = path_obj
    assert plot.filename == path_obj
    
    # Test with Path object containing subdirectories
    path_with_subdir = Path("subdir") / "test_plot"
    plot.filename = path_with_subdir
    assert plot.filename == path_with_subdir


if __name__ == "__main__":
    test_pathlike_type_checking()
    print("[STATUS] PathLike type checking works correctly")
    
    test_plot_filename_pathlike()
    print("[STATUS] Plot filename PathLike functionality works")
    
    print("All tests passed!")
