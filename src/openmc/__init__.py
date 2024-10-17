import os
import glob
import importlib.metadata
from openmc.arithmetic import *
from openmc.bounding_box import *
from openmc.cell import *
from openmc.checkvalue import *
from openmc.mesh import *
from openmc.element import *
from openmc.geometry import *
from openmc.nuclide import *
from openmc.macroscopic import *
from openmc.material import *
from openmc.plots import *
from openmc.region import *
from openmc.volume import *
from openmc.weight_windows import *
from openmc.surface import *
from openmc.universe import *
from openmc.source import *
from openmc.settings import *
from openmc.lattice import *
from openmc.filter import *
from openmc.filter_expansion import *
from openmc.trigger import *
from openmc.tally_derivative import *
from openmc.tallies import *
from openmc.mgxs_library import *
from openmc.executor import *
from openmc.statepoint import *
from openmc.summary import *
from openmc.particle_restart import *
from openmc.mixin import *
from openmc.plotter import *
from openmc.search import *
from openmc.polynomial import *
from openmc.tracks import *
from .config import *
from .openmc_exec import main

# Import a few names from the model module
from openmc.model import Model

from . import examples


__version__ = importlib.metadata.version("openmc")

def get_path(subdir, pattern=""):
    """Helper function to return paths that match a given pattern within a subdirectory."""
    path = os.path.join(__path__[0], "core", subdir)
    return glob.glob(os.path.join(path, pattern)) if os.path.exists(path) else []

def get_include_path():
    """Return the include directory path for OpenMC headers."""
    return get_path("include")

def get_core_libraries():
    """Return library paths and library directory paths."""
    lib_paths = [libs for lib in ["lib", "lib64"] for libs in get_path(lib)]
    lib = [libs for lib in ["lib", "lib64"] for libs in get_path(lib, pattern="*openmc*")]
    return lib, lib_paths

def get_extra_libraries():
    """List all the extra libraries of OpenMC."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.join(__path__[0], "..", "pyne.libs")
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]
