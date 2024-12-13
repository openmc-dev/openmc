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

# Import a few names from the model module
from openmc.model import Model

from . import examples


__version__ = importlib.metadata.version("openmc")

try:
    OPENMC_CORE_BASE_PATH = os.path.join(__path__[0], "core")
except NameError:
    OPENMC_CORE_BASE_PATH = None

if not OPENMC_CORE_BASE_PATH or not os.path.exists(OPENMC_CORE_BASE_PATH):
    import sysconfig
    OPENMC_CORE_BASE_PATH = os.path.join(sysconfig.get_path("platlib"), "openmc", "core")
    if not os.path.exists(OPENMC_CORE_BASE_PATH):
        raise ImportError("OpenMC is not installed. Please run 'pip install openmc'.")
    warnings.warn(
        "It seems OpenMC is being run from its source directory. "
        "This setup is not recommended as it may lead to unexpected behavior, "
        "such as conflicts between source and installed versions. "
        "Please run your script from outside the OpenMC source tree.",
        RuntimeWarning
    )

def get_paths(subdir, pattern="*", recursive=False):
    """
    Helper function to return paths that match a given pattern within a subdirectory.

    Args:
        subdir (str): The subdirectory within the 'core' directory.
        pattern (str): The pattern to match files or directories.
        recursive (bool): Whether to search recursively in subdirectories.

    Returns:
        list: A list of matched paths.
    """
    search_pattern = os.path.join(OPENMC_CORE_BASE_PATH, subdir, "**", pattern) if recursive else os.path.join(OPENMC_CORE_BASE_PATH, subdir, pattern)
    return glob.glob(search_pattern, recursive=recursive) if os.path.exists(search_pattern) else []

def get_include_path():
    """Return includes and include path for OpenMC headers."""
    include = get_paths("include", "*", recursive=True)
    include_path = get_paths("include", "", recursive=False)
    return include, include_path

def get_core_libraries():
    """Return libraries and library paths for OpenMC."""
    lib = [lib_file for lib in ["lib", "lib64"] for lib_file in get_paths(lib, "libopenmc*", recursive=True)]
    lib_path = [lib_file for lib in ["lib", "lib64"] for lib_file in get_paths(lib, "", recursive=False)]
    return lib, lib_path

def get_extra_libraries():
    """Return the extra libraries installed by auditwheel or delocate."""
    libs_path = os.path.join(__path__[0], ".dylibs") if sys.platform == "darwin" else os.path.normpath(os.path.join(__path__[0], "..", "openmc.libs"))
    return (glob.glob(os.path.join(libs_path, "*")), libs_path) if os.path.exists(libs_path) else ([], [])

# Setup variables
include, include_path = get_include_path()
lib, lib_path = get_core_libraries()
extra_lib, extra_lib_path = get_extra_libraries()

# Export variables for easy access
__all__ = ["include", "include_path", "lib", "lib_path", "extra_lib", "extra_lib_path"]