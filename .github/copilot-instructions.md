# OpenMC AI Coding Agent Instructions

## Project Overview

OpenMC is a Monte Carlo particle transport code for nuclear reactor physics simulations. It's a hybrid C++17/Python codebase where:
- **C++ core** (`src/`, `include/openmc/`) handles the computationally intensive transport simulation
- **Python API** (`openmc/`) provides user-facing model building, post-processing, and depletion capabilities
- **C API bindings** (`openmc/lib/`) wrap the C++ library via ctypes for runtime control

## Architecture & Key Components

### C++ Component Structure
- **Global vectors of unique_ptrs**: Core objects like `model::cells`, `model::universes`, `nuclides` are stored as `vector<unique_ptr<T>>` in nested namespaces (`openmc::model`, `openmc::simulation`, `openmc::settings`, `openmc::data`)
- **Custom container types**: OpenMC provides its own `vector`, `array`, `unique_ptr`, and `make_unique` in the `openmc::` namespace (defined in `vector.h`, `array.h`, `memory.h`). These are currently typedefs to `std::` equivalents but may become custom implementations for accelerator support. Always use `openmc::vector`, not `std::vector`.
- **Geometry systems**: 
  - **CSG (default)**: Arbitrarily complex Constructive Solid Geometry using `Cell`, `Surface`, `Universe`, `Lattice`
  - **DAGMC**: CAD-based geometry via Direct Accelerated Geometry Monte Carlo (optional, requires `OPENMC_USE_DAGMC`)
  - **Unstructured mesh**: libMesh-based geometry (optional, requires `OPENMC_USE_LIBMESH`)
- **Particle tracking**: `Particle` class with `GeometryState` manages particle transport through geometry
- **Tallies**: Score quantities during simulation via `Filter` and `Tally` objects
- **Random ray solver**: Alternative deterministic method in `src/random_ray/`
- **Optional features**: DAGMC (CAD geometry), libMesh (unstructured mesh), MPI, all controlled by `#ifdef OPENMC_MPI`, etc.

### Python Component Structure
- **ID management**: All geometry objects (Cell, Surface, Material, etc.) inherit from `IDManagerMixin` which auto-assigns unique integer IDs and tracks them via class-level `used_ids` and `next_id`
- **Input validation**: Extensive use of `openmc.checkvalue` module functions (`check_type`, `check_value`, `check_length`) for all setters
- **XML I/O**: Most classes implement `to_xml_element()` and `from_xml_element()` for serialization to OpenMC's XML input format
- **HDF5 output**: Post-simulation data in statepoint files read via `openmc.StatePoint`
- **Depletion**: `openmc.deplete` implements burnup via operator-splitting with various integrators (Predictor, CECM, etc.)

## Critical Build & Test Workflows

### Build Dependencies
- **CMake** (3.16+): Required for configuring and building the C++ library
- **HDF5**: Required for cross section data and output file formats
- **C++17 compiler**: GCC, Clang, or Intel

Without CMake and HDF5, OpenMC cannot be compiled.

### Building the C++ Library
```bash
# Configure with CMake (from build/ directory)
cmake .. -DOPENMC_USE_MPI=ON -DOPENMC_USE_OPENMP=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo

# Available CMake options (all default OFF except OPENMC_USE_OPENMP and OPENMC_BUILD_TESTS):
# -DOPENMC_USE_OPENMP=ON/OFF        # OpenMP parallelism
# -DOPENMC_USE_MPI=ON/OFF           # MPI support
# -DOPENMC_USE_DAGMC=ON/OFF         # CAD geometry support
# -DOPENMC_USE_LIBMESH=ON/OFF       # Unstructured mesh
# -DOPENMC_ENABLE_PROFILE=ON/OFF    # Profiling flags
# -DOPENMC_ENABLE_COVERAGE=ON/OFF   # Coverage analysis

# Build
make -j

# C++ unit tests (uses Catch2)
ctest
```

### Python Development
```bash
# Install in development mode (requires building C++ library first)
pip install -e .

# Python tests (uses pytest)
pytest tests/unit_tests/           # Fast unit tests
pytest tests/regression_tests/     # Full regression suite (requires nuclear data)
```

### Nuclear Data Setup (CRITICAL for Testing)
Most tests require the NNDC HDF5 nuclear cross-section library.

**Important**: Check if `OPENMC_CROSS_SECTIONS` is already set in the user's environment before downloading, as many users already have nuclear data installed.

**If not already configured, download and setup:**
```bash
# Download NNDC HDF5 cross section library (~800 MB compressed)
wget -q -O - https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz | tar -C $HOME -xJ

# Set environment variable (add to ~/.bashrc or ~/.zshrc for persistence)
export OPENMC_CROSS_SECTIONS=$HOME/nndc_hdf5/cross_sections.xml
```

**Alternative**: Use the provided download script (checks if data exists before downloading):
```bash
bash tools/ci/download-xs.sh  # Downloads both NNDC HDF5 and ENDF/B-VII.1 data
```

Without this data, regression tests will fail with "No cross_sections.xml file found" errors. The `cross_sections.xml` file is an index listing paths to individual HDF5 nuclear data files for each nuclide.

## Code Style & Conventions

### C++ Style (enforced by .clang-format)
- **Naming**: 
  - Classes: `CamelCase` (e.g., `HexLattice`)
  - Functions/methods: `snake_case` (e.g., `get_indices`)
  - Variables: `snake_case` with trailing underscore for class members (e.g., `n_particles_`, `energy_`)
  - Constants: `UPPER_SNAKE_CASE` (e.g., `SQRT_PI`)
- **Namespaces**: All code in `openmc::` namespace, global state in sub-namespaces
- **Include order**: Related header first, then C/C++ stdlib, third-party libs, local headers
- **Comments**: C++-style (`//`) only, never C-style (`/* */`)
- **Standard**: C++17 features allowed
- **Formatting**: Run `clang-format` (version 15) before committing; install via `tools/dev/install-commit-hooks.sh`

### Python Style
- **PEP8** compliant
- **Docstrings**: numpydoc format for all public functions/methods
- **Type hints**: Use sparingly, primarily for complex signatures
- **Path handling**: Use `pathlib.Path` for filesystem operations, accept `str | os.PathLike` in function arguments
- **Dependencies**: Core dependencies only (numpy, scipy, h5py, pandas, matplotlib, lxml, ipython, uncertainties, setuptools, endf). Other packages must be optional
- **Python version**: Minimum 3.11 (as of Nov 2025)

### ID Management Pattern (Python)
When creating geometry objects, IDs can be auto-assigned or explicit:
```python
# Auto-assigned ID
cell = openmc.Cell()  # Gets next available ID

# Explicit ID
cell = openmc.Cell(id=10)  # Warning if ID already used

# Reset all IDs (useful in test fixtures)
openmc.reset_auto_ids()
```

### Input Validation Pattern (Python)
All setters use checkvalue functions:
```python
import openmc.checkvalue as cv

@property
def temperature(self):
    return self._temperature

@temperature.setter
def temperature(self, temp):
    cv.check_type('temperature', temp, Real)
    cv.check_greater_than('temperature', temp, 0.0)
    self._temperature = temp
```

### Working with HDF5 Files
C++ uses custom HDF5 wrappers in `src/hdf5_interface.cpp`. Python uses h5py directly. Statepoint format version is `VERSION_STATEPOINT` in `include/openmc/constants.h`.

### Conditional Compilation
Check for optional features:
```cpp
#ifdef OPENMC_MPI
  // MPI-specific code
#endif

#ifdef OPENMC_DAGMC
  // DAGMC-specific code
#endif
```

## Testing Expectations

### C++ Tests
Located in `tests/cpp_unit_tests/`, use Catch2 framework. Run via `ctest` after building with `-DOPENMC_BUILD_TESTS=ON`.

### Python Unit Tests
Located in `tests/unit_tests/`, these are fast, standalone tests that verify Python API functionality without running full simulations. Use standard pytest patterns:

**Categories**:
- **API validation**: Test object creation, property setters/getters, XML serialization (e.g., `test_material.py`, `test_cell.py`, `test_source.py`)
- **Data processing**: Test nuclear data handling, cross sections, depletion chains (e.g., `test_data_neutron.py`, `test_deplete_chain.py`)
- **Library bindings**: Test `openmc.lib` ctypes interface with `model.init_lib()`/`model.finalize_lib()` (e.g., `test_lib.py`)
- **Geometry operations**: Test bounding boxes, containment, lattice generation (e.g., `test_bounding_box.py`, `test_lattice.py`)

**Common patterns**:
- Use fixtures from `tests/unit_tests/conftest.py` (e.g., `uo2`, `water`, `sphere_model`)
- Test invalid inputs with `pytest.raises(ValueError)` or `pytest.raises(TypeError)`
- Use `run_in_tmpdir` fixture for tests that create files
- Tests with `openmc.lib` require calling `model.init_lib()` in try/finally with `model.finalize_lib()`

**Example**:
```python
def test_material_properties():
    m = openmc.Material()
    m.add_nuclide('U235', 1.0)
    assert 'U235' in m.nuclides
    
    with pytest.raises(TypeError):
        m.add_nuclide('H1', '1.0')  # Invalid type
```

Unit tests should be fast and not require nuclear data or full OpenMC runs. For tests requiring simulation output, use regression tests instead.

### Python Regression Tests
Regression tests compare OpenMC output against reference data. **Prefer using existing models from `openmc.examples`** (like `pwr_pin_cell()`, `pwr_assembly()`, `slab_mg()`) rather than building from scratch.

**Test Harness Types** (in `tests/testing_harness.py`):
- **PyAPITestHarness**: Standard harness for Python API tests. Compares `inputs_true.dat` (XML hash) and `results_true.dat` (statepoint k-eff and tally values). Requires `model.xml` generation.
- **HashedPyAPITestHarness**: Like PyAPITestHarness but hashes the results for compact comparison
- **TolerantPyAPITestHarness**: For tests with floating-point non-associativity (e.g., random ray solver with single precision). Uses relative tolerance comparisons.
- **WeightWindowPyAPITestHarness**: Compares weight window bounds from `weight_windows.h5`
- **CollisionTrackTestHarness**: Compares collision track data from `collision_track.h5` against `collision_track_true.h5`
- **TestHarness**: Base harness for XML-based tests (no Python model building)
- **PlotTestHarness**: Compares plot output files (PNG or voxel HDF5)
- **CMFDTestHarness**: Specialized for CMFD acceleration tests
- **ParticleRestartTestHarness**: Tests particle restart functionality

**Example Test**:
```python
from openmc.examples import pwr_pin_cell
from tests.testing_harness import PyAPITestHarness

def test_my_feature():
    model = pwr_pin_cell()
    model.settings.particles = 1000  # Modify to exercise feature
    harness = PyAPITestHarness('statepoint.10.h5', model)
    harness.main()
```

**Workflow**: Create `test.py` and `__init__.py` in `tests/regression_tests/my_test/`, run `pytest --update` to generate reference files (`inputs_true.dat`, `results_true.dat`, etc.), then verify with `pytest` without `--update`.

**Critical**: When modifying OpenMC code, regenerate affected test references with `pytest --update` and commit updated reference files.

### Test Configuration
`pytest.ini` sets: `python_files = test*.py`, `python_classes = NoThanks` (disables class-based test collection).

## Cross-Language Boundaries

The C API (defined in `include/openmc/capi.h`) exposes C++ functionality to Python via ctypes bindings in `openmc/lib/`. Example:
```cpp
// C++ API in capi.h
extern "C" int openmc_run();

// Python binding in openmc/lib/core.py
_dll.openmc_run.restype = c_int
def run():
    _dll.openmc_run()
```

When modifying C++ public APIs, update corresponding ctypes signatures in `openmc/lib/*.py`.

## Documentation

- **User docs**: Sphinx documentation in `docs/source/` hosted at docs.openmc.org
- **C++ docs**: Doxygen-style comments with `\brief`, `\param` tags
- **Python docs**: numpydoc format docstrings

## Common Pitfalls

1. **Forgetting nuclear data**: Tests fail without `OPENMC_CROSS_SECTIONS` environment variable
2. **ID conflicts**: Python objects with duplicate IDs trigger `IDWarning`, use `reset_auto_ids()` between tests
3. **MPI builds**: Code must work with and without MPI; use `#ifdef OPENMC_MPI` guards
4. **Path handling**: Use `pathlib.Path` in new Python code, not `os.path`
5. **Clang-format version**: CI uses version 15; other versions may produce different formatting
