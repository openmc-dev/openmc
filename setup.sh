#!/bin/bash
################################################################################
# OpenMC Installation and Setup Script
################################################################################
# This script automates the installation and configuration of OpenMC,
# including setting up cross section data paths, building the code, setting up
# a Python virtual environment, and installing the Python API.
#
# Author: William Zywiec (willzywiec@gmail.com)
#
# Usage: ./setup.sh [OPTIONS]
#
# Options:
#   --xs-dir PATH       Path to cross section data directory (default: ../endfb80-hdf5)
#   --skip-xs           Skip cross section data setup
#   --with-mpi          Build with MPI support
#   --build-type        Set build type (Debug|Release|RelWithDebInfo)
#   --force-submodules  Force re-download of vendor submodules (fixes corrupted state)
#   --help              Show this help message
################################################################################

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
SKIP_XS=false
WITH_MPI=false
BUILD_TYPE="RelWithDebInfo"
INSTALL_PREFIX="${HOME}/.local"
XS_DIR=""  # Will be set to default after SCRIPT_DIR is determined
FORCE_SUBMODULES=false

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default cross section directory (sibling to openmc directory)
DEFAULT_XS_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)/endfb80-hdf5"

################################################################################
# Helper functions
################################################################################

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_usage() {
    sed -n '2,16p' "$0" | sed 's/^# //'
}

################################################################################
# Parse command line arguments
################################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        --xs-dir)
            XS_DIR="$2"
            shift 2
            ;;
        --skip-xs)
            SKIP_XS=true
            shift
            ;;
        --with-mpi)
            WITH_MPI=true
            shift
            ;;
        --build-type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        --install-prefix)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --force-submodules)
            FORCE_SUBMODULES=true
            shift
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

################################################################################
# Check prerequisites
################################################################################

log_info "Checking prerequisites..."

# Check if we're in the OpenMC directory
if [[ ! -f "${SCRIPT_DIR}/CMakeLists.txt" ]] || [[ ! -d "${SCRIPT_DIR}/src" ]]; then
    log_error "This script must be run from the OpenMC root directory"
    exit 1
fi

# Check for CMake
if ! command -v cmake &> /dev/null; then
    log_error "CMake is required but not found"
    exit 1
fi
CMAKE_VERSION=$(cmake --version | head -n1 | awk '{print $3}')
log_info "Found CMake ${CMAKE_VERSION}"

# Check for Python
if ! command -v python3 &> /dev/null; then
    log_error "Python 3 is required but not found"
    exit 1
fi
PYTHON_VERSION=$(python3 --version | awk '{print $2}')
log_info "Found Python ${PYTHON_VERSION}"

# Check for HDF5 development libraries
HDF5_FOUND=false
if pkg-config --exists hdf5 2>/dev/null; then
    HDF5_FOUND=true
elif [[ -f "/usr/include/hdf5.h" ]] || [[ -f "/usr/local/include/hdf5.h" ]] || \
     [[ -d "/usr/include/hdf5" ]] || [[ -d "/usr/local/include/hdf5" ]]; then
    HDF5_FOUND=true
elif ldconfig -p 2>/dev/null | grep -q libhdf5; then
    HDF5_FOUND=true
fi

if [[ "${HDF5_FOUND}" == true ]]; then
    log_info "Found HDF5 development libraries"
else
    log_warning "HDF5 development libraries not detected (CMake will check during configuration)"
fi

################################################################################
# Initialize git submodules
################################################################################

log_info "Checking vendor dependencies..."

# Check if git is available
if ! command -v git &> /dev/null; then
    log_error "Git is required but not found"
    exit 1
fi

# Function to clone a vendor dependency directly (fallback for corrupted submodules)
clone_vendor_dep() {
    local name="$1"
    local url="$2"
    local tag="$3"
    local dest="${SCRIPT_DIR}/vendor/${name}"

    log_info "Cloning ${name} from ${url}..."

    # Remove existing directory with multiple fallback approaches
    if [[ -d "${dest}" ]]; then
        # First, try to fix permissions (git submodules often have restrictive perms)
        chmod -R u+rwX "${dest}" 2>/dev/null || true
        # Remove any git index lock files
        rm -f "${dest}/.git/index.lock" 2>/dev/null || true
        # Try standard removal
        rm -rf "${dest}" 2>/dev/null
        # If still exists, try with sudo or more aggressive approach
        if [[ -d "${dest}" ]]; then
            log_warning "Standard removal failed for ${dest}, trying alternative methods..."
            # Try removing contents first, then directory
            find "${dest}" -type f -exec rm -f {} \; 2>/dev/null || true
            find "${dest}" -type d -empty -delete 2>/dev/null || true
            rm -rf "${dest}" 2>/dev/null || true
        fi
        # Final check
        if [[ -d "${dest}" ]]; then
            log_error "Cannot remove existing directory ${dest}"
            log_error "Please manually remove it with: sudo rm -rf ${dest}"
            return 1
        fi
    fi

    if git clone --depth 1 --branch "${tag}" "${url}" "${dest}"; then
        rm -rf "${dest}/.git"  # Remove .git to avoid submodule conflicts
        log_success "Successfully cloned ${name} ${tag}"
        return 0
    else
        log_error "Failed to clone ${name}"
        return 1
    fi
}

# Initialize and update git submodules
cd "${SCRIPT_DIR}"
if [[ -d ".git" ]]; then
    # Force re-download if requested
    if [[ "${FORCE_SUBMODULES}" == true ]]; then
        log_warning "Force submodule mode: removing and re-downloading vendor dependencies..."
        for dep in xtl xtensor pugixml fmt Catch2; do
            rm -rf "${SCRIPT_DIR}/vendor/${dep}"
        done
        git submodule deinit -f --all 2>/dev/null || true
    fi

    log_info "Initializing git submodules..."
    if ! git submodule update --init --recursive; then
        log_warning "Standard submodule update failed, trying with force..."
        if ! git submodule update --init --recursive --force; then
            log_warning "Submodule update failed. Attempting direct clone fallback..."
            SUBMODULE_FAILED=true
        fi
    fi

    # Check if submodules are still empty and use direct clone as fallback
    # Must check for actual header files, not just directories
    NEED_FALLBACK=false
    if [[ ! -f "${SCRIPT_DIR}/vendor/xtl/include/xtl/xbasic_fixed_string.hpp" ]]; then
        log_warning "xtl headers not found after submodule update"
        NEED_FALLBACK=true
    fi
    if [[ ! -f "${SCRIPT_DIR}/vendor/xtensor/include/xtensor/xtensor.hpp" ]]; then
        log_warning "xtensor headers not found after submodule update"
        NEED_FALLBACK=true
    fi

    if [[ "${NEED_FALLBACK}" == true ]] || [[ "${SUBMODULE_FAILED:-false}" == true ]]; then
        log_warning "Submodules appear corrupted. Using direct clone fallback..."

        # Clone with specific compatible versions
        clone_vendor_dep "xtl" "https://github.com/xtensor-stack/xtl.git" "0.8.1" || exit 1
        clone_vendor_dep "xtensor" "https://github.com/xtensor-stack/xtensor.git" "0.27.1" || exit 1
        clone_vendor_dep "pugixml" "https://github.com/zeux/pugixml.git" "latest" || \
            clone_vendor_dep "pugixml" "https://github.com/zeux/pugixml.git" "v1.14" || exit 1
        clone_vendor_dep "fmt" "https://github.com/fmtlib/fmt.git" "11.0.2" || exit 1
        clone_vendor_dep "Catch2" "https://github.com/catchorg/Catch2.git" "v3.3.2" || exit 1

        log_success "Vendor dependencies cloned via fallback method"
    else
        log_success "Git submodules initialized"
    fi
else
    log_warning "Not a git repository, skipping submodule initialization"
fi

# Verify critical vendor dependencies exist and have include directories
VENDOR_DEPS=("xtl" "xtensor" "pugixml" "fmt" "Catch2")
for dep in "${VENDOR_DEPS[@]}"; do
    DEP_DIR="${SCRIPT_DIR}/vendor/${dep}"
    # Check if directory exists and has content (more than just . and ..)
    if [[ ! -d "${DEP_DIR}" ]] || [[ -z "$(ls -A "${DEP_DIR}" 2>/dev/null)" ]]; then
        log_error "Vendor dependency '${dep}' not found or empty at ${DEP_DIR}"
        log_error "Try running with --force-submodules to fix:"
        log_error "  ./setup.sh --force-submodules --skip-xs"
        exit 1
    fi
    # Extra check for header-only libraries - verify actual header files exist
    if [[ "${dep}" == "xtl" ]]; then
        if [[ ! -f "${DEP_DIR}/include/xtl/xbasic_fixed_string.hpp" ]]; then
            log_error "Vendor dependency 'xtl' is missing header files"
            log_error "Try running with --force-submodules to fix:"
            log_error "  ./setup.sh --force-submodules --skip-xs"
            exit 1
        fi
    fi
    if [[ "${dep}" == "xtensor" ]]; then
        if [[ ! -f "${DEP_DIR}/include/xtensor/xtensor.hpp" ]]; then
            log_error "Vendor dependency 'xtensor' is missing header files"
            log_error "Try running with --force-submodules to fix:"
            log_error "  ./setup.sh --force-submodules --skip-xs"
            exit 1
        fi
    fi
    log_info "Found ${dep} (verified)"
done

log_success "Vendor dependencies ready"

################################################################################
# Setup cross section data
################################################################################

if [[ "$SKIP_XS" == false ]]; then
    log_info "Setting up cross section data..."

    # Use specified XS_DIR or default to sibling directory
    if [[ -z "${XS_DIR}" ]]; then
        XS_DATA_DIR="${DEFAULT_XS_DIR}"
    else
        # Resolve to absolute path
        XS_DATA_DIR="$(cd "$(dirname "${XS_DIR}")" && pwd)/$(basename "${XS_DIR}")"
    fi

    # Verify cross section data exists
    if [[ -d "${XS_DATA_DIR}" ]]; then
        if [[ -f "${XS_DATA_DIR}/cross_sections.xml" ]]; then
            log_success "Cross section data found at ${XS_DATA_DIR}"
        else
            log_error "cross_sections.xml not found in ${XS_DATA_DIR}"
            exit 1
        fi
    else
        log_error "Cross section directory not found: ${XS_DATA_DIR}"
        log_error "Please ensure the endfb80-hdf5 data is available, or specify with --xs-dir"
        exit 1
    fi
else
    log_info "Skipping cross section data setup"
    # Still need XS_DATA_DIR for environment script
    if [[ -z "${XS_DIR}" ]]; then
        XS_DATA_DIR="${DEFAULT_XS_DIR}"
    else
        XS_DATA_DIR="$(cd "$(dirname "${XS_DIR}")" && pwd)/$(basename "${XS_DIR}")"
    fi
fi

################################################################################
# Build OpenMC
################################################################################

log_info "Building OpenMC..."

# Create build directory
BUILD_DIR="${SCRIPT_DIR}/build"
if [[ -d "${BUILD_DIR}" ]]; then
    log_info "Removing existing build directory"
    rm -rf "${BUILD_DIR}"
fi
mkdir -p "${BUILD_DIR}"

# Configure CMake
log_info "Configuring CMake (build type: ${BUILD_TYPE})..."
cd "${BUILD_DIR}"

CMAKE_OPTIONS=(
    "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
    "-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}"
)

if [[ "$WITH_MPI" == true ]]; then
    CMAKE_OPTIONS+=("-DOPENMC_USE_MPI=ON")
    log_info "MPI support enabled"
fi

cmake "${CMAKE_OPTIONS[@]}" ..

# Build
log_info "Compiling OpenMC..."
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

# Install
log_info "Installing OpenMC to ${INSTALL_PREFIX}..."
make install

log_success "OpenMC compiled and installed"

################################################################################
# Setup Python virtual environment and install OpenMC Python API
################################################################################

log_info "Setting up Python virtual environment..."

cd "${SCRIPT_DIR}"

# Create virtual environment if it doesn't exist
VENV_DIR="${SCRIPT_DIR}/.env"
if [[ ! -d "${VENV_DIR}" ]]; then
    log_info "Creating Python virtual environment at ${VENV_DIR}..."
    python3 -m venv "${VENV_DIR}"
    log_success "Virtual environment created"
else
    log_info "Virtual environment already exists at ${VENV_DIR}"
fi

# Activate virtual environment
log_info "Activating virtual environment..."
source "${VENV_DIR}/bin/activate"

# Upgrade pip in virtual environment
log_info "Upgrading pip..."
python -m pip install --upgrade pip --quiet --timeout=120 --retries=5 || {
    log_warning "Failed to upgrade pip, continuing with existing version..."
}

# Install OpenMC Python API in development mode
log_info "Installing OpenMC Python package in development mode..."
python -m pip install -e . --timeout=120 --retries=5

log_success "Python virtual environment setup complete"
log_success "OpenMC Python bindings installed"

################################################################################
# Set up environment
################################################################################

log_info "Configuring environment..."

# Create environment setup script with absolute paths
ENV_SCRIPT="${SCRIPT_DIR}/openmc_env.sh"
cat > "${ENV_SCRIPT}" << EOF
#!/bin/bash
# OpenMC environment setup script
# Source this file to set up the OpenMC environment:
#   source openmc_env.sh

# Get the directory where this script is located
SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"

# Activate Python virtual environment
if [[ -d "\${SCRIPT_DIR}/.env" ]]; then
    source "\${SCRIPT_DIR}/.env/bin/activate"
    echo "Python virtual environment activated"
else
    echo "Warning: Virtual environment not found at \${SCRIPT_DIR}/.env"
fi

# Add OpenMC binary to PATH (using absolute path)
export PATH="${INSTALL_PREFIX}/bin:\${PATH}"

# Set cross section data path (using absolute path)
export OPENMC_CROSS_SECTIONS="${XS_DATA_DIR}/cross_sections.xml"

echo "OpenMC environment configured:"
echo "  OpenMC executable: \$(which openmc 2>/dev/null || echo 'not found in PATH')"
echo "  Cross sections: \${OPENMC_CROSS_SECTIONS}"
echo "  Python: \$(which python)"
EOF

chmod +x "${ENV_SCRIPT}"

log_success "Environment script created: ${ENV_SCRIPT}"

################################################################################
# Verify installation
################################################################################

log_info "Verifying installation..."

# Check if openmc binary exists
if command -v openmc &> /dev/null; then
    OPENMC_VERSION=$(openmc --version 2>&1 || echo "unknown")
    log_success "OpenMC executable found: $(which openmc)"
    log_info "Version: ${OPENMC_VERSION}"
else
    log_warning "OpenMC executable not found in PATH"
fi

# Check Python module (should still be in activated venv)
if python -c "import openmc; print(f'OpenMC Python API version: {openmc.__version__}')" 2>/dev/null; then
    log_success "OpenMC Python module imported successfully"
else
    log_error "Failed to import OpenMC Python module"
    exit 1
fi

################################################################################
# Print summary
################################################################################

echo ""
echo "================================================================================"
log_success "OpenMC installation completed successfully!"
echo "================================================================================"
echo ""
echo "Installation summary:"
echo "  - OpenMC installed to: ${INSTALL_PREFIX}"
echo "  - Cross section data: ${XS_DATA_DIR}"
echo "  - Python virtual environment: ${VENV_DIR}"
echo "  - Build type: ${BUILD_TYPE}"
echo "  - MPI support: $([ "$WITH_MPI" == true ] && echo "Enabled" || echo "Disabled")"
echo ""
echo "To use OpenMC in a new terminal, run:"
echo "  cd ${SCRIPT_DIR}"
echo "  source openmc_env.sh"
echo ""
echo "This will:"
echo "  - Activate the Python virtual environment"
echo "  - Add OpenMC binary to PATH"
echo "  - Set OPENMC_CROSS_SECTIONS environment variable"
echo ""
echo "To test the installation, try running one of the examples:"
echo "  source openmc_env.sh"
echo "  cd examples"
echo "  python kinetics_benchmark_problem1.py"
echo ""
echo "================================================================================"
