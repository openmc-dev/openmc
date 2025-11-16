#!/bin/bash
################################################################################
# OpenMC Installation and Setup Script
################################################################################
# This script automates the installation and configuration of OpenMC,
# including downloading cross-section data, installing dependencies,
# building the code, and setting up the Python API.
#
# Usage: ./setup.sh [OPTIONS]
#
# Options:
#   --skip-deps       Skip system dependency installation
#   --skip-xs         Skip cross-section data download
#   --with-mpi        Build with MPI support
#   --build-type      Set build type (Debug|Release|RelWithDebInfo)
#   --help            Show this help message
################################################################################

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
SKIP_DEPS=false
SKIP_XS=false
WITH_MPI=false
BUILD_TYPE="RelWithDebInfo"
INSTALL_PREFIX="${HOME}/.local"

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
    sed -n '2,17p' "$0" | sed 's/^# //'
}

################################################################################
# Parse command line arguments
################################################################################

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-deps)
            SKIP_DEPS=true
            shift
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

# Check for Python
if ! command -v python3 &> /dev/null; then
    log_error "Python 3 is required but not found"
    exit 1
fi

PYTHON_VERSION=$(python3 --version | awk '{print $2}')
log_info "Found Python ${PYTHON_VERSION}"

################################################################################
# Install system dependencies
################################################################################

if [[ "$SKIP_DEPS" == false ]]; then
    log_info "Installing system dependencies..."

    # Detect package manager
    if command -v apt-get &> /dev/null; then
        PKG_MANAGER="apt-get"
        PACKAGES="g++ cmake libhdf5-dev libpng-dev git"

        if [[ "$WITH_MPI" == true ]]; then
            PACKAGES="${PACKAGES} mpich libmpich-dev libhdf5-mpich-dev"
        fi

        log_info "Using apt-get to install: ${PACKAGES}"
        sudo apt-get update
        sudo apt-get install -y ${PACKAGES}

    elif command -v dnf &> /dev/null; then
        PKG_MANAGER="dnf"
        PACKAGES="gcc-c++ cmake hdf5-devel libpng-devel git"

        if [[ "$WITH_MPI" == true ]]; then
            PACKAGES="${PACKAGES} mpich-devel hdf5-mpich-devel"
        fi

        log_info "Using dnf to install: ${PACKAGES}"
        sudo dnf install -y ${PACKAGES}

    elif command -v brew &> /dev/null; then
        PKG_MANAGER="brew"
        PACKAGES="cmake hdf5 libpng"

        if [[ "$WITH_MPI" == true ]]; then
            PACKAGES="${PACKAGES} mpich"
        fi

        log_info "Using Homebrew to install: ${PACKAGES}"
        brew install ${PACKAGES}

    else
        log_warning "No supported package manager found. Please install dependencies manually:"
        log_warning "  - C++ compiler (g++)"
        log_warning "  - CMake"
        log_warning "  - HDF5 library"
        log_warning "  - libpng (optional, for plotting)"
        if [[ "$WITH_MPI" == true ]]; then
            log_warning "  - MPI implementation (MPICH or OpenMPI)"
        fi
    fi

    log_success "System dependencies installed"
else
    log_info "Skipping system dependency installation"
fi

################################################################################
# Download cross-section data
################################################################################

if [[ "$SKIP_XS" == false ]]; then
    log_info "Downloading cross-section data..."

    # Create data directory
    XS_DATA_DIR="${HOME}/openmc_data"
    mkdir -p "${XS_DATA_DIR}"

    # Download NNDC HDF5 cross sections
    if [[ ! -e "${XS_DATA_DIR}/nndc_hdf5/cross_sections.xml" ]]; then
        log_info "Downloading NNDC HDF5 cross-section library..."
        wget -q --show-progress -O - \
            https://anl.box.com/shared/static/teaup95cqv8s9nn56hfn7ku8mmelr95p.xz \
            | tar -C "${XS_DATA_DIR}" -xJ
        log_success "NNDC HDF5 data downloaded"
    else
        log_info "NNDC HDF5 data already exists, skipping download"
    fi

    # Download ENDF/B-VII.1 distribution
    ENDF_DIR="${XS_DATA_DIR}/endf-b-vii.1"
    if [[ ! -d "${ENDF_DIR}/neutrons" ]] || \
       [[ ! -d "${ENDF_DIR}/photoat" ]] || \
       [[ ! -d "${ENDF_DIR}/atomic_relax" ]]; then
        log_info "Downloading ENDF/B-VII.1 distribution..."
        wget -q --show-progress -O - \
            https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz \
            | tar -C "${XS_DATA_DIR}" -xJ
        log_success "ENDF/B-VII.1 data downloaded"
    else
        log_info "ENDF/B-VII.1 data already exists, skipping download"
    fi

    # Copy ENDF/B-VIII.0 HDF5 data if available
    if [[ -d "endfb80-hdf5" ]]; then
        log_info "Copying ENDF/B-VIII.0 HDF5 cross-section data..."
        mkdir -p "${SCRIPT_DIR}/openmc"
        cp -r endfb80-hdf5 "${SCRIPT_DIR}/openmc/endfb80-hdf5"
        log_success "ENDF/B-VIII.0 HDF5 data copied to ${SCRIPT_DIR}/openmc/endfb80-hdf5"
    elif [[ -d "${SCRIPT_DIR}/openmc/endfb80-hdf5" ]]; then
        log_info "ENDF/B-VIII.0 HDF5 data already exists at ${SCRIPT_DIR}/openmc/endfb80-hdf5"
    else
        log_info "ENDF/B-VIII.0 HDF5 data not found, skipping (optional)"
    fi

    log_success "Cross-section data ready at ${XS_DATA_DIR}"
else
    log_info "Skipping cross-section data download"
    XS_DATA_DIR="${HOME}/openmc_data"
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
# Install Python API
################################################################################

log_info "Installing Python API..."

cd "${SCRIPT_DIR}"

# Upgrade pip
python3 -m pip install --upgrade pip

# Install Python dependencies
log_info "Installing Python dependencies..."
python3 -m pip install \
    numpy \
    scipy \
    h5py \
    matplotlib \
    pandas \
    lxml \
    uncertainties \
    ipython \
    setuptools \
    endf

# Install optional dependencies
log_info "Installing optional Python packages..."
python3 -m pip install mcpl ncrystal || log_warning "Some optional packages failed to install"

# Install OpenMC Python API in development mode
log_info "Installing OpenMC Python package..."
python3 -m pip install -e .

log_success "Python API installed"

################################################################################
# Set up environment
################################################################################

log_info "Configuring environment..."

# Create environment setup script
ENV_SCRIPT="${SCRIPT_DIR}/openmc_env.sh"
cat > "${ENV_SCRIPT}" << EOF
#!/bin/bash
# OpenMC environment setup script
# Source this file to set up the OpenMC environment:
#   source ${ENV_SCRIPT}

# Add OpenMC binary to PATH
export PATH="${INSTALL_PREFIX}/bin:\${PATH}"

# Set cross-section data paths
export OPENMC_CROSS_SECTIONS="${XS_DATA_DIR}/nndc_hdf5/cross_sections.xml"
export OPENMC_ENDF_DATA="${XS_DATA_DIR}/endf-b-vii.1"

echo "OpenMC environment configured:"
echo "  OpenMC executable: \$(which openmc 2>/dev/null || echo 'not found in PATH')"
echo "  Cross sections: \${OPENMC_CROSS_SECTIONS}"
echo "  ENDF data: \${OPENMC_ENDF_DATA}"
EOF

chmod +x "${ENV_SCRIPT}"

log_success "Environment script created: ${ENV_SCRIPT}"

################################################################################
# Verify installation
################################################################################

log_info "Verifying installation..."

# Source the environment
source "${ENV_SCRIPT}"

# Check if openmc binary exists
if command -v openmc &> /dev/null; then
    OPENMC_VERSION=$(openmc --version 2>&1 || echo "unknown")
    log_success "OpenMC executable found: $(which openmc)"
    log_info "Version: ${OPENMC_VERSION}"
else
    log_warning "OpenMC executable not found in PATH"
fi

# Check Python module
if python3 -c "import openmc; print(f'OpenMC Python API version: {openmc.__version__}')" 2>/dev/null; then
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
echo "  - Cross-section data: ${XS_DATA_DIR}"
echo "  - Build type: ${BUILD_TYPE}"
echo "  - MPI support: $([ "$WITH_MPI" == true ] && echo "Enabled" || echo "Disabled")"
echo ""
echo "To use OpenMC, run:"
echo "  source ${ENV_SCRIPT}"
echo ""
echo "Or add the following to your ~/.bashrc or ~/.zshrc:"
echo "  source ${ENV_SCRIPT}"
echo ""
echo "To test the installation, try running one of the examples:"
echo "  cd ${SCRIPT_DIR}/examples/pincell"
echo "  python build_xml.py"
echo "  openmc"
echo ""
echo "================================================================================"
