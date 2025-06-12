# Install VCPKG
git clone https://github.com/microsoft/vcpkg.git "$Env:GITHUB_WORKSPACE\..\vcpkg"

cd "$Env:GITHUB_WORKSPACE\..\vcpkg"
.\bootstrap-vcpkg.bat

[Environment]::SetEnvironmentVariable("VCPKG_ROOT", "$Env:GITHUB_WORKSPACE\..\vcpkg", [System.EnvironmentVariableTarget]::User)
$Env:VCPKG_ROOT = "$Env:GITHUB_WORKSPACE\..\vcpkg"
$Env:Path += ";$Env:VCPKG_ROOT"

# Install HDF5
vcpkg install hdf5:x64-windows-static

# Go back to main directory for install
cd "$Env:GITHUB_WORKSPACE"

# Upgrade pip, pytest, numpy before doing anything else.
pip install --upgrade pip
pip install --upgrade pytest
pip install --upgrade numpy

# Build and install OpenMC executable
python tools/ci/gha-install.py

# Install Python API in editable mode
pip install -e .[test,vtk,ci]