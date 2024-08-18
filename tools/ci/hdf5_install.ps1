git clone https://github.com/microsoft/vcpkg.git "$Env:GITHUB_WORKSPACE\vcpkg"

cd "$Env:GITHUB_WORKSPACE\vcpkg"
.\bootstrap-vcpkg.bat

[Environment]::SetEnvironmentVariable("VCPKG_ROOT", "$Env:GITHUB_WORKSPACE\vcpkg", [System.EnvironmentVariableTarget]::User)
$Env:VCPKG_ROOT = "$Env:GITHUB_WORKSPACE\vcpkg"
$Env:Path += ";$Env:VCPKG_ROOT"

vcpkg install hdf5:x64-windows-static