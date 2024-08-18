git clone https://github.com/microsoft/vcpkg.git "$Env:USERPROFILE\vcpkg"

cd "$Env:USERPROFILE\vcpkg"
.\bootstrap-vcpkg.bat

[Environment]::SetEnvironmentVariable("VCPKG_ROOT", "$Env:USERPROFILE\vcpkg", [System.EnvironmentVariableTarget]::User)
$Env:Path += ";$Env:VCPKG_ROOT"

vcpkg install hdf5:x64-windows-static