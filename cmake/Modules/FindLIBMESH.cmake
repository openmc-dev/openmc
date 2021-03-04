# Finds the libMesh installation using CMake's PkgConfig
# module and creates a libmesh imported target

if(${CMAKE_VERSION} VERSION_LESS 3.12.0)
    message(FATAL_ERROR "OpenMC builds with libMesh support require CMake version 3.12.0 or greater.")
endif()

set(LIBMESH_PC_FILE libmesh)

# if the METHOD variable is present, check specifically for
# the libMesh .pc file for that build type
if(DEFINED ENV{METHOD})
  set(LIBMESH_PC_FILE libmesh-$ENV{METHOD})
  message(STATUS "Using environment variable METHOD to determine libMesh build: ${LIBMESH_PC_FILE}")
endif()

include(FindPkgConfig)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH True)
pkg_check_modules(LIBMESH REQUIRED ${LIBMESH_PC_FILE}>=1.6.0 IMPORTED_TARGET)
