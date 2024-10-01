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

find_package(PkgConfig REQUIRED)

# Set hints to search for libmesh.pc in CMAKE_PREFIX_PATH
find_path(PKGCONFIG_DIR
    NAMES ${LIBMESH_PC_FILE}.pc
    HINTS ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES lib/pkgconfig lib64/pkgconfig lib/x86_64-linux-gnu/pkgconfig share/pkgconfig
    NO_DEFAULT_PATH
    DOC "The path to libMesh pkgconfig directory"
)

if(PKGCONFIG_DIR)
    message(STATUS "Found ${LIBMESH_PC_FILE}.pc in: ${PKGCONFIG_DIR}")
    set(ENV{PKG_CONFIG_PATH} "${PKGCONFIG_DIR}:$ENV{PKG_CONFIG_PATH}")
else()
    message(FATAL_ERROR "Could not find ${LIBMESH_PC_FILE}.pc in CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
endif()

set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)
pkg_check_modules(LIBMESH REQUIRED ${LIBMESH_PC_FILE}>=1.7.0 IMPORTED_TARGET)
pkg_get_variable(LIBMESH_PREFIX ${LIBMESH_PC_FILE} prefix)
