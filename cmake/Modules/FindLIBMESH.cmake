# Finds the libMesh installation using CMake's PkgConfig
# module and creates a libmesh imported target

if (${CMAKE_VERSION} VERSION_LESS 3.12.0)
    message(FATAL_ERROR "OpenMC builds with libMesh support require CMake version 3.12.0 or greater.")
endif()

include(FindPkgConfig)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH True)
pkg_check_modules(LIBMESH REQUIRED libmesh>=1.6.0 IMPORTED_TARGET)

find_path(LIBMESH_PC NAMES libmesh.pc
          HINTS ${LIBMESH_LIBDIR}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib/pkgconfig pkgconfig
          NO_DEFAULT_PATH)
