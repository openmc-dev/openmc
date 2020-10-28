# Finds the libMesh installation using CMake's PkgConfig
# module and creates a libmesh imported target

include(FindPkgConfig)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
pkg_check_modules(LIBMESH REQUIRED libmesh IMPORTED_TARGET)

find_path(LIBMESH_PC NAMES libmesh.pc
          HINTS ${LIBMESH_LIBDIR}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib/pkgconfig pkgconfig
          NO_DEFAULT_PATH)
