# Finds the libMesh installation using CMake's PkgConfig
# module and creates a libmesh imported target

find_path(LIBMESH_PC NAMES libmesh-opt.pc
          HINTS ${LIBMESH_DIR} $ENV{LIBMESH_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib/pkgconfig pkgconfig
          NO_DEFAULT_PATH)

include(FindPkgConfig)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
pkg_check_modules(LIBMESH REQUIRED libmesh IMPORTED_TARGET)
pkg_check_modules(LIBMESH_DEVEL REQUIRED libmesh-devel IMPORTED_TARGET)

message(STATUS "Found LIBMESH in ${LIBMESH_PC}")
