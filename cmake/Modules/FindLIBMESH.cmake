# Try to find LIBMESH
#
# Once done this will define
#
#  LIBMESH_FOUND - system has LIBMESH
#  LIBMESH_INCLUDE_DIRS - the LIBMESH include directory
#  LIBMESH_LIBRARIES - Link these to use LIBMESH
#  LIBMESH_DEFINITIONS - Compiler switches required for using LIBMESH

find_path(LIBMESH_PC NAMES libmesh-opt.pc
          HINTS ${LIBMESH_ROOT} $ENV{LIBMESH_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib/pkgconfig pkgconfig
          NO_DEFAULT_PATH)

include(FindPkgConfig)
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBMESH_PC}")
pkg_check_modules(LIBMESH REQUIRED libmesh IMPORTED_TARGET)

message(STATUS "Found LIBMESH in ${LIBMESH_PC}")
