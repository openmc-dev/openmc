# Try to find LIBMESH
#
# Once done this will define
#
#  LIBMESH_FOUND - system has LIBMESH
#  LIBMESH_INCLUDE_DIRS - the LIBMESH include directory
#  LIBMESH_LIBRARIES - Link these to use LIBMESH
#  LIBMESH_DEFINITIONS - Compiler switches required for using LIBMESH

find_path(LIBMESH NAMES libmesh_opt.so
          HINTS ${LIBMESH_ROOT} $ENV{LIBMESH_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib
          NO_DEFAULT_PATH)

if(DEFINED LIBMESH)
  set(LIBMESH_FOUND TRUE)
endif()


message(STATUS "Found LIBMESH in ${LIBMESH}")
