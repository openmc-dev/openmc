# Try to find DAGMC
#
# Once done this will define
#
#  DAGMC_FOUND - system has DAGMC
#  DAGMC_INCLUDE_DIRS - the DAGMC include directory
#  DAGMC_LIBRARIES - Link these to use DAGMC
#  DAGMC_DEFINITIONS - Compiler switches required for using DAGMC

find_path(DAGMC_CMAKE_CONFIG NAMES DAGMCConfig.cmake
          HINTS ${DAGMC_ROOT} $ENV{DAGMC_ROOT}
          PATHS ENV LD_LIBRARY_PATH
          PATH_SUFFIXES lib Lib cmake lib/cmake
          NO_DEFAULT_PATH)
if(DAGMC_CMAKE_CONFIG)
  message(STATUS "Found DAGMC in ${DAGMC_CMAKE_CONFIG}")
  include(${DAGMC_CMAKE_CONFIG}/DAGMCConfig.cmake)
else()
  message(WARNING "Cound not find DAGMC")
endif()
