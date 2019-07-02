cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(openmc C CXX)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

#===============================================================================
# Command line options
#===============================================================================

option(openmp   "Enable shared-memory parallelism with OpenMP"   ON)
option(profile  "Compile with profiling flags"                   OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)
option(coverage "Compile with coverage analysis flags"           OFF)
option(dagmc    "Enable support for DAGMC (CAD) geometry"        OFF)

#===============================================================================
# MPI for distributed-memory parallelism
#===============================================================================

set(MPI_ENABLED FALSE)
if($ENV{CXX} MATCHES "(mpi[^/]*|CC)$")
  message(STATUS "Detected MPI wrapper: $ENV{CXX}")
  set(MPI_ENABLED TRUE)
endif()

#===============================================================================
# DAGMC Geometry Support - need DAGMC/MOAB
#===============================================================================
if(dagmc)
  find_package(DAGMC REQUIRED)
  link_directories(${DAGMC_LIBRARY_DIRS})
endif()

#===============================================================================
# HDF5 for binary output
#===============================================================================

# Allow user to specify HDF5_ROOT
if (NOT (CMAKE_VERSION VERSION_LESS 3.12))
  cmake_policy(SET CMP0074 NEW)
endif()

# Unfortunately FindHDF5.cmake will always prefer a serial HDF5 installation
# over a parallel installation if both appear on the user's PATH. To get around
# this, we check for the environment variable HDF5_ROOT and if it exists, use it
# to check whether its a parallel version.

if(NOT DEFINED HDF5_PREFER_PARALLEL)
  if(DEFINED ENV{HDF5_ROOT} AND EXISTS $ENV{HDF5_ROOT}/bin/h5pcc)
    set(HDF5_PREFER_PARALLEL TRUE)
  else()
    set(HDF5_PREFER_PARALLEL FALSE)
  endif()
endif()

find_package(HDF5 REQUIRED COMPONENTS C HL)
if(HDF5_IS_PARALLEL)
  if(NOT MPI_ENABLED)
    message(FATAL_ERROR "Parallel HDF5 must be used with MPI.")
  endif()
  message(STATUS "Using parallel HDF5")
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

# Skip for Visual Stduio which has its own configurations through GUI
if(NOT MSVC)

if(openmp)
  # Requires CMake 3.1+
  find_package(OpenMP)
  if(OPENMP_FOUND)
    list(APPEND cxxflags ${OpenMP_CXX_FLAGS})
    list(APPEND ldflags ${OpenMP_CXX_FLAGS})
  endif()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

list(APPEND cxxflags -O2)
if(debug)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -g -O0)
endif()
if(profile)
  list(APPEND cxxflags -g -fno-omit-frame-pointer)
endif()
if(optimize)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -O3)
endif()
if(coverage)
  list(APPEND cxxflags --coverage)
  list(APPEND ldflags --coverage)
endif()

# Show flags being used
message(STATUS "OpenMC C++ flags: ${cxxflags}")
message(STATUS "OpenMC Linker flags: ${ldflags}")

endif()

#===============================================================================
# pugixml library
#===============================================================================

add_library(pugixml vendor/pugixml/pugixml.cpp)
target_include_directories(pugixml PUBLIC vendor/pugixml/)

#===============================================================================
# xtensor header-only library
#===============================================================================

add_subdirectory(vendor/xtl)
add_subdirectory(vendor/xtensor)
target_link_libraries(xtensor INTERFACE xtl)

#===============================================================================
# GSL header-only library
#===============================================================================

add_library(gsl INTERFACE)
target_include_directories(gsl INTERFACE vendor/gsl/include)

# Make sure contract violations throw exceptions
target_compile_definitions(gsl INTERFACE GSL_THROW_ON_CONTRACT_VIOLATION)

#===============================================================================
# RPATH information
#===============================================================================

# This block of code ensures that dynamic libraries can be found via the RPATH
# whether the executable is the original one from the build directory or the
# installed one in CMAKE_INSTALL_PREFIX. Ref:
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#===============================================================================
# faddeeva library
#===============================================================================

add_library(faddeeva STATIC vendor/faddeeva/Faddeeva.cc)
target_include_directories(faddeeva PUBLIC vendor/faddeeva/)
target_compile_options(faddeeva PRIVATE ${cxxflags})

#===============================================================================
# libopenmc
#===============================================================================

list(APPEND libopenmc_SOURCES
  src/bank.cpp
  src/bremsstrahlung.cpp
  src/dagmc.cpp
  src/cell.cpp
  src/cmfd_solver.cpp
  src/cross_sections.cpp
  src/distribution.cpp
  src/distribution_angle.cpp
  src/distribution_energy.cpp
  src/distribution_multi.cpp
  src/distribution_spatial.cpp
  src/eigenvalue.cpp
  src/endf.cpp
  src/error.cpp
  src/initialize.cpp
  src/finalize.cpp
  src/geometry.cpp
  src/geometry_aux.cpp
  src/hdf5_interface.cpp
  src/lattice.cpp
  src/material.cpp
  src/math_functions.cpp
  src/mesh.cpp
  src/message_passing.cpp
  src/mgxs.cpp
  src/mgxs_interface.cpp
  src/nuclide.cpp
  src/output.cpp
  src/particle.cpp
  src/particle_restart.cpp
  src/photon.cpp
  src/physics.cpp
  src/physics_common.cpp
  src/physics_mg.cpp
  src/plot.cpp
  src/position.cpp
  src/progress_bar.cpp
  src/random_lcg.cpp
  src/reaction.cpp
  src/reaction_product.cpp
  src/scattdata.cpp
  src/secondary_correlated.cpp
  src/secondary_kalbach.cpp
  src/secondary_nbody.cpp
  src/secondary_thermal.cpp
  src/secondary_uncorrelated.cpp
  src/settings.cpp
  src/simulation.cpp
  src/source.cpp
  src/state_point.cpp
  src/string_utils.cpp
  src/summary.cpp
  src/surface.cpp
  src/tallies/derivative.cpp
  src/tallies/filter.cpp
  src/tallies/filter_azimuthal.cpp
  src/tallies/filter_cellborn.cpp
  src/tallies/filter_cellfrom.cpp
  src/tallies/filter_cell.cpp
  src/tallies/filter_delayedgroup.cpp
  src/tallies/filter_distribcell.cpp
  src/tallies/filter_energyfunc.cpp
  src/tallies/filter_energy.cpp
  src/tallies/filter_legendre.cpp
  src/tallies/filter_material.cpp
  src/tallies/filter_mesh.cpp
  src/tallies/filter_meshsurface.cpp
  src/tallies/filter_mu.cpp
  src/tallies/filter_particle.cpp
  src/tallies/filter_polar.cpp
  src/tallies/filter_sph_harm.cpp
  src/tallies/filter_sptl_legendre.cpp
  src/tallies/filter_surface.cpp
  src/tallies/filter_universe.cpp
  src/tallies/filter_zernike.cpp
  src/tallies/tally.cpp
  src/tallies/tally_scoring.cpp
  src/tallies/trigger.cpp
  src/timer.cpp
  src/thermal.cpp
  src/track_output.cpp
  src/urr.cpp
  src/volume_calc.cpp
  src/wmp.cpp
  src/xml_interface.cpp
  src/xsdata.cpp)

# For Visual Studio compilers
if(MSVC)
  # Use static library (otherwise explicit symbol portings are needed)
  add_library(libopenmc STATIC ${libopenmc_SOURCES})

  # To use the shared HDF5 libraries on Windows, the H5_BUILT_AS_DYNAMIC_LIB
  # compile definition must be specified.
  target_compile_definitions(libopenmc PRIVATE -DH5_BUILT_AS_DYNAMIC_LIB)
else()
  add_library(libopenmc SHARED ${libopenmc_SOURCES})
endif()

set_target_properties(libopenmc PROPERTIES
  OUTPUT_NAME openmc)

target_include_directories(libopenmc
  PUBLIC include ${HDF5_INCLUDE_DIRS})

# Set compile flags
target_compile_options(libopenmc PRIVATE ${cxxflags})

if (HDF5_IS_PARALLEL)
  target_compile_definitions(libopenmc PRIVATE -DPHDF5)
endif()
if (MPI_ENABLED)
  target_compile_definitions(libopenmc PUBLIC -DOPENMC_MPI)
endif()

# Set git SHA1 hash as a compile definition
execute_process(COMMAND git rev-parse HEAD
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SHA1_SUCCESS
                OUTPUT_VARIABLE GIT_SHA1
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if(GIT_SHA1_SUCCESS EQUAL 0)
  target_compile_definitions(libopenmc PRIVATE -DGIT_SHA1="${GIT_SHA1}")
endif()

# target_link_libraries treats any arguments starting with - but not -l as
# linker flags. Thus, we can pass both linker flags and libraries together.
target_link_libraries(libopenmc ${ldflags} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
                      pugixml faddeeva xtensor gsl)

if(dagmc)
  target_compile_definitions(libopenmc PRIVATE DAGMC)
  target_link_libraries(libopenmc ${DAGMC_LIBRARIES})
  target_include_directories(libopenmc PRIVATE ${DAGMC_INCLUDE_DIRS})
endif()

#===============================================================================
# openmc executable
#===============================================================================
add_executable(openmc src/main.cpp)
target_compile_options(openmc PRIVATE ${cxxflags})
target_link_libraries(openmc libopenmc)

# Ensure C++14 standard is used. Starting with CMake 3.8, another way this could
# be done is using the cxx_std_14 compiler feature.
set_target_properties(
    openmc libopenmc faddeeva pugixml
    PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)

#===============================================================================
# Python package
#===============================================================================

add_custom_command(TARGET libopenmc POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  $<TARGET_FILE:libopenmc>
  ${CMAKE_CURRENT_SOURCE_DIR}/openmc/capi/$<TARGET_FILE_NAME:libopenmc>
  COMMENT "Copying libopenmc to Python module directory")

#===============================================================================
# Install executable, scripts, manpage, license
#===============================================================================

install(TARGETS openmc libopenmc
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  )
install(DIRECTORY src/relaxng DESTINATION share/openmc)
install(FILES man/man1/openmc.1 DESTINATION share/man/man1)
install(FILES LICENSE DESTINATION "share/doc/openmc" RENAME copyright)
install(DIRECTORY include/ DESTINATION include)
