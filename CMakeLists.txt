cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(openmc Fortran C CXX)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Make sure Fortran module directory is included when building
include_directories(${CMAKE_BINARY_DIR}/include)

#===============================================================================
# Architecture specific definitions
#===============================================================================

if (${UNIX})
  add_definitions(-DUNIX)
endif()

#===============================================================================
# Command line options
#===============================================================================

option(openmp   "Enable shared-memory parallelism with OpenMP"   ON)
option(profile  "Compile with profiling flags"                   OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)
option(coverage "Compile with coverage analysis flags"           OFF)
option(mpif08   "Use Fortran 2008 MPI interface"                 OFF)

# Maximum number of nested coordinates levels
set(maxcoord 10 CACHE STRING "Maximum number of nested coordinate levels")
add_definitions(-DMAX_COORD=${maxcoord})

#===============================================================================
# MPI for distributed-memory parallelism
#===============================================================================

set(MPI_ENABLED FALSE)
if($ENV{FC} MATCHES "(mpi[^/]*|ftn)$")
  message("-- Detected MPI wrapper: $ENV{FC}")
  add_definitions(-DOPENMC_MPI)
  set(MPI_ENABLED TRUE)

  # Get directory containing MPI wrapper
  get_filename_component(MPI_DIR $ENV{FC} DIRECTORY)
endif()

# Check for Fortran 2008 MPI interface
if(MPI_ENABLED AND mpif08)
  message("-- Using Fortran 2008 MPI bindings")
  add_definitions(-DOPENMC_MPIF08)
endif()

#===============================================================================
# HDF5 for binary output
#===============================================================================

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

find_package(HDF5 COMPONENTS Fortran_HL)
if(NOT HDF5_FOUND)
  message(FATAL_ERROR "Could not find HDF5")
endif()
if(HDF5_IS_PARALLEL)
  if(NOT MPI_ENABLED)
    message(FATAL_ERROR "Parallel HDF5 must be used with MPI.")
  endif()
  add_definitions(-DPHDF5)
  message("-- Using parallel HDF5")
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

# Support for Fortran in FindOpenMP was added in CMake 3.1. To support lower
# versions, we manually add the flags. However, at some point in time, the
# manual logic can be removed in favor of the block below

#if(NOT (CMAKE_VERSION VERSION_LESS 3.1))
#  if(openmp)
#    find_package(OpenMP)
#    if(OPENMP_FOUND)
#      list(APPEND f90flags ${OpenMP_Fortran_FLAGS})
#      list(APPEND ldflags ${OpenMP_Fortran_FLAGS})
#    endif()
#  endif()
#endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Make sure version is sufficient
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)
  if(GCC_VERSION VERSION_LESS 4.8)
    message(FATAL_ERROR "gfortran version must be 4.8 or higher")
  endif()

  # GCC compiler options
  list(APPEND f90flags -cpp -std=f2008ts -fbacktrace -O2 -fstack-arrays)
  if(debug)
    list(REMOVE_ITEM f90flags -O2 -fstack-arrays)
    list(APPEND f90flags -g -Wall -Wno-unused-dummy-argument -pedantic
      -fbounds-check -ffpe-trap=invalid,overflow,underflow)
    list(APPEND ldflags -g)
  endif()
  if(profile)
    list(APPEND f90flags -pg)
    list(APPEND ldflags -pg)
  endif()
  if(optimize)
    list(REMOVE_ITEM f90flags -O2)
    list(APPEND f90flags -O3)
  endif()
  if(openmp)
    list(APPEND f90flags -fopenmp)
    list(APPEND ldflags -fopenmp)
  endif()
  if(coverage)
    list(APPEND f90flags -coverage)
    list(APPEND ldflags -coverage)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  # Intel compiler options
  list(APPEND f90flags -fpp -std08 -assume byterecl -traceback)
  if(debug)
    list(APPEND f90flags -g -warn -ftrapuv -fp-stack-check
      "-check all" -fpe0 -O0)
    list(APPEND ldflags -g)
  endif()
  if(profile)
    list(APPEND f90flags -pg)
    list(APPEND ldflags -pg)
  endif()
  if(optimize)
    list(APPEND f90flags -O3)
  endif()
  if(openmp)
    list(APPEND f90flags -qopenmp)
    list(APPEND ldflags -qopenmp)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
  # PGI Fortran compiler options
  list(APPEND f90flags -Mpreprocess -Minform=inform -traceback)
  add_definitions(-DNO_F2008)
  if(debug)
    list(APPEND f90flags -g -Mbounds -Mchkptr -Mchkstk)
    list(APPEND ldflags -g)
  endif()
  if(profile)
    list(APPEND f90flags -pg)
    list(APPEND ldflags -pg)
  endif()
  if(optimize)
    list(APPEND f90flags -fast -Mipa)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL XL)
  # IBM XL compiler options
  list(APPEND f90flags -O2)
  add_definitions(-DNO_F2008)
  if(debug)
    list(REMOVE_ITEM f90flags -O2)
    list(APPEND f90flags -g -C -qflag=i:i -u -O0)
    list(APPEND ldflags -g)
  endif()
  if(profile)
    list(APPEND f90flags -p)
    list(APPEND ldflags -p)
  endif()
  if(optimize)
    list(REMOVE_ITEM f90flags -O2)
    list(APPEND f90flags -O3)
  endif()
  if(openmp)
    list(APPEND f90flags -qsmp=omp)
    list(APPEND ldflags -qsmp=omp)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  # Cray Fortran compiler options
  list(APPEND f90flags -e Z -m 0)
  if(debug)
    list(APPEND f90flags -g -R abcnsp -O0)
    list(APPEND ldflags -g)
  endif()

endif()

if(CMAKE_C_COMPILER_ID STREQUAL GNU)
  # GCC compiler options
  list(APPEND cflags -std=c99 -O2)
  if(debug)
    list(REMOVE_ITEM cflags -O2)
    list(APPEND cflags -g -Wall -pedantic -fbounds-check)
  endif()
  if(profile)
    list(APPEND cflags -pg)
  endif()
  if(optimize)
    list(REMOVE_ITEM cflags -O2)
    list(APPEND cflags -O3)
  endif()
  if(coverage)
    list(APPEND cflags -coverage)
  endif()

elseif(CMAKE_C_COMPILER_ID STREQUAL Intel)
  # Intel compiler options
  list(APPEND cflags -std=c99)
  if(debug)
    list(APPEND cflags -g -w3 -ftrapuv -fp-stack-check -O0)
  endif()
  if(profile)
    list(APPEND cflags -pg)
  endif()
  if(optimize)
    list(APPEND cflags -O3)
  endif()

elseif(CMAKE_C_COMPILER_ID MATCHES Clang)
  # Clang options
  list(APPEND cflags -std=c99)
  if(debug)
    list(APPEND cflags -g -O0 -ftrapv)
  endif()
  if(optimize)
    list(APPEND cflags -O3)
  endif()

endif()

list(APPEND cxxflags -std=c++11 -O2)
if(debug)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -g -O0)
endif()
if(profile)
  list(APPEND cxxflags -pg)
endif()
if(optimize)
  list(REMOVE_ITEM cxxflags -O2)
  list(APPEND cxxflags -O3)
endif()
if(openmp)
  list(APPEND cxxflags -fopenmp)
endif()

# Show flags being used
message(STATUS "Fortran flags: ${f90flags}")
message(STATUS "C flags: ${cflags}")
message(STATUS "C++ flags: ${cxxflags}")
message(STATUS "Linker flags: ${ldflags}")

#===============================================================================
# git SHA1 hash
#===============================================================================

execute_process(COMMAND git rev-parse HEAD
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SHA1_SUCCESS
                OUTPUT_VARIABLE GIT_SHA1
                ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
if(GIT_SHA1_SUCCESS EQUAL 0)
  add_definitions(-DGIT_SHA1="${GIT_SHA1}")
endif()

#===============================================================================
# pugixml library
#===============================================================================

add_library(pugixml src/pugixml/pugixml_c.cpp src/pugixml/pugixml.cpp)
add_library(pugixml_fortran src/pugixml/pugixml_f.F90)
target_link_libraries(pugixml_fortran pugixml)

#===============================================================================
# RPATH information
#===============================================================================

# This block of code ensures that dynamic libraries can be found via the RPATH
# whether the executable is the original one from the build directory or the
# installed one in CMAKE_INSTALL_PREFIX. Ref:
# https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH

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
# Build faddeeva library
#===============================================================================

add_library(faddeeva STATIC src/faddeeva/Faddeeva.c)

#===============================================================================
# List source files.  Define the libopenmc and the OpenMC executable
#===============================================================================

set(program "openmc")
set(LIBOPENMC_FORTRAN_SRC
  src/algorithm.F90
  src/angle_distribution.F90
  src/angleenergy_header.F90
  src/bank_header.F90
  src/api.F90
  src/cmfd_data.F90
  src/cmfd_execute.F90
  src/cmfd_header.F90
  src/cmfd_input.F90
  src/cmfd_loss_operator.F90
  src/cmfd_prod_operator.F90
  src/cmfd_solver.F90
  src/constants.F90
  src/dict_header.F90
  src/distribution_multivariate.F90
  src/distribution_univariate.F90
  src/doppler.F90
  src/eigenvalue.F90
  src/endf.F90
  src/endf_header.F90
  src/energy_distribution.F90
  src/error.F90
  src/geometry.F90
  src/geometry_header.F90
  src/hdf5_interface.F90
  src/initialize.F90
  src/input_xml.F90
  src/list_header.F90
  src/material_header.F90
  src/math.F90
  src/matrix_header.F90
  src/mesh.F90
  src/mesh_header.F90
  src/message_passing.F90
  src/mgxs_data.F90
  src/mgxs_header.F90
  src/multipole.F90
  src/multipole_header.F90
  src/nuclide_header.F90
  src/output.F90
  src/particle_header.F90
  src/particle_restart.F90
  src/physics_common.F90
  src/physics.F90
  src/physics_mg.F90
  src/plot.F90
  src/plot_header.F90
  src/product_header.F90
  src/progress_header.F90
  src/random_lcg.F90
  src/reaction_header.F90
  src/relaxng
  src/sab_header.F90
  src/scattdata_header.F90
  src/secondary_correlated.F90
  src/secondary_kalbach.F90
  src/secondary_nbody.F90
  src/secondary_uncorrelated.F90
  src/set_header.F90
  src/settings.F90
  src/simulation_header.F90
  src/simulation.F90
  src/source.F90
  src/source_header.F90
  src/state_point.F90
  src/stl_vector.F90
  src/string.F90
  src/summary.F90
  src/surface_header.F90
  src/timer_header.F90
  src/tracking.F90
  src/track_output.F90
  src/urr_header.F90
  src/vector_header.F90
  src/volume_calc.F90
  src/volume_header.F90
  src/xml_interface.F90
  src/tallies/tally.F90
  src/tallies/tally_derivative_header.F90
  src/tallies/tally_filter.F90
  src/tallies/tally_filter_header.F90
  src/tallies/tally_filter_azimuthal.F90
  src/tallies/tally_filter_cell.F90
  src/tallies/tally_filter_cellborn.F90
  src/tallies/tally_filter_cellfrom.F90
  src/tallies/tally_filter_delayedgroup.F90
  src/tallies/tally_filter_distribcell.F90
  src/tallies/tally_filter_energy.F90
  src/tallies/tally_filter_energyfunc.F90
  src/tallies/tally_filter_material.F90
  src/tallies/tally_filter_mesh.F90
  src/tallies/tally_filter_mu.F90
  src/tallies/tally_filter_polar.F90
  src/tallies/tally_filter_surface.F90
  src/tallies/tally_filter_universe.F90
  src/tallies/tally_header.F90
  src/tallies/trigger.F90
  src/tallies/trigger_header.F90
)
set(LIBOPENMC_CXX_SRC
  src/error.h
  src/hdf5_interface.h
  src/random_lcg.cpp
  src/random_lcg.h
  src/surface.cpp
  src/surface.h
  src/xml_interface.h
  src/pugixml/pugixml.cpp
  src/pugixml/pugixml.hpp)
add_library(libopenmc SHARED ${LIBOPENMC_FORTRAN_SRC} ${LIBOPENMC_CXX_SRC})
set_target_properties(libopenmc PROPERTIES OUTPUT_NAME openmc)
add_executable(${program} src/main.F90)

#===============================================================================
# Add compiler/linker flags
#===============================================================================

set_property(TARGET ${program} libopenmc pugixml_fortran
  PROPERTY LINKER_LANGUAGE Fortran)

target_include_directories(libopenmc PUBLIC ${HDF5_INCLUDE_DIRS})

# The executable and the faddeeva package use only one language.  They can be
# set via target_compile_options which accepts a list.
target_compile_options(${program} PUBLIC ${f90flags})
target_compile_options(faddeeva PRIVATE ${cflags})

# The libopenmc library has both F90 and C++ so the compile flags must be set
# file-by-file via set_source_file_properties.  The compile flags must first be
# converted from lists to strings.
string(REPLACE ";" " " f90flags "${f90flags}")
string(REPLACE ";" " " cxxflags "${cxxflags}")
set_source_files_properties(${LIBOPENMC_FORTRAN_SRC} PROPERTIES COMPILE_FLAGS
                            ${f90flags})
set_source_files_properties(${LIBOPENMC_CXX_SRC} PROPERTIES COMPILE_FLAGS
                            ${cxxflags})

# Add HDF5 library directories to link line with -L
foreach(LIBDIR ${HDF5_LIBRARY_DIRS})
  list(APPEND ldflags "-L${LIBDIR}")
endforeach()

# target_link_libraries treats any arguments starting with - but not -l as
# linker flags. Thus, we can pass both linker flags and libraries together.
target_link_libraries(libopenmc ${ldflags} ${HDF5_LIBRARIES} pugixml_fortran
                      faddeeva)
target_link_libraries(${program} ${ldflags} libopenmc)

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

install(TARGETS ${program} libopenmc
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib)
install(DIRECTORY src/relaxng DESTINATION share/openmc)
install(FILES man/man1/openmc.1 DESTINATION share/man/man1)
install(FILES LICENSE DESTINATION "share/doc/${program}" RENAME copyright)
