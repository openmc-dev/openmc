cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(openmc C CXX)

# Set version numbers
set(OPENMC_VERSION_MAJOR 0)
set(OPENMC_VERSION_MINOR 14)
set(OPENMC_VERSION_RELEASE 1)
set(OPENMC_VERSION ${OPENMC_VERSION_MAJOR}.${OPENMC_VERSION_MINOR}.${OPENMC_VERSION_RELEASE})
configure_file(include/openmc/version.h.in "${CMAKE_BINARY_DIR}/include/openmc/version.h" @ONLY)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set module path
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

# Allow user to specify <project>_ROOT variables
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

# Enable correct usage of CXX_EXTENSIONS
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
  cmake_policy(SET CMP0128 NEW)
endif()

#===============================================================================
# Command line options
#===============================================================================

option(OPENMC_USE_OPENMP      "Enable shared-memory parallelism with OpenMP"         ON)
option(OPENMC_BUILD_TESTS     "Build tests"                                          ON)
option(OPENMC_ENABLE_PROFILE  "Compile with profiling flags"                         OFF)
option(OPENMC_ENABLE_COVERAGE "Compile with coverage analysis flags"                 OFF)
option(OPENMC_USE_DAGMC       "Enable support for DAGMC (CAD) geometry"              OFF)
option(OPENMC_USE_LIBMESH     "Enable support for libMesh unstructured mesh tallies" OFF)
option(OPENMC_USE_MPI         "Enable MPI"                                           OFF)
option(OPENMC_USE_MCPL        "Enable MCPL"                                          OFF)
option(OPENMC_USE_NCRYSTAL    "Enable support for NCrystal scattering"               OFF)

# Warnings for deprecated options
foreach(OLD_OPT IN ITEMS "openmp" "profile" "coverage" "dagmc" "libmesh")
  if(DEFINED ${OLD_OPT})
    string(TOUPPER ${OLD_OPT} OPT_UPPER)
    if ("${OLD_OPT}" STREQUAL "profile" OR "${OLD_OPT}" STREQUAL "coverage")
      set(NEW_OPT_PREFIX "OPENMC_ENABLE")
    else()
      set(NEW_OPT_PREFIX "OPENMC_USE")
    endif()
    message(WARNING "The OpenMC CMake option '${OLD_OPT}' has been deprecated. "
      "Its value will be ignored. "
      "Please use '-D${NEW_OPT_PREFIX}_${OPT_UPPER}=${${OLD_OPT}}' instead.")
    unset(${OLD_OPT} CACHE)
  endif()
endforeach()

foreach(OLD_BLD in ITEMS "debug" "optimize")
  if(DEFINED ${OLD_BLD})
    if("${OLD_BLD}" STREQUAL "debug")
      set(BLD_VAR "Debug")
    else()
      set(BLD_VAR "Release")
    endif()
    message(WARNING "The OpenMC CMake option '${OLD_BLD}' has been deprecated. "
      "Its value will be ignored. "
      "OpenMC now uses the CMAKE_BUILD_TYPE variable to set the build mode. "
      "Please use '-DCMAKE_BUILD_TYPE=${BLD_VAR}' instead.")
    unset(${OLD_BLD} CACHE)
  endif()
endforeach()

#===============================================================================
# Set a default build configuration if not explicitly specified
#===============================================================================

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, defaulting to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
endif()

#===============================================================================
# MPI for distributed-memory parallelism
#===============================================================================

if(OPENMC_USE_MPI)
  find_package(MPI REQUIRED)
endif()

#===============================================================================
# Helper macro for finding a dependency
#===============================================================================

macro(find_package_write_status pkg)
  find_package(${pkg} QUIET NO_SYSTEM_ENVIRONMENT_PATH)
  if(${pkg}_FOUND)
    message(STATUS "Found ${pkg}: ${${pkg}_DIR} (version ${${pkg}_VERSION})")
  else()
    message(STATUS "Did not find ${pkg}, will use submodule instead")
  endif()
endmacro()

#===============================================================================
# NCrystal Scattering Support
#===============================================================================

if(OPENMC_USE_NCRYSTAL)
  find_package(NCrystal REQUIRED)
  message(STATUS "Found NCrystal: ${NCrystal_DIR} (version ${NCrystal_VERSION})")
endif()

#===============================================================================
# DAGMC Geometry Support - need DAGMC/MOAB
#===============================================================================

if(OPENMC_USE_DAGMC)
  find_package(DAGMC REQUIRED PATH_SUFFIXES lib/cmake)
  if (${DAGMC_VERSION} VERSION_LESS 3.2.0)
    message(FATAL_ERROR "Discovered DAGMC Version: ${DAGMC_VERSION}. \
    Please update DAGMC to version 3.2.0 or greater.")
  endif()
  message(STATUS "Found DAGMC: ${DAGMC_DIR} (version ${DAGMC_VERSION})")
endif()

#===============================================================================
# libMesh Unstructured Mesh Support
#===============================================================================

if(OPENMC_USE_LIBMESH)
  find_package(LIBMESH REQUIRED)
endif()

#===============================================================================
# libpng
#===============================================================================

find_package(PNG)

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

find_package(HDF5 REQUIRED COMPONENTS C HL)
if(HDF5_IS_PARALLEL)
  if(NOT OPENMC_USE_MPI)
    message(FATAL_ERROR "Parallel HDF5 was detected, but MPI was not enabled.\
      To use parallel HDF5, OpenMC needs to be built with MPI support by passing\
      -DOPENMC_USE_MPI=ON when calling cmake.")
  endif()
  message(STATUS "Using parallel HDF5")
endif()

# Version 1.12 of HDF5 deprecates the H5Oget_info_by_idx() interface.
# Thus, we give these flags to allow usage of the old interface in newer
# versions of HDF5.
if(${HDF5_VERSION} VERSION_GREATER_EQUAL 1.12.0)
  list(APPEND cxxflags -DH5Oget_info_by_idx_vers=1 -DH5O_info_t_vers=1)
endif()

#===============================================================================
# MCPL
#===============================================================================

if (OPENMC_USE_MCPL)
  find_package(MCPL REQUIRED)
  message(STATUS "Found MCPL: ${MCPL_DIR} (found version \"${MCPL_VERSION}\")")
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

# Skip for Visual Studio which has its own configurations through GUI
if(NOT MSVC)

if(OPENMC_USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    # In CMake 3.9+, can use the OpenMP::OpenMP_CXX imported target
    list(APPEND cxxflags ${OpenMP_CXX_FLAGS})
    list(APPEND ldflags ${OpenMP_CXX_FLAGS})
  endif()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(OPENMC_ENABLE_PROFILE)
  list(APPEND cxxflags -g -fno-omit-frame-pointer)
endif()

if(OPENMC_ENABLE_COVERAGE)
  list(APPEND cxxflags --coverage)
  list(APPEND ldflags --coverage)
endif()

# Show flags being used
message(STATUS "OpenMC C++ flags: ${cxxflags}")
message(STATUS "OpenMC Linker flags: ${ldflags}")

endif()

#===============================================================================
# Update git submodules as needed
#===============================================================================

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "git submodule update --init failed with \
        ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

# Check to see if submodules exist (by checking one)
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/vendor/pugixml/CMakeLists.txt")
  message(FATAL_ERROR "The git submodules were not downloaded! GIT_SUBMODULE was \
    turned off or failed. Please update submodules and try again.")
endif()

#===============================================================================
# pugixml library
#===============================================================================

find_package_write_status(pugixml)
if (NOT pugixml_FOUND)
  add_subdirectory(vendor/pugixml)
  set_target_properties(pugixml PROPERTIES CXX_STANDARD 14 CXX_EXTENSIONS OFF)
endif()

#===============================================================================
# {fmt} library
#===============================================================================

find_package_write_status(fmt)
if (NOT fmt_FOUND)
  set(FMT_INSTALL ON CACHE BOOL "Generate the install target.")
  add_subdirectory(vendor/fmt)
endif()

#===============================================================================
# xtensor header-only library
#===============================================================================

# CMake 3.13+ will complain about policy CMP0079 unless it is set explicitly
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
  cmake_policy(SET CMP0079 NEW)
endif()

find_package_write_status(xtensor)
if (NOT xtensor_FOUND)
  add_subdirectory(vendor/xtl)
  set(xtl_DIR ${CMAKE_CURRENT_BINARY_DIR}/vendor/xtl)
  add_subdirectory(vendor/xtensor)
endif()

#===============================================================================
# GSL header-only library
#===============================================================================

find_package_write_status(gsl-lite)
if (NOT gsl-lite_FOUND)
  add_subdirectory(vendor/gsl-lite)

  # Make sure contract violations throw exceptions
  target_compile_definitions(gsl-lite-v1 INTERFACE GSL_THROW_ON_CONTRACT_VIOLATION)
  target_compile_definitions(gsl-lite-v1 INTERFACE gsl_CONFIG_ALLOWS_NONSTRICT_SPAN_COMPARISON=1)
endif()

#===============================================================================
# Catch2 library
#===============================================================================

if(OPENMC_BUILD_TESTS)
  find_package_write_status(Catch2)
  if (NOT Catch2_FOUND)
    add_subdirectory(vendor/Catch2)
  endif()
endif()

#===============================================================================
# RPATH information
#===============================================================================

# Provide install directory variables as defined by GNU coding standards
include(GNUInstallDirs)

# This block of code ensures that dynamic libraries can be found via the RPATH
# whether the executable is the original one from the build directory or the
# installed one in CMAKE_INSTALL_PREFIX. Ref:
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
endif()

#===============================================================================
# libopenmc
#===============================================================================

list(APPEND libopenmc_SOURCES
  src/bank.cpp
  src/boundary_condition.cpp
  src/bremsstrahlung.cpp
  src/cell.cpp
  src/cmfd_solver.cpp
  src/cross_sections.cpp
  src/dagmc.cpp
  src/distribution.cpp
  src/distribution_angle.cpp
  src/distribution_energy.cpp
  src/distribution_multi.cpp
  src/distribution_spatial.cpp
  src/eigenvalue.cpp
  src/endf.cpp
  src/error.cpp
  src/event.cpp
  src/file_utils.cpp
  src/finalize.cpp
  src/geometry.cpp
  src/geometry_aux.cpp
  src/hdf5_interface.cpp
  src/initialize.cpp
  src/lattice.cpp
  src/material.cpp
  src/math_functions.cpp
  src/mcpl_interface.cpp
  src/mesh.cpp
  src/message_passing.cpp
  src/mgxs.cpp
  src/mgxs_interface.cpp
  src/ncrystal_interface.cpp
  src/nuclide.cpp
  src/output.cpp
  src/particle.cpp
  src/particle_data.cpp
  src/particle_restart.cpp
  src/photon.cpp
  src/physics.cpp
  src/physics_common.cpp
  src/physics_mg.cpp
  src/plot.cpp
  src/position.cpp
  src/progress_bar.cpp
  src/random_dist.cpp
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
  src/tallies/filter_cell.cpp
  src/tallies/filter_cell_instance.cpp
  src/tallies/filter_cellborn.cpp
  src/tallies/filter_cellfrom.cpp
  src/tallies/filter_collision.cpp
  src/tallies/filter_delayedgroup.cpp
  src/tallies/filter_distribcell.cpp
  src/tallies/filter_energy.cpp
  src/tallies/filter_energyfunc.cpp
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
  src/tallies/filter_time.cpp
  src/tallies/filter_universe.cpp
  src/tallies/filter_zernike.cpp
  src/tallies/tally.cpp
  src/tallies/tally_scoring.cpp
  src/tallies/trigger.cpp
  src/thermal.cpp
  src/timer.cpp
  src/track_output.cpp
  src/universe.cpp
  src/urr.cpp
  src/volume_calc.cpp
  src/weight_windows.cpp
  src/wmp.cpp
  src/xml_interface.cpp
  src/xsdata.cpp)

# Add bundled external dependencies
list(APPEND libopenmc_SOURCES
  src/external/quartic_solver.cpp
  src/external/Faddeeva.cc)

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

add_library(OpenMC::libopenmc ALIAS libopenmc)

# Avoid vs error lnk1149 :output filename matches input filename
if(NOT MSVC)
  set_target_properties(libopenmc PROPERTIES OUTPUT_NAME openmc)
endif()

target_include_directories(libopenmc
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    ${HDF5_INCLUDE_DIRS}
)

# Set compile flags
target_compile_options(libopenmc PRIVATE ${cxxflags})

# Add include directory for configured version file
target_include_directories(libopenmc
  PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>)

if (HDF5_IS_PARALLEL)
  target_compile_definitions(libopenmc PRIVATE -DPHDF5)
endif()
if (OPENMC_USE_MPI)
  target_compile_definitions(libopenmc PUBLIC -DOPENMC_MPI)
endif()

# Set git SHA1 hash as a compile definition
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  RESULT_VARIABLE GIT_SHA1_SUCCESS
                  OUTPUT_VARIABLE GIT_SHA1
                  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(GIT_SHA1_SUCCESS EQUAL 0)
    target_compile_definitions(libopenmc PRIVATE -DGIT_SHA1="${GIT_SHA1}")
  endif()
endif()

# target_link_libraries treats any arguments starting with - but not -l as
# linker flags. Thus, we can pass both linker flags and libraries together.
target_link_libraries(libopenmc ${ldflags} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
                      xtensor gsl::gsl-lite-v1 fmt::fmt)

if(TARGET pugixml::pugixml)
  target_link_libraries(libopenmc pugixml::pugixml)
else()
  target_link_libraries(libopenmc pugixml)
endif()

if(OPENMC_USE_DAGMC)
  target_compile_definitions(libopenmc PRIVATE DAGMC)
  target_link_libraries(libopenmc dagmc-shared uwuw-shared)
endif()

if(OPENMC_USE_LIBMESH)
  target_compile_definitions(libopenmc PRIVATE LIBMESH)
  target_link_libraries(libopenmc PkgConfig::LIBMESH)
endif()

if (PNG_FOUND)
  target_compile_definitions(libopenmc PRIVATE USE_LIBPNG)
  target_link_libraries(libopenmc PNG::PNG)
endif()

if (OPENMC_USE_MPI)
  target_link_libraries(libopenmc MPI::MPI_CXX)
endif()

if (OPENMC_BUILD_TESTS)
  # Add cpp tests directory
  include(CTest)
  add_subdirectory(tests/cpp_unit_tests)
endif()

if (OPENMC_USE_MCPL)
  target_compile_definitions(libopenmc PUBLIC OPENMC_MCPL)
  target_link_libraries(libopenmc MCPL::mcpl)
endif()

if(OPENMC_USE_NCRYSTAL)
  target_compile_definitions(libopenmc PRIVATE NCRYSTAL)
  target_link_libraries(libopenmc NCrystal::NCrystal)
endif()

#===============================================================================
# Log build info that this executable can report later
#===============================================================================
target_compile_definitions(libopenmc PRIVATE BUILD_TYPE=${CMAKE_BUILD_TYPE})
target_compile_definitions(libopenmc PRIVATE COMPILER_ID=${CMAKE_CXX_COMPILER_ID})
target_compile_definitions(libopenmc PRIVATE COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION})
if (OPENMC_ENABLE_PROFILE)
  target_compile_definitions(libopenmc PRIVATE PROFILINGBUILD)
endif()
if (OPENMC_ENABLE_COVERAGE)
  target_compile_definitions(libopenmc PRIVATE COVERAGEBUILD)
endif()

#===============================================================================
# openmc executable
#===============================================================================
add_executable(openmc src/main.cpp)
add_executable(OpenMC::openmc ALIAS openmc)
target_compile_options(openmc PRIVATE ${cxxflags})
target_include_directories(openmc PRIVATE ${CMAKE_BINARY_DIR}/include)
target_link_libraries(openmc libopenmc)

# Ensure C++14 standard is used and turn off GNU extensions
target_compile_features(openmc PUBLIC cxx_std_14)
target_compile_features(libopenmc PUBLIC cxx_std_14)
set_target_properties(openmc libopenmc PROPERTIES CXX_EXTENSIONS OFF)

#===============================================================================
# Python package
#===============================================================================

add_custom_command(TARGET libopenmc POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  $<TARGET_FILE:libopenmc>
  ${CMAKE_CURRENT_SOURCE_DIR}/openmc/lib/$<TARGET_FILE_NAME:libopenmc>
  COMMENT "Copying libopenmc to Python module directory")

#===============================================================================
# Install executable, scripts, manpage, license
#===============================================================================

configure_file(cmake/OpenMCConfig.cmake.in "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMCConfig.cmake" @ONLY)
configure_file(cmake/OpenMCConfigVersion.cmake.in "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMCConfigVersion.cmake" @ONLY)

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/OpenMC)
install(TARGETS openmc libopenmc
  EXPORT openmc-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(EXPORT openmc-targets
  FILE OpenMCTargets.cmake
  NAMESPACE OpenMC::
  DESTINATION ${INSTALL_CONFIGDIR})

install(FILES
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMCConfig.cmake"
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/OpenMCConfigVersion.cmake"
  DESTINATION ${INSTALL_CONFIGDIR})
install(FILES man/man1/openmc.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
install(FILES LICENSE DESTINATION "${CMAKE_INSTALL_DOCDIR}" RENAME copyright)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${CMAKE_BINARY_DIR}/include/openmc/version.h" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/openmc)
