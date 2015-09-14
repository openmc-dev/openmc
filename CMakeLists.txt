cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(openmc Fortran)

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

option(openmp   "Enable shared-memory parallelism with OpenMP"   OFF)
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
if($ENV{FC} MATCHES "mpi[^/]*$")
  message("-- Detected MPI wrapper: $ENV{FC}")
  add_definitions(-DMPI)
  set(MPI_ENABLED TRUE)
endif()

# Check for Fortran 2008 MPI interface
if(MPI_ENABLED AND mpif08)
  message("-- Using Fortran 2008 MPI bindings")
  add_definitions(-DMPIF08)
endif()

#===============================================================================
# HDF5 for binary output
#===============================================================================

# Unfortunately FindHDF5.cmake will always prefer a serial HDF5 installation
# over a parallel installation if both appear on the user's PATH. To get around
# this, we check for the environment variable HDF5_ROOT and if it exists, use it
# to check whether its a parallel version.

if(DEFINED ENV{HDF5_ROOT} AND EXISTS $ENV{HDF5_ROOT}/bin/h5pcc)
  set(HDF5_PREFER_PARALLEL TRUE)
else()
  set(HDF5_PREFER_PARALLEL FALSE)
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

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  # Make sure version is sufficient
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)
  if(GCC_VERSION VERSION_LESS 4.6)
    message(FATAL_ERROR "gfortran version must be 4.6 or higher")
  endif()

  # GNU Fortran compiler options
  list(APPEND f90flags -cpp -std=f2008 -fbacktrace)
  if(debug)
    if(NOT (GCC_VERSION VERSION_LESS 4.7))
      list(APPEND f90flags -Wall)
    endif()
    list(APPEND f90flags -g -pedantic -fbounds-check
      -ffpe-trap=invalid,overflow,underflow)
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
    list(APPEND f90flags -fopenmp)
    list(APPEND ldflags -fopenmp)
  endif()
  if(coverage)
    list(APPEND f90flags -coverage)
    list(APPEND ldflags -coverage)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  # Intel Fortran compiler options
  list(APPEND f90flags -fpp -std08 -assume byterecl -traceback)
  if(debug)
    list(APPEND f90flags -g -warn -ftrapuv -fp-stack-check
      "-check all" -fpe0)
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
    list(APPEND f90flags -openmp)
    list(APPEND ldflags -openmp)
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
    list(APPEND f90flags -g -C -qflag=i:i -u)
    list(APPEND ldflags -g)
  endif()
  if(profile)
    list(APPEND f90flags -p)
    list(APPEND ldflags -p)
  endif()
  if(optimize)
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
# FoX Fortran XML Library
#===============================================================================

# Only initialize git submodules if it is not there. User is responsible
# for future updates of fox xml submodule.
if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/src/xml/fox/.git)
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    message("-- Cloning FoX XML git repository...")
    execute_process(COMMAND git clone https://github.com/mit-crpg/fox.git src/xml/fox
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    execute_process(COMMAND git checkout bdc852f4f43d969fb1b179cba79295c1e095a455
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/xml/fox)
  else()
    message("-- Initializing/Updating FoX XML submodule...")
    execute_process(COMMAND git submodule init
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    execute_process(COMMAND git submodule update
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
endif()
add_subdirectory(src/xml/fox)

#===============================================================================
# Build OpenMC executable
#===============================================================================

set(program "openmc")
file(GLOB source src/*.F90 src/xml/openmc_fox.F90)
add_executable(${program} ${source})

# target_include_directories was added in CMake 2.8.11 and is the recommended
# way to set include directories. For lesser versions, we revert to set_property
if(CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories(${HDF5_INCLUDE_DIRS})
else()
  target_include_directories(${program} PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

# target_compile_options was added in CMake 2.8.12 and is the recommended way to
# set compile flags. Note that this sets the COMPILE_OPTIONS property (also
# available only in 2.8.12+) rather than the COMPILE_FLAGS property, which is
# deprecated. The former can handle lists whereas the latter cannot.
if(CMAKE_VERSION VERSION_LESS 4.8.12)
  string(REPLACE ";" " " f90flags "${f90flags}")
  set_property(TARGET ${program} PROPERTY COMPILE_FLAGS "${f90flags}")
else()
  target_compile_options(${program} PUBLIC ${f90flags})
endif()

# Add HDF5 library directories to link line with -L
foreach(LIBDIR ${HDF5_LIBRARY_DIRS})
  list(APPEND ldflags "-L${LIBDIR}")
endforeach()

# target_link_libraries treats any arguments starting with - but not -l as
# linker flags. Thus, we can pass both linker flags and libraries together.
target_link_libraries(${program} ${ldflags} ${HDF5_LIBRARIES} fox_dom)

#===============================================================================
# Install executable, scripts, manpage, license
#===============================================================================

install(TARGETS ${program} RUNTIME DESTINATION bin)
install(DIRECTORY src/relaxng DESTINATION share/openmc)
install(FILES man/man1/openmc.1 DESTINATION share/man/man1)
install(FILES LICENSE DESTINATION "share/doc/${program}/copyright")

find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
  install(CODE "execute_process(
                  COMMAND ${PYTHON_EXECUTABLE} setup.py install
                  --prefix=${CMAKE_INSTALL_PREFIX}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})")
endif()

#===============================================================================
# Regression tests
#===============================================================================

# This allows for dashboard configuration
include(CTest)

# Get a list of all the tests to run
file(GLOB_RECURSE TESTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.py)

# Check for MEM_CHECK and COVERAGE variables
if (DEFINED ENV{MEM_CHECK})
  set(MEM_CHECK $ENV{MEM_CHECK})
else(DEFINED ENV{MEM_CHECK})
  set(MEM_CHECK FALSE)
endif(DEFINED ENV{MEM_CHECK})
if (DEFINED ENV{COVERAGE})
  set(COVERAGE $ENV{COVERAGE})
else(DEFINED ENV{COVERAGE})
  set(COVERAGE FALSE)
endif(DEFINED ENV{COVERAGE})

# Loop through all the tests
foreach(test ${TESTS})

  # Get test information
  get_filename_component(TEST_NAME ${test} NAME)
  get_filename_component(TEST_PATH ${test} PATH)

  # Check for running standard tests (no valgrind, no gcov)
  if(NOT ${MEM_CHECK} AND NOT ${COVERAGE})

    # Check serial/parallel
    if (${MPI_ENABLED})

      # Preform a parallel test
      add_test(NAME ${TEST_NAME}
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND ${PYTHON_EXECUTABLE} ${TEST_NAME} --exe $<TARGET_FILE:openmc>
        --mpi_exec $ENV{MPI_DIR}/bin/mpiexec)

    else(${MPI_ENABLED})

      # Perform a serial test
      add_test(NAME ${TEST_NAME}
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND ${PYTHON_EXECUTABLE} ${TEST_NAME} --exe $<TARGET_FILE:openmc>)

    endif(${MPI_ENABLED})

  # Handle special case for valgrind and gcov (run openmc directly, no python)
  else(NOT ${MEM_CHECK} AND NOT ${COVERAGE})

    # If a plot test is encountered, run with "-p"
    if (${test} MATCHES "test_plot")

      # Perform serial valgrind and coverage test with plot flag
      add_test(NAME ${TEST_NAME}
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND $<TARGET_FILE:openmc> -p ${TEST_PATH})

    elseif(${test} MATCHES "test_filter_distribcell")

      # Add each case for distribcell tests
      add_test(NAME ${TEST_NAME}_case-1
        WORKING_DIRECTORY ${TEST_PATH}/case-1
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH}/case-1)
      add_test(NAME ${TEST_NAME}_case-2
        WORKING_DIRECTORY ${TEST_PATH}/case-2
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH}/case-2)
      add_test(NAME ${TEST_NAME}_case-3
        WORKING_DIRECTORY ${TEST_PATH}/case-3
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH}/case-3)
      add_test(NAME ${TEST_NAME}_case-4
        WORKING_DIRECTORY ${TEST_PATH}/case-4
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH}/case-4)

    # If a restart test is encounted, need to run with -r and restart file(s)
    elseif(${test} MATCHES "restart")

      # Handle restart tests separately
      if(${test} MATCHES "test_statepoint_restart")
        set(RESTART_FILE statepoint.07.h5)
      elseif(${test} MATCHES "test_sourcepoint_restart")
        set(RESTART_FILE statepoint.07.h5 source.07.h5)
      elseif(${test} MATCHES "test_particle_restart_eigval")
        set(RESTART_FILE particle_9_555.h5)
      elseif(${test} MATCHES "test_particle_restart_fixed")
        set(RESTART_FILE particle_7_928.h5)
      else(${test} MATCHES "test_statepoint_restart")
        message(FATAL_ERROR "Restart test ${test} not recognized")
      endif(${test} MATCHES "test_statepoint_restart")

      # Perform serial valgrind and coverage test
      add_test(NAME ${TEST_NAME}
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH})

      # Perform serial valgrind and coverage restart test
      add_test(NAME ${TEST_NAME}_restart
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND $<TARGET_FILE:openmc> -r ${RESTART_FILE} ${TEST_PATH})

      # Set test dependency
      set_tests_properties(${TEST_NAME}_restart PROPERTIES DEPENDS ${TEST_NAME})


    # Handle standard tests for valgrind and gcov
    else(${test} MATCHES "test_plot")

      # Perform serial valgrind and coverage test
      add_test(NAME ${TEST_NAME}
        WORKING_DIRECTORY ${TEST_PATH}
        COMMAND $<TARGET_FILE:openmc> ${TEST_PATH})

    endif(${test} MATCHES "test_plot")

  endif(NOT ${MEM_CHECK} AND NOT ${COVERAGE})

endforeach(test)
