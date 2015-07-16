cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(openmc Fortran)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

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
option(verbose  "Create verbose Makefiles"                       OFF)
option(coverage "Compile with coverage analysis flags"           OFF)
option(mpif08   "Use Fortran 2008 MPI interface"                 OFF)

if (verbose)
  set(CMAKE_VERBOSE_MAKEFILE on)
endif()

# Maximum number of nested coordinates levels
set(maxcoord 10 CACHE STRING "Maximum number of nested coordinate levels")
add_definitions(-DMAX_COORD=${maxcoord})

#===============================================================================
# MPI for distributed-memory parallelism / HDF5 for binary output
#===============================================================================

set(MPI_ENABLED FALSE)
set(HDF5_ENABLED FALSE)
if($ENV{FC} MATCHES "mpi[^/]*$")
  message("-- Detected MPI wrapper: $ENV{FC}")
  add_definitions(-DMPI)
  set(MPI_ENABLED TRUE)
elseif($ENV{FC} MATCHES "h5fc$")
  message("-- Detected HDF5 wrapper: $ENV{FC}")
  add_definitions(-DHDF5)
  set(HDF5_ENABLED TRUE)
elseif($ENV{FC} MATCHES "h5pfc$")
  message("-- Detected parallel HDF5 wrapper: $ENV{FC}")
  add_definitions(-DMPI -DHDF5)
  set(MPI_ENABLED TRUE)
  set(HDF5_ENABLED TRUE)
endif()

# Check for Fortran 2008 MPI interface
if(MPI_ENABLED AND mpif08)
  message("-- Using Fortran 2008 MPI bindings")
  add_definitions(-DMPIF08)
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Make sure version is sufficient
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)
  if(GCC_VERSION VERSION_LESS 4.6)
    message(FATAL_ERROR "gfortran version must be 4.6 or higher")
  endif()

  # GNU Fortran compiler options
  set(f90flags  "-cpp -std=f2008 -fbacktrace")
  if(debug)
    set(f90flags  "-g -Wall -pedantic -fbounds-check -ffpe-trap=invalid,overflow,underflow ${f90flags}")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-fopenmp ${f90flags}")
    set(ldflags   "-fopenmp ${ldflags}")
  endif()
  if(coverage)
    set(f90flags "-coverage ${f90flags}")
    set(ldflags  "-coverage ${ldflags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Intel Fortran compiler options
  set(f90flags  "-fpp -std08 -assume byterecl -traceback")
  if(debug)
    set(f90flags "-g -warn -ftrapuv -fp-stack-check -check all -fpe0 ${f90flags}")
    set(ldflags  "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-openmp ${f90flags}")
    set(ldflags   "-openmp ${ldflags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  # PGI Fortran compiler options
  set(f90flags  "-Mpreprocess -Minform=inform -traceback")
  add_definitions(-DNO_F2008)
  if(debug)
    set(f90flags  "-g -Mbounds -Mchkptr -Mchkstk ${f90flags}")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-fast -Mipa ${f90flags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
  # IBM XL compiler options
  set(f90flags  "-O2")
  add_definitions(-DNO_F2008)
  if(debug)
    set(f90flags  "-g -C -qflag=i:i -u")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-p ${f90flags}")
    set(ldflags   "-p ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-qsmp=omp ${f90flags}")
    set(ldflags   "-qsmp=omp ${ldflags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  # Cray Fortran compiler options
  set(f90flags  "-e Z -m 0")
  if(debug)
    set(f90flags  "-g -R abcnsp -O0 ${f90flags}")
    set(ldflags   "-g")
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
target_link_libraries(${program} ${libraries} fox_dom)
set_target_properties(${program} PROPERTIES
  COMPILE_FLAGS "${f90flags}"
  LINK_FLAGS "${ldflags}")

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

    # If a restart test is encounted, need to run with -r and restart file(s)
    elseif(${test} MATCHES "restart")

      # Set restart file names
      if (${HDF5_ENABLED})

        # Handle restart tests separately
        if(${test} MATCHES "test_statepoint_restart")
          set(RESTART_FILE statepoint.07.h5)
        elseif(${test} MATCHES "test_sourcepoint_restart")
          set(RESTART_FILE statepoint.07.h5 source.07.h5)
        elseif(${test} MATCHES "test_particle_restart_eigval")
          set(RESTART_FILE particle_12_616.h5)
        elseif(${test} MATCHES "test_particle_restart_fixed")
          set(RESTART_FILE particle_7_6144.h5)
        else(${test} MATCHES "test_statepoint_restart")
          message(FATAL_ERROR "Restart test ${test} not recognized")
        endif(${test} MATCHES "test_statepoint_restart")

      else(${HDF5_ENABLED})

        # Handle restart tests separately
        if(${test} MATCHES "test_statepoint_restart")
          set(RESTART_FILE statepoint.07.binary)
        elseif(${test} MATCHES "test_sourcepoint_restart")
          set(RESTART_FILE statepoint.07.binary source.07.binary)
        elseif(${test} MATCHES "test_particle_restart_eigval")
          set(RESTART_FILE particle_12_616.binary)
        elseif(${test} MATCHES "test_particle_restart_fixed")
          set(RESTART_FILE particle_7_6144.binary)
        else(${test} MATCHES "test_statepoint_restart")
          message(FATAL_ERROR "Restart test ${test} not recognized")
        endif(${test} MATCHES "test_statepoint_restart")

      endif(${HDF5_ENABLED})

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
