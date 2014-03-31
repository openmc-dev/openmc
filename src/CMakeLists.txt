cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(openmc Fortran)

# Setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

#===============================================================================
# Command line options
#===============================================================================

option(openmp   "Enable shared-memory parallelism with OpenMP"   OFF)
option(profile  "Compile with profiling flags"                   OFF)
option(petsc    "Enable PETSC for use in CMFD acceleration"      OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)
option(verbose  "Create verbose Makefiles"                       OFF)

if (verbose)
  set(CMAKE_VERBOSE_MAKEFILE on)
endif()

#===============================================================================
# MPI for distributed-memory parallelism / HDF5 for binary output
#===============================================================================

if($ENV{FC} MATCHES "mpi.*")
  message("-- Detected MPI wrapper: $ENV{FC}")
  add_definitions(-DMPI)
elseif($ENV{FC} MATCHES "h5fc$")
  message("-- Detected HDF5 wrapper: $ENV{FC}")
  add_definitions(-DHDF5)
elseif($ENV{FC} MATCHES "h5pfc$")
  message("-- Detected parallel HDF5 wrapper: $ENV{FC}")
  add_definitions(-DMPI -DHDF5)
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
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
  set(f90flags  "-WF,-DNO_F2008 -O2")
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
# PETSc for CMFD functionality
#===============================================================================

if(petsc)
  find_package(PETSc REQUIRED HINTS $ENV{PETSC_DIR}/conf)
  find_library(libpetsc petsc $ENV{PETSC_DIR}/lib)

  # If libfblas wasn't found, search the PETSc lib directory
  if(PETSC_FBLAS_LIB STREQUAL "PETSC_FBLAS_LIB-NOTFOUND")
    find_library(PETSC_FBLAS_LIB fblas $ENV{PETSC_DIR}/lib)
    list(REMOVE_ITEM PETSC_PACKAGE_LIBS PETSC_FBLAS_LIB-NOTFOUND)
    list(INSERT PETSC_PACKAGE_LIBS 0 ${PETSC_FBLAS_LIB})
  endif()

  # If libflapack wasn't found, search the PETSc lib directory
  if(PETSC_FLAPACK_LIB STREQUAL "PETSC_FLAPACK_LIB-NOTFOUND")
    find_library(PETSC_FLAPACK_LIB flapack $ENV{PETSC_DIR}/lib)
    list(REMOVE_ITEM PETSC_PACKAGE_LIBS PETSC_FLAPACK_LIB-NOTFOUND)
    list(INSERT PETSC_PACKAGE_LIBS 0 ${PETSC_FLAPACK_LIB})
  endif()

  message("-- Using PETSC: ${libpetsc}")
  add_definitions(-DPETSC)
  include_directories($ENV{PETSC_DIR}/include)
  set(libraries "${libpetsc};${PETSC_PACKAGE_LIBS};${libraries}")
endif()

#===============================================================================
# git SHA1 hash
#===============================================================================

execute_process(COMMAND git rev-parse -q HEAD
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

message("-- Initializing/Updating FoX XML submodule...")
execute_process(COMMAND git submodule init
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
execute_process(COMMAND git submodule update
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
file(GLOB_RECURSE source_fox xml/*.F90 xml/*.f90)
add_library(fox STATIC ${source_fox})

#===============================================================================
# Build OpenMC executable
#===============================================================================

set(program "openmc")
file(GLOB source *.F90)
add_executable(${program} ${source})
target_link_libraries(${program} ${libraries} fox)
set_target_properties(${program} PROPERTIES
  COMPILE_FLAGS "${f90flags}"
  LINK_FLAGS "${ldflags}")

#===============================================================================
# Install executable, scripts, manpage, license
#===============================================================================

install(TARGETS ${program} RUNTIME DESTINATION bin)
install(PROGRAMS utils/statepoint_cmp.py
  DESTINATION bin
  RENAME statepoint_cmp)
install(PROGRAMS utils/statepoint_histogram.py
  DESTINATION bin
  RENAME statepoint_histogram)
install(PROGRAMS utils/statepoint_meshplot.py
  DESTINATION bin
  RENAME statepoint_meshplot)
install(FILES ../man/man1/openmc.1 DESTINATION share/man/man1)
install(FILES ../LICENSE DESTINATION "share/doc/${program}/copyright")

find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
  install(CODE "execute_process(
                  COMMAND ${PYTHON_EXECUTABLE} setup.py install --user
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/utils)")
endif()
