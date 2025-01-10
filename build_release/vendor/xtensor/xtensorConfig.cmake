############################################################################
# Copyright (c) Johan Mabille, Sylvain Corlay and Wolf Vollprecht          #
# Copyright (c) QuantStack
#                                                                          #
# Distributed under the terms of the BSD 3-Clause License.                 #
#                                                                          #
# The full license is in the file LICENSE, distributed with this software. #
############################################################################

# xtensor cmake module
# This module sets the following variables in your project::
#
#   xtensor_FOUND - true if xtensor found on the system
#   xtensor_INCLUDE_DIRS - the directory containing xtensor headers
#   xtensor_LIBRARY - empty


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was xtensorConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)
find_dependency(xtl 0.7.5)

if(NOT TARGET xtensor)
    include("${CMAKE_CURRENT_LIST_DIR}/xtensorTargets.cmake")
    get_target_property(xtensor_INCLUDE_DIRS xtensor INTERFACE_INCLUDE_DIRECTORIES)
endif()

if(XTENSOR_USE_XSIMD)
    find_dependency(xsimd )
    target_link_libraries(xtensor INTERFACE xsimd)
    target_compile_definitions(xtensor INTERFACE XTENSOR_USE_XSIMD)
endif()

if(XTENSOR_USE_TBB)
    find_dependency(TBB)
    target_link_libraries(xtensor INTERFACE TBB::tbb)
    target_compile_definitions(xtensor INTERFACE XTENSOR_USE_TBB)
endif()

if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER_EQUAL 3.11)
    if(NOT TARGET xtensor::optimize)
        add_library(xtensor::optimize INTERFACE IMPORTED)
        # Microsoft compiler
        if(CMAKE_${COMPILER_LANGUAGE}_COMPILER_IS_MSVC)
            target_compile_options(xtensor::optimize INTERFACE /EHsc /MP /bigobj)
        # gcc, clang, ...
        else()
            include(CheckCXXCompilerFlag)
            CHECK_CXX_COMPILER_FLAG(-march=native arch_native_supported)
            if(arch_native_supported)
              target_compile_options(xtensor::optimize INTERFACE -march=native)
          endif()
        endif()
    endif()

    if(NOT TARGET xtensor::use_xsimd)
        find_package(xsimd QUIET)
        if (xsimd_FOUND)
            add_library(xtensor::use_xsimd INTERFACE IMPORTED)
            target_link_libraries(xtensor::use_xsimd INTERFACE xsimd)
            target_compile_definitions(xtensor::use_xsimd INTERFACE XTENSOR_USE_XSIMD)
        endif()
    endif()

    if(NOT TARGET xtensor::use_TBB)
        find_package(TBB QUIET)
        if (TBB_FOUND)
            add_library(xtensor::use_TBB INTERFACE IMPORTED)
            target_compile_definitions(xtensor::use_TBB INTERFACE XTENSOR_USE_TBB)
        endif()
    endif()
endif()
