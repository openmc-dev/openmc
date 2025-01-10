############################################################################
# Copyright (c) 2017, Sylvain Corlay and Johan Mabille                     #
#                                                                          #
# Distributed under the terms of the BSD 3-Clause License.                 #
#                                                                          #
# The full license is in the file LICENSE, distributed with this software. #
############################################################################

# xtl cmake module
# This module sets the following variables in your project::
#
#   xtl_FOUND - true if xtl found on the system
#   xtl_INCLUDE_DIRS - the directory containing xtl headers
#   xtl_LIBRARY - empty


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was xtlConfig.cmake.in                            ########

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

if(NOT TARGET xtl)
    include("${CMAKE_CURRENT_LIST_DIR}/xtlTargets.cmake")
    get_target_property(xtl_INCLUDE_DIRS xtl INTERFACE_INCLUDE_DIRECTORIES)
endif()
