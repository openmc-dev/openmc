# Install script for directory: /mnt/c/Users/bamidele/openmc/vendor/xtl

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xtl" TYPE FILE FILES
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xany.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xbasic_fixed_string.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xbase64.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xclosure.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xcompare.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xcomplex.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xcomplex_sequence.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xspan.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xspan_impl.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xdynamic_bitset.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xfunctional.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xhalf_float.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xhalf_float_impl.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xhash.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xhierarchy_generator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xiterator_base.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xjson.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xmasked_value_meta.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xmasked_value.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xmeta_utils.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xmultimethods.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xoptional_meta.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xoptional.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xoptional_sequence.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xplatform.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xproxy_wrapper.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xsequence.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xsystem.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xtl_config.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xtype_traits.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xvariant.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xvariant_impl.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtl/include/xtl/xvisitor.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtl" TYPE FILE FILES
    "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtl/xtlConfig.cmake"
    "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtl/xtlConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake"
         "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtl/CMakeFiles/Export/2fc63ec57839ed115fc15a5438bb5aec/xtlTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtl/xtlTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtl" TYPE FILE FILES "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtl/CMakeFiles/Export/2fc63ec57839ed115fc15a5438bb5aec/xtlTargets.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pkgconfig" TYPE FILE FILES "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtl/xtl.pc")
endif()

