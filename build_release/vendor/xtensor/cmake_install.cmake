# Install script for directory: /mnt/c/Users/bamidele/openmc/vendor/xtensor

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/xtensor" TYPE FILE FILES
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xaccessible.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xaccumulator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xadapt.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xarray.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xassign.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xaxis_iterator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xaxis_slice_iterator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xblockwise_reducer.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xblockwise_reducer_functors.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xbroadcast.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xbuffer_adaptor.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xbuilder.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xchunked_array.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xchunked_assign.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xchunked_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xcomplex.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xcontainer.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xcsv.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xdynamic_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xeval.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xexception.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xexpression.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xexpression_holder.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xexpression_traits.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xfixed.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xfunction.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xfunctor_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xgenerator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xhistogram.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xindex_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xinfo.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xio.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xiterable.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xiterator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xjson.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xlayout.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xmanipulation.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xmasked_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xmath.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xmime.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xmultiindex_iterator.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xnoalias.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xnorm.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xnpy.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoffset_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoperation.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoptional.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoptional_assembly.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoptional_assembly_base.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xoptional_assembly_storage.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xpad.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xrandom.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xreducer.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xrepeat.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xscalar.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xsemantic.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xset_operation.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xshape.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xslice.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xsort.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xstorage.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xstrided_view.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xstrided_view_base.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xstrides.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xtensor.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xtensor_config.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xtensor_forward.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xtensor_simd.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xutils.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xvectorize.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xview.hpp"
    "/mnt/c/Users/bamidele/openmc/vendor/xtensor/include/xtensor/xview_utils.hpp"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor" TYPE FILE FILES
    "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtensor/xtensorConfig.cmake"
    "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtensor/xtensorConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake"
         "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtensor/CMakeFiles/Export/0f2a327e949144b0b747c96acb1eb12f/xtensorTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor/xtensorTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/xtensor" TYPE FILE FILES "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtensor/CMakeFiles/Export/0f2a327e949144b0b747c96acb1eb12f/xtensorTargets.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/pkgconfig" TYPE FILE FILES "/mnt/c/Users/bamidele/openmc/build_release/vendor/xtensor/xtensor.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/mnt/c/Users/bamidele/openmc/build_release/xtensor.hpp")
endif()

