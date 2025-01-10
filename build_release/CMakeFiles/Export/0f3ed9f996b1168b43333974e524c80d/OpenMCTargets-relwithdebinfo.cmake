#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "OpenMC::openmc" for configuration "RelWithDebInfo"
set_property(TARGET OpenMC::openmc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(OpenMC::openmc PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/bin/openmc"
  )

list(APPEND _cmake_import_check_targets OpenMC::openmc )
list(APPEND _cmake_import_check_files_for_OpenMC::openmc "${_IMPORT_PREFIX}/bin/openmc" )

# Import target "OpenMC::libopenmc" for configuration "RelWithDebInfo"
set_property(TARGET OpenMC::libopenmc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(OpenMC::libopenmc PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libopenmc.so"
  IMPORTED_SONAME_RELWITHDEBINFO "libopenmc.so"
  )

list(APPEND _cmake_import_check_targets OpenMC::libopenmc )
list(APPEND _cmake_import_check_files_for_OpenMC::libopenmc "${_IMPORT_PREFIX}/lib/libopenmc.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
