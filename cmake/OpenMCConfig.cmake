get_filename_component(OpenMC_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)

find_package(xtl REQUIRED HINTS ${OpenMC_CMAKE_DIR}/../xtl)
find_package(xtensor REQUIRED HINTS ${OpenMC_CMAKE_DIR}/../xtensor)

if(NOT TARGET OpenMC::libopenmc)
  include("${OpenMC_CMAKE_DIR}/OpenMCTargets.cmake")
endif()
