cmake_minimum_required(VERSION 3.8.2)
# Search for rocm in common locations
if(NOT DEFINED HIP_PATH)
    if(NOT DEFINED ENV{HIP_PATH})
        set(HIP_PATH "/opt/rocm/hip" CACHE PATH "Path to which HIP has been installed")
    else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
    endif()
endif()
set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
enable_language(HIP)

project(openmc_thrust_sort
  DESCRIPTION "Sort EventQueueItem arrays in OpenMC on AMD devices with Thrust"
  LANGUAGES CXX HIP)

list(APPEND openmc_thrust_sort_SOURCES
  src/openmc_thrust_sort.hip)

add_library(openmc_thrust_sort SHARED ${openmc_thrust_sort_SOURCES})
  
target_compile_definitions(openmc_thrust_sort PRIVATE CUDA_THRUST_SORT)

target_include_directories(openmc_thrust_sort PRIVATE ${OPENMC_HEADER_DIR})

set_target_properties(
    openmc_thrust_sort
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    COMPILE_FLAGS -O3
)
