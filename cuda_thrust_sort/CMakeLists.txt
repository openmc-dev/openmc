cmake_minimum_required(VERSION 3.8.2)

project(openmc_thrust_sort
  DESCRIPTION "Sort EventQueueItem arrays in OpenMC on NVIDIA devices with Thrust"
    LANGUAGES CXX CUDA)

include(FindCUDA) 

list(APPEND openmc_thrust_sort_SOURCES
  src/openmc_thrust_sort.cu)

add_library(openmc_thrust_sort SHARED ${openmc_thrust_sort_SOURCES})
  
target_compile_definitions(openmc_thrust_sort PRIVATE CUDA_THRUST_SORT)

target_include_directories(openmc_thrust_sort PRIVATE ${OPENMC_HEADER_DIR})
#target_include_directories(openmc_thrust_sort PRIVATE ${OPENMC_VERSION_HEADER_DIR})
#target_include_directories(openmc_thrust_sort PRIVATE ${PUGI_DIR})

set_target_properties(
    openmc_thrust_sort
    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    CUDA_ARCHITECTURES "70;80"
)

set_target_properties(openmc_thrust_sort PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
