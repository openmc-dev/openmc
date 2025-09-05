
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(openmc_sources CXX)
add_library(source SHARED parameterized_source_sampling.cpp)
find_package(OpenMC REQUIRED HINTS /Users/grego01/OpenMC/openmc/build)
target_link_libraries(source OpenMC::libopenmc)
