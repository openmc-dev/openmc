
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(openmc_sources CXX)
add_library(source SHARED parameterized_source_sampling.cpp)
find_package(OpenMC REQUIRED HINTS /home/intern/schumann/src/openmc/openmc-ojs/build)
target_link_libraries(source OpenMC::libopenmc)
