
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(openmc_sources CXX)
add_library(source SHARED source_sampling.cpp)
find_package(OpenMC REQUIRED HINTS /home/cpf/openmc/build)
target_link_libraries(source OpenMC::libopenmc)
