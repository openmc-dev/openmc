
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(openmc_cpp_driver CXX)
add_executable(cpp_driver driver.cpp)
find_package(OpenMC REQUIRED HINTS /home/intern/schumann/src/openmc/openmc-ojs/build)
target_link_libraries(cpp_driver OpenMC::libopenmc)
