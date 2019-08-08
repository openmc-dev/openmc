# Create a library called "pyne_source_sampling".
# The extension is already found. Any number of sources could be listed here.

# Find MOAB
find_package(MOAB REQUIRED)

add_library(pyne_source_sampling STATIC
	include/pyne_source_sampling/source_sampling.cpp
	include/pyne_source_sampling/measure.cpp)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to pyne_source_sampling
target_include_directories (pyne_source_sampling PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (pyne_source_sampling PUBLIC ${MOAB_INCLUDE_DIRS})

# link MOAB
target_link_libraries(pyne_source_sampling PUBLIC ${MOAB_LIBRARIES})
