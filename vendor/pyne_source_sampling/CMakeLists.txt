# Create a library called "pyne_source_sampling".

# Find DAGMC/MOAB
find_package(DAGMC REQUIRED)

add_library(pyne_source_sampling STATIC
	include/pyne_source_sampling/source_sampling.cpp
	include/pyne_source_sampling/measure.cpp)

# Make sure the compiler can find include files for our pyne_source_sampling library
# when other libraries or executables link to pyne_source_sampling
target_include_directories (pyne_source_sampling PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (pyne_source_sampling PUBLIC ${DAGMC_INCLUDE_DIRS})

# link MOAB. MOAB is in DAGMC_LIBRARIES
target_link_libraries(pyne_source_sampling PUBLIC ${DAGMC_LIBRARIES})
