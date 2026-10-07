set(TEST_NAMES
  test_distribution
  test_file_utils
  test_tally
  test_interpolate
  test_math
  test_mcpl_stat_sum
  test_mesh
  test_photon
  test_ray
  test_region
  test_surface
  test_tensor
  test_geometry
  # Add additional unit test files here
)

foreach(test ${TEST_NAMES})
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} Catch2::Catch2WithMain libopenmc_test)
  add_test(NAME ${test} COMMAND ${test} WORKING_DIRECTORY ${UNIT_TEST_BIN_OUTPUT_DIR})
endforeach()

if(OPENMC_USE_MPI)
  add_executable(test_message_passing test_message_passing.cpp)
  target_link_libraries(test_message_passing PRIVATE Catch2::Catch2 libopenmc_test)
  add_test(NAME test_message_passing
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_message_passing>
            ${MPIEXEC_POSTFLAGS})
  set_tests_properties(test_message_passing PROPERTIES TIMEOUT 60 PROCESSORS 2)
endif()
