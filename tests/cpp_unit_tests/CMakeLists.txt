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
  test_tensor
  test_geometry
  # Add additional unit test files here
)

foreach(test ${TEST_NAMES})
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} Catch2::Catch2WithMain libopenmc_test)
  add_test(NAME ${test} COMMAND ${test} WORKING_DIRECTORY ${UNIT_TEST_BIN_OUTPUT_DIR})

  # Catch2 uses exit code 4 when all test cases in a binary are skipped. Treat
  # that result as a skipped CTest test rather than a failure when an optional
  # dependency is unavailable.
  set_tests_properties(${test} PROPERTIES SKIP_RETURN_CODE 4)
endforeach()
