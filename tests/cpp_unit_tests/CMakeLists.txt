set(TEST_NAMES
  test_atomic_mass
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
  test_particle
  test_secondary_bank
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

# Also run the shared secondary bank test on two MPI ranks, so that banked
# sites are sent between ranks
if(OPENMC_USE_MPI)
  add_test(NAME test_secondary_bank_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_secondary_bank> ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${UNIT_TEST_BIN_OUTPUT_DIR})
endif()

if(OPENMC_USE_XDG)
  target_compile_definitions(test_mesh PRIVATE OPENMC_XDG_ENABLED)
endif()
