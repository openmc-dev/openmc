## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)

# Generic information about CDASH site
set(CTEST_PROJECT_NAME "OpenMC")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "openmc.mit.edu")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=OpenMC")
set(CTEST_DROP_SITE_CDASH TRUE)

# Set file size larger to see more output
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "20000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "20000")

# User/password to CDASH site
# Please contact Nick Horelik <nhorelik@mit.edu> or
# Bryan Herman <bherman@mit.edu> if you want to push
# test suite information to our CDASH site.
set(CTEST_DROP_SITE_USER "")
set(CTEST_DROP_SITE_PASSWORD "")
