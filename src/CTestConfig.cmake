## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "OpenMC")
set(CTEST_NIGHTLY_START_TIME "05:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "neutronbalance.mit.edu")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=OpenMC")
set(CTEST_DROP_SITE_CDASH TRUE)
