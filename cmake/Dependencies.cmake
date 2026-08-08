include(FetchContent)

# Use the download time for extracted files when supported. This avoids stale
# dependency build products if a URL changes while retaining CMake 3.22 support.
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Declare all dependencies before making any of them available. Parent projects
# can then override these declarations according to FetchContent's
# first-to-declare behavior.
FetchContent_Declare(
  pugixml
  URL https://github.com/zeux/pugixml/archive/ee86beb30e4973f5feffe3ce63bfa4fbadf72f38.tar.gz
  URL_HASH SHA256=51c102d4187fac99daa38af281b0772c5e6c586f65004cdc63f8f2e011a21492
)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/0c9fce2ffefecfdce794e1859584e25877b7b592.tar.gz
  URL_HASH SHA256=f94052c10b611fd374194ca6e0dc4d159459c0b370abfe9002c13058863b7039
)
FetchContent_Declare(
  Catch2
  URL https://github.com/catchorg/Catch2/archive/5a40b2275caa05cf809bf04df848764a9d7df2e2.tar.gz
  URL_HASH SHA256=be038aac877893ea0fa02cdb5f24a46db03085b7f053c8b78ef7bd437c8c6c22
)

function(openmc_find_or_fetch name target)
  set(one_value_args LEGACY_TARGET VERSION)
  cmake_parse_arguments(DEPENDENCY "" "${one_value_args}" "" ${ARGN})

  if(TARGET "${target}")
    message(STATUS "Using existing ${name} target ${target}")
    return()
  endif()

  if(NOT OPENMC_FORCE_FETCHCONTENT)
    if(DEPENDENCY_VERSION)
      find_package(${name} ${DEPENDENCY_VERSION} CONFIG QUIET
        NO_SYSTEM_ENVIRONMENT_PATH)
    else()
      find_package(${name} CONFIG QUIET NO_SYSTEM_ENVIRONMENT_PATH)
    endif()
  endif()

  if(DEPENDENCY_LEGACY_TARGET
     AND TARGET "${DEPENDENCY_LEGACY_TARGET}"
     AND NOT TARGET "${target}")
    add_library("${target}" ALIAS "${DEPENDENCY_LEGACY_TARGET}")
  endif()

  if(TARGET "${target}")
    set(version_variable "${name}_VERSION")
    if(DEFINED ${version_variable})
      message(STATUS "Found ${name} ${${version_variable}}")
    else()
      message(STATUS "Found ${name}")
    endif()
    return()
  endif()

  message(STATUS "Fetching ${name}")
  FetchContent_MakeAvailable(${name})

  if(DEPENDENCY_LEGACY_TARGET
     AND TARGET "${DEPENDENCY_LEGACY_TARGET}"
     AND NOT TARGET "${target}")
    add_library("${target}" ALIAS "${DEPENDENCY_LEGACY_TARGET}")
  endif()

  if(NOT TARGET "${target}")
    message(FATAL_ERROR "${name} did not provide expected target ${target}")
  endif()
endfunction()

openmc_find_or_fetch(pugixml pugixml::pugixml LEGACY_TARGET pugixml)
openmc_find_or_fetch(fmt fmt::fmt)

if(OPENMC_BUILD_TESTS)
  openmc_find_or_fetch(Catch2 Catch2::Catch2WithMain VERSION 3)
endif()
