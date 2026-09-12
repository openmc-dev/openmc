include(FetchContent)

# Use the download time for extracted files when supported. This avoids stale
# dependency build products if a URL changes while retaining CMake 3.22 support.
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Declare every dependency before making any of them available. FetchContent
# honors the first declaration for a given name, so a project that pulls OpenMC
# in as a subproject can override any of these by declaring them beforehand.
FetchContent_Declare(
  pugixml
  URL https://github.com/zeux/pugixml/archive/refs/tags/v1.15.tar.gz
  URL_HASH SHA256=b39647064d9e28297a34278bfb897092bf33b7c487906ddfc094c9e8868bddcb
)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/11.0.2.tar.gz
  URL_HASH SHA256=6cb1e6d37bdcb756dbbe59be438790db409cdb4868c66e888d5df9f13f7c027f
)
FetchContent_Declare(
  Catch2
  URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.16.0.tar.gz
  URL_HASH SHA256=0957cae5821b17ce07f0833aaa52b5137643a8382203221f363a8303c109af34
)

# Make a dependency available as ${target}, preferring an installed package and
# falling back to the pinned sources declared above. LEGACY_TARGET names an
# unnamespaced target exported by older releases of a dependency, which is
# aliased to ${target} so that callers only ever refer to the namespaced name.
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
  endif()

  message(STATUS "Fetching ${name}")
  FetchContent_MakeAvailable(${name})

  if(NOT TARGET "${target}")
    message(FATAL_ERROR "${name} did not provide expected target ${target}")
  endif()
endfunction()

openmc_find_or_fetch(pugixml pugixml::pugixml LEGACY_TARGET pugixml)
openmc_find_or_fetch(fmt fmt::fmt)

if(OPENMC_BUILD_TESTS)
  openmc_find_or_fetch(Catch2 Catch2::Catch2WithMain VERSION 3)
endif()
