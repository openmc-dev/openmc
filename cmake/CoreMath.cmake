# Correctly rounded transcendental functions from the CORE-MATH project
# (https://core-math.gitlabpages.inria.fr/, MIT license), used when
# OPENMC_USE_COREMATH is ON. See include/openmc/math.h and the developer's
# guide for why.
#
# CORE-MATH has no releases and no build system of its own; each function is a
# self-contained C file. The files OpenMC needs are fetched from the upstream
# repository at a pinned commit, verified against the SHA-256 hashes recorded
# here, and compiled into a static library. Point OPENMC_COREMATH_SOURCE_DIR at
# a checkout of the repository to build without network access.

set(COREMATH_COMMIT dae5dd8c0d5bd2691b45de28c9edbe12f067c3d3)
set(COREMATH_RAW_URL
  "https://gitlab.inria.fr/core-math/core-math/-/raw/${COREMATH_COMMIT}")

set(OPENMC_COREMATH_SOURCE_DIR "" CACHE PATH
  "Local checkout of the CORE-MATH repository to use instead of downloading")

# Relative path below src/binary64 and SHA-256 of each file at COREMATH_COMMIT
set(COREMATH_FILES
  acos/acos.c     ecd6173837214367f0140c038c615a4f37d27795c0cc74e16bdaa59f7df8cff1
  asin/asin.c     08c4d31ad429e4aac3cc1adf9dbc6e31d59b6ff997626fbfcd3b60ebe2cbffe8
  atan/atan.c     2d62a4d780f233fd0460b8e6076880fcb9d2f55280ead5799180983b6b77c863
  atan2/atan2.c   404b0db90e536f8248101544d0fcad9141fd76b2a6cc3d189dd80dbfc865bd91
  atan2/tint.h    41264d54812a8ce2b3f8ca4c8c033094df7cbf2d214da1c7555783328171b921
  cbrt/cbrt.c     0a86b7a8f433dfa4610023767102ea19eb3642ea6982871c61c956148837dde4
  cos/cos.c       04bffbc8b483c886f78f2c9ab3f2a1ca318dda78ab637a2983a4a9d82374e9d9
  cosh/cosh.c     aa9e80b3fe3fa4aa0b3b2abf11a9fe184f0170c5faa042ad80dcd4698417842b
  erf/erf.c       a984176d6b175699c0da69827776f6dab98fad9cf17fd8c8af53558f2c85b39f
  erfc/erfc.c     53b96001983317d267b3d9867995a74fe21fffe44931d53cc51d0fd688fe9964
  exp/exp.c       55ad39a6610e0f46794be6b70a1ee4d76b10dd4f64d47b61c4e918d711e3041d
  expm1/expm1.c   f33b7d36cdfe319d52d57cb048d8f72ccaf1e9b7038dd74672e0c8ee6a99c033
  hypot/hypot.c   71d47330651a699e1b2791d7654cd4683f55b203bcc425b5b48f2fa2a5b379d0
  lgamma/lgamma.c edda41084ffb82fce0e6be72aac5fe7d977332e739ce92338cec6c9677e1f6de
  log/log.c       922af9ae4b4bd009bf13238ac8a20c1be25d1fee912b0c35dbfd4962c1572497
  log/dint.h      0da80a15d311bb1e550035f98edc2266ef31828cd423aec3a2e4507d2dad5e0e
  log1p/log1p.c   068da268a928d8f18792e6a9ae8a0c54291925788b4e8c3b3c3c33daac220666
  log2/log2.c     2097e97c8955de82c0c5df22997ba13617e148ef691c2a86f01ece4cc2310ace
  pow/pow.c       b14c65a1d1367f49a0d3e00e65370c3e9dd0fd004fe00e30cdad2a9887954562
  pow/pow.h       8e28930b514f1f439011cb201c738c6b8940f321f03f4679362e7f6b7893fd6a
  pow/dint.h      c62a9788d872c314c743a821602d6819e5c13113b70f6d4a5d0076695dc413ef
  pow/qint.h      32de9c10388056e2294339e3302c3695cdf4abc87ca790ad48917d9794bfb6e4
  sin/sin.c       cd772a06d75d1bb8bb4374847065acfab28d5a06c7fb81be3b36127f65d8b1c2
  sinh/sinh.c     bafb6b0e42c69607ee36d00b5d72724c4662881e1152cdc0adc02c4bdbf3edb6
  tan/tan.c       dd662b7af88a86d85e01b88c57f9df4b6019ddfd02120ca381b5ba42283987ae
  tgamma/tgamma.c 61536c08333b468f7702e8b3bf27039e5fa7bf9b6c78c53c8875778018d0f496
)
set(COREMATH_LICENSE_SHA256
  81830875362f4d477637de58eecd6d458db108cc90c40904fca3694807b86de3)

if(MSVC)
  message(FATAL_ERROR "OPENMC_USE_COREMATH requires a compiler with 128-bit "
    "integer support (GCC or Clang); it is not available with MSVC.")
endif()

if(OPENMC_COREMATH_SOURCE_DIR)
  # A local checkout must hold the pinned content: the hashes are checked
  # exactly as they are for downloaded files.
  set(COREMATH_DIR "${OPENMC_COREMATH_SOURCE_DIR}")
  message(STATUS "CORE-MATH: using sources in ${COREMATH_DIR}")
  list(LENGTH COREMATH_FILES _n)
  math(EXPR _last "${_n} - 1")
  foreach(_i RANGE 0 ${_last} 2)
    list(GET COREMATH_FILES ${_i} _path)
    math(EXPR _j "${_i} + 1")
    list(GET COREMATH_FILES ${_j} _sha)
    set(_file "${COREMATH_DIR}/src/binary64/${_path}")
    if(NOT EXISTS "${_file}")
      message(FATAL_ERROR "CORE-MATH: ${_file} not found")
    endif()
    file(SHA256 "${_file}" _actual)
    if(NOT _actual STREQUAL _sha)
      message(FATAL_ERROR "CORE-MATH: ${_file} does not match the pinned "
        "commit ${COREMATH_COMMIT} (SHA-256 ${_actual}, expected ${_sha})")
    endif()
  endforeach()
else()
  set(COREMATH_DIR "${CMAKE_BINARY_DIR}/core-math-${COREMATH_COMMIT}")
  message(STATUS "CORE-MATH: fetching sources at commit ${COREMATH_COMMIT}")
  list(LENGTH COREMATH_FILES _n)
  math(EXPR _last "${_n} - 1")
  foreach(_i RANGE 0 ${_last} 2)
    list(GET COREMATH_FILES ${_i} _path)
    math(EXPR _j "${_i} + 1")
    list(GET COREMATH_FILES ${_j} _sha)
    set(_dest "${COREMATH_DIR}/src/binary64/${_path}")
    file(DOWNLOAD "${COREMATH_RAW_URL}/src/binary64/${_path}" "${_dest}"
      EXPECTED_HASH SHA256=${_sha} STATUS _status TLS_VERIFY ON)
    list(GET _status 0 _code)
    if(NOT _code EQUAL 0)
      list(GET _status 1 _msg)
      message(FATAL_ERROR "CORE-MATH: failed to fetch ${_path}: ${_msg}. "
        "Set OPENMC_COREMATH_SOURCE_DIR to a local checkout to build offline.")
    endif()
  endforeach()
  file(DOWNLOAD "${COREMATH_RAW_URL}/LICENSE" "${COREMATH_DIR}/LICENSE"
    EXPECTED_HASH SHA256=${COREMATH_LICENSE_SHA256} STATUS _status
    TLS_VERIFY ON)
  list(GET _status 0 _code)
  if(NOT _code EQUAL 0)
    message(FATAL_ERROR "CORE-MATH: failed to fetch LICENSE")
  endif()
endif()

set(COREMATH_SOURCES)
list(LENGTH COREMATH_FILES _n)
math(EXPR _last "${_n} - 1")
foreach(_i RANGE 0 ${_last} 2)
  list(GET COREMATH_FILES ${_i} _path)
  if(_path MATCHES "\\.c$")
    list(APPEND COREMATH_SOURCES "${COREMATH_DIR}/src/binary64/${_path}")
  endif()
endforeach()

add_library(coremath STATIC ${COREMATH_SOURCES})
set_target_properties(coremath PROPERTIES
  C_STANDARD 11
  C_EXTENSIONS ON
  POSITION_INDEPENDENT_CODE ON)
# Each function includes headers from its own directory
foreach(_src ${COREMATH_SOURCES})
  get_filename_component(_dir "${_src}" DIRECTORY)
  set_property(SOURCE "${_src}" APPEND PROPERTY INCLUDE_DIRECTORIES "${_dir}")
endforeach()

# The error-free transformations in CORE-MATH assume no FMA contraction; the
# explicit fused multiply-adds it uses are unaffected by this flag.
include(CheckCCompilerFlag)
check_c_compiler_flag(-ffp-contract=off COREMATH_SUPPORTS_FP_CONTRACT_OFF)
if(COREMATH_SUPPORTS_FP_CONTRACT_OFF)
  target_compile_options(coremath PRIVATE -ffp-contract=off)
endif()

# CORE-MATH calls fma() in its inner loops. On x86-64 the baseline ISA has no
# FMA instruction, so unless the compiler is allowed to emit it every call goes
# to a software routine and the functions are several times slower. The results
# are identical either way, since fma is correctly rounded in both cases.
#
# The instruction set of the rest of the build is never raised here: AUTO adds
# -mfma only when the compiler flags in use already target FMA (for example
# -march=native or -march=x86-64-v3), so a binary built for the baseline ISA
# stays runnable on every x86-64 CPU. ON forces the flag for the CORE-MATH
# sources; the resulting binary then requires an FMA-capable CPU (2013 or
# later).
set(OPENMC_COREMATH_FMA "AUTO" CACHE STRING
  "Compile CORE-MATH with the FMA instruction on x86-64 (AUTO, ON or OFF)")
set_property(CACHE OPENMC_COREMATH_FMA PROPERTY STRINGS AUTO ON OFF)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
  check_c_compiler_flag(-mfma COREMATH_SUPPORTS_MFMA)
  if(COREMATH_SUPPORTS_MFMA)
    if(OPENMC_COREMATH_FMA STREQUAL "AUTO")
      include(CheckCSourceCompiles)
      set(CMAKE_REQUIRED_FLAGS "${CMAKE_C_FLAGS}")
      check_c_source_compiles("
        #ifndef __FMA__
        #error FMA is not part of the target instruction set
        #endif
        int main(void) { return 0; }" COREMATH_TARGET_HAS_FMA)
      unset(CMAKE_REQUIRED_FLAGS)
      set(_use_fma ${COREMATH_TARGET_HAS_FMA})
    else()
      set(_use_fma ${OPENMC_COREMATH_FMA})
    endif()
    if(_use_fma)
      target_compile_options(coremath PRIVATE -mfma)
      message(STATUS "CORE-MATH: compiling with -mfma")
    else()
      message(STATUS "CORE-MATH: compiling without -mfma (software fma); set "
        "OPENMC_COREMATH_FMA=ON or build with -march=native for FMA hardware")
    endif()
  endif()
endif()

target_link_libraries(coremath PUBLIC m)
