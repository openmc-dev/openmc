#============================================#
# determine if and how FLUSH intrinsic works #
#============================================#

MESSAGE(STATUS "Determining method to call flush intrinsic")
SET(TEST_FLUSH_OK FALSE)

# try standard (bare) flush call
IF(NOT DEFINED FLUSH)
  TRY_COMPILE(TEST_FLUSH_OK ${fox_BINARY_DIR} ${fox_SOURCE_DIR}/cmake/flush_bare.f90
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_FLUSH_OK} MATCHES TRUE)
    SET(FLUSH "bare")
    SET(FLUSH_METHOD "default")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_FLUSH")
  ENDIF(${TEST_FLUSH_OK} MATCHES TRUE)
ENDIF(NOT DEFINED FLUSH)

# try standard flush with -Vaxlib as command line option
IF(NOT DEFINED FLUSH)
  TRY_COMPILE(TEST_FLUSH_OK ${fox_BINARY_DIR} ${fox_SOURCE_DIR}/cmake/flush_bare.f90
    COMPILE_DEFINITIONS " -Vaxlib "
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_FLUSH_OK} MATCHES TRUE)
    SET(FLUSH "INTEL")
    SET(FLUSH_METHOD "with -Vaxlib")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_FLUSH")
    SET(LDFLAGS "${LD_FLAGS} -Vaxlib")
  ENDIF(${TEST_FLUSH_OK} MATCHES TRUE)
ENDIF(NOT DEFINED FLUSH)

# NAG flush call
IF(NOT DEFINED FLUSH)
  TRY_COMPILE(TEST_FLUSH_OK ${fox_BINARY_DIR} ${fox_SOURCE_DIR}/cmake/flush_nag.f90
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_FLUSH_OK} MATCHES TRUE)
    SET(FLUSH "NAG")
    SET(FLUSH_METHOD "with f90_unix_io")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_FLUSH")
  ENDIF(${TEST_FLUSH_OK} MATCHES TRUE)
ENDIF(NOT DEFINED FLUSH)

# XLF flush call
IF(NOT DEFINED FLUSH)
  TRY_COMPILE(TEST_FLUSH_OK ${fox_BINARY_DIR} ${fox_SOURCE_DIR}/cmake/flush_xlf.f90
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_FLUSH_OK} MATCHES TRUE)
    SET(FLUSH "xlf")
    SET(FLUSH_METHOD "with underscore")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_FLUSH")
  ENDIF(${TEST_FLUSH_OK} MATCHES TRUE)
ENDIF(NOT DEFINED FLUSH)

IF(DEFINED FLUSH)
  MESSAGE("   -> flush intrinsic method is ${FLUSH_METHOD}")
ELSE(DEFINED FLUSH)
  MESSAGE(FATAL_ERROR "   -> ERROR : could not determine how to call FLUSH intrinsic")
ENDIF(DEFINED FLUSH)


