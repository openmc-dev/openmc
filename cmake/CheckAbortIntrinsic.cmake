#============================================#
# determine if and how ABORT intrinsic works #
#============================================#

MESSAGE(STATUS "Determining method to call abort intrinsic")
SET(TEST_ABORT_OK FALSE) 

# intel abort
IF(NOT DEFINED ABORT)
  TRY_COMPILE(TEST_ABORT_INTEL ${fox_BINARY_DIR}/abort_intel ${fox_SOURCE_DIR}/cmake
    fox_config abort_intel
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_ABORT_INTEL} MATCHES TRUE)
    MESSAGE(" abort : intel works")
    SET(ABORT "INTEL")
    SET(ABORT_METHOD "with argument")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_ABORT -DFC_ABORT_ARG")
  ENDIF(${TEST_ABORT_INTEL} MATCHES TRUE)
ENDIF(NOT DEFINED ABORT)

# xlf abort
IF(NOT DEFINED ABORT)
  TRY_COMPILE(TEST_ABORT_XLF ${fox_BINARY_DIR}/abort_xlf ${fox_SOURCE_DIR}/cmake
    fox_config abort_xlf
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_ABORT_XLF} MATCHES TRUE)
    MESSAGE(" abort : xlf works")
    SET(ABORT "XLF")
    SET(ABORT_METHOD "with underscore")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_ABORT -DFC_ABORT_UNDERSCORE")
  ENDIF(${TEST_ABORT_XLF} MATCHES TRUE)
ENDIF(NOT DEFINED ABORT)

# bare abort
IF(NOT DEFINED ABORT)
  TRY_COMPILE(TEST_ABORT_BARE ${fox_BINARY_DIR}/abort_bare ${fox_SOURCE_DIR}/cmake
    fox_config abort_bare
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_ABORT_BARE} MATCHES TRUE)
    MESSAGE(" abort : bare works")
    SET(ABORT "BARE")
    SET(ABORT_METHOD "default")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_ABORT")
  ENDIF(${TEST_ABORT_BARE} MATCHES TRUE)
ENDIF(NOT DEFINED ABORT)

# nag abort
IF(NOT DEFINED ABORT)
  TRY_COMPILE(TEST_ABORT_NAG ${fox_BINARY_DIR}/abort_nag ${fox_SOURCE_DIR}/cmake
    fox_config abort_nag
    OUTPUT_VARIABLE BUILD_OUTPUT 
  )
  IF(${TEST_ABORT_NAG} MATCHES TRUE)
    MESSAGE(" abort : nag works")
    SET(ABORT "NAG")
    SET(ABORT_METHOD "with f90_unix_proc")
    SET(FPPFLAGS "${FPPFLAGS} -DFC_HAVE_ABORT")
  ENDIF(${TEST_ABORT_NAG} MATCHES TRUE)
ENDIF(NOT DEFINED ABORT)


