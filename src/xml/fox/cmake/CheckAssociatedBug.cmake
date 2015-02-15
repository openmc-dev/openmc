#=====================================================#
# check for 'associated in restricted expression' bug #
#=====================================================#

MESSAGE(STATUS "Checking for 'associated in restricted expression' bug")
TRY_COMPILE(TEST_ASSOCIATED ${fox_BINARY_DIR}/associated ${fox_SOURCE_DIR}/cmake
  fox_config expr_bug
  OUTPUT_VARIABLE BUILD_OUTPUT 
)

IF(${TEST_ASSOCIATED} MATCHES FALSE)
  MESSAGE("   -> yes")
  SET(FPPFLAGS "${FPPFLAGS} -DRESTRICTED_ASSOCIATED_BUG")
ENDIF(${TEST_ASSOCIATED} MATCHES FALSE)


