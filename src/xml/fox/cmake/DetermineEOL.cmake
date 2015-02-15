#===========================================================#
# determine end-of-line character based on operating system #
#===========================================================#

MESSAGE(STATUS "Determining end-of-line character by the host name")
IF(${CMAKE_HOST_APPLE})
  SET(FPPFLAGS ${FPPFLAGS} "-DFC_EOR_CR")
  SET(EOLCHAR "CR")
ELSEIF(${CMAKE_HOST_WIN32})
  SET(FPPFLAGS ${FPPFLAGS} "-DFC_EOR_CRLF")
  SET(EOLCHAR "CRLF")
ELSEIF(${CMAKE_HOST_UNIX})
  SET(FPPFLAGS ${FPPFLAGS} "-DFC_EOR_LF")
  SET(EOLCHAR "LF")
ELSE(${CMAKE_HOST_APPLE})
  SET(FPPFLAGS ${FPPFLAGS})
  MESSAGE(STATUS "warning: could not determine host system: not unix, not win32 and not mac")
  SET(EOLCHAR "")
ENDIF(${CMAKE_HOST_APPLE})
MESSAGE("   -> end-of-line character is ${EOLCHAR}")


