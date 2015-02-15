project(fox_config Fortran)

SET(CMAKE_BUILD_TYPE Debug)

SET(ASSOCIATED_SRCS 
  associated_bug.f90
)

ADD_LIBRARY(
  expr_bug
  ${ASSOCIATED_SRCS}
)

# abort: bare
SET(ABORT_BARE_SRC
  abort_bare.f90
)

ADD_EXECUTABLE(
  abort_bare
  ${ABORT_BARE_SRC}
)

# abort: nag
SET(ABORT_NAG_SRC
  abort_nag.f90
)

ADD_EXECUTABLE(
  abort_nag
  ${ABORT_NAG_SRC}
)

# abort: xlf
SET(ABORT_XLF_SRC
  abort_xlf.f90
)

ADD_EXECUTABLE(
  abort_xlf
  ${ABORT_XLF_SRC}
)

