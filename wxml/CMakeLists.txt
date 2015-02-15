# m_wxml_overloads.F90: m_wxml_overloads.m4
#         m4 -I ../m4 $< > $@
IF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_wxml_overloads.F90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_wxml_overloads.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_wxml_overloads.F90
    DEPENDS m_wxml_overloads.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
ENDIF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")

set(fox_wxml_srcs
  FoX_wxml.f90
)

set(fox_wxml_srcs_pp
  m_wxml_core.F90
  m_wxml_escape.F90
  m_wxml_overloads.F90
)
 
add_library(
  fox_wxml STATIC
  ${fox_wxml_srcs}
  ${fox_wxml_srcs_pp}
)

target_link_libraries(
  fox_wxml
  fox_common
  fox_fsys
)

include_directories(
  ${fox_BINARY_DIR}/fsys
  ${fox_BINARY_DIR}/common
  ${fox_BINARY_DIR}/utils
)

