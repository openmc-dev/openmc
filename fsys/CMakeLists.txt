# fox_m_fsys_parse_input.F90: fox_m_fsys_parse_input.m4
#         m4 -I ../m4 $< > $@
IF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/fox_m_fsys_parse_input.F90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/fox_m_fsys_parse_input.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/fox_m_fsys_parse_input.F90
    DEPENDS fox_m_fsys_parse_input.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
ENDIF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")

set(fox_fsys_srcs
  fox_m_fsys_realtypes.f90
#  m_ieee.f90
)

set(fox_fsys_srcs_pp
  fox_m_fsys_abort_flush.F90
  fox_m_fsys_array_str.F90
  fox_m_fsys_count_parse_input.F90
  fox_m_fsys_format.F90
  fox_m_fsys_parse_input.F90
  fox_m_fsys_string.F90
  fox_m_fsys_string_list.F90
  fox_m_fsys_varstr.F90
)

add_library(
  fox_fsys STATIC
  ${fox_fsys_srcs}
  ${fox_fsys_srcs_pp}
)

