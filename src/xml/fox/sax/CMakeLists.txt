set(fox_sax_srcs
  FoX_sax.f90
)

set(fox_sax_srcs_pp
  m_sax_operate.F90
  m_sax_parser.F90
  m_sax_reader.F90
  m_sax_tokenizer.F90
  m_sax_types.F90
  m_sax_xml_source.F90
)
  

add_library(
  fox_sax STATIC
  ${fox_sax_srcs}
  ${fox_sax_srcs_pp}
)

target_link_libraries(
  fox_sax
  fox_common
  fox_utils
  fox_fsys
)

include_directories(
  ${fox_BINARY_DIR}/common
  ${fox_BINARY_DIR}/utils
  ${fox_BINARY_DIR}/fsys
)

