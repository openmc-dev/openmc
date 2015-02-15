

set(m_dom_dom_srcs
  m_dom_types.m4 
  m_dom_node.m4 
  m_dom_nodelist.m4
  m_dom_namednodemap.m4 
  m_dom_implementation.m4
  m_dom_document.m4 
  m_dom_document_type.m4
  m_dom_element.m4 
  m_dom_attribute.m4
  m_dom_character_data.m4 
  m_dom_entity.m4
  m_dom_processing_instruction.m4 
  m_dom_text.m4
  m_dom_common.m4 
  m_dom_exception.m4 
  m_dom_dom.m4
  m_dom_namespaces.m4 
  m_dom_configuration.m4
  m_dom_treewalk.m4 
  m_dom_exception.m4 
  m_dom_object.m4
)

# m_dom_extras.F90: m_dom_extras.m4
#         m4 -I ../m4 $< > $@
IF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_extras.F90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_extras.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_extras.F90
    DEPENDS m_dom_extras.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
  
  # m_dom_parse.f90: m_dom_parse.m4
  #         m4 -I ../m4 $< > $@
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_parse.f90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_parse.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_parse.f90
    DEPENDS m_dom_parse.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
  
  # m_dom_utils.f90: m_dom_utils.m4
  #         m4 -I ../m4 $< > $@
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_utils.f90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_utils.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_utils.f90
    DEPENDS m_dom_utils.m4
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
  
  # m_dom_dom.F90: $(M4FILES)
  #         m4 -I ../m4 m_dom_dom.m4 > $@ 
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_dom.F90
    COMMAND ${M4_PROGRAM} -I ../m4 ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_dom.m4 > ${CMAKE_CURRENT_SOURCE_DIR}/m_dom_dom.F90
    DEPENDS ${m_dom_dom_srcs}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )
ENDIF(NOT ${M4_PROGRAM} MATCHES "NOTFOUND")


set(fox_dom_srcs_pp
  m_dom_dom.F90
  m_dom_extras.F90
)
  
set(fox_dom_srcs
  FoX_dom.f90
  m_dom_error.f90
  m_dom_parse.f90
  m_dom_utils.f90
)

add_library(
  fox_dom STATIC
  ${fox_dom_srcs}
  ${fox_dom_srcs_pp}
)

## add_dependencies(
##   fox_dom
##   dom_generate
## )

target_link_libraries(
  fox_dom
  fox_wxml
  fox_sax
)

include_directories(
  ${fox_BINARY_DIR}/dom
  ${fox_BINARY_DIR}/wxml
  ${fox_BINARY_DIR}/common
  ${fox_BINARY_DIR}/utils
  ${fox_BINARY_DIR}/fsys
  ${fox_BINARY_DIR}/sax
)



