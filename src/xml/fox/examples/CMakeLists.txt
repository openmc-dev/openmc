add_executable(
  wxml_example
  wxml_example.f90)
target_link_libraries(
  wxml_example
  fox_wxml)

add_executable(
  wcml_example
  wcml_example.f90)
target_link_libraries(
  wcml_example
  fox_wcml)

add_executable(
  wkml_example
  wkml_example.f90)
target_link_libraries(
  wkml_example
  fox_wkml)

add_executable(
  sax_example
  sax_example.f90)
target_link_libraries(
  sax_example
  fox_sax)

add_executable(
  sax_example_2
  sax_example_2.f90)
target_link_libraries(
  sax_example_2
  fox_sax)

add_executable(
  dom_example_2
  dom_example_2.f90)
target_link_libraries(
  dom_example_2
  fox_dom)

add_executable(
  dom_example_3
  dom_example_3.f90)
target_link_libraries(
  dom_example_3
  fox_dom)


include_directories(
  ${fox_BINARY_DIR}/dom
  ${fox_BINARY_DIR}/wxml
  ${fox_BINARY_DIR}/wcml
  ${fox_BINARY_DIR}/wkml
  ${fox_BINARY_DIR}/common
  ${fox_BINARY_DIR}/utils
  ${fox_BINARY_DIR}/fsys
  ${fox_BINARY_DIR}/sax
)



