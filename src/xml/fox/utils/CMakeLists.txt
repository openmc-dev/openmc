set(fox_utils_srcs
  FoX_utils.f90
)

set(fox_utils_srcs_pp
  fox_m_utils_mtprng.F90
  fox_m_utils_uri.F90
  fox_m_utils_uuid.F90
)

add_library(
  fox_utils STATIC
  ${fox_utils_srcs}
  ${fox_utils_srcs_pp}
)

target_link_libraries(
  fox_utils
  fox_fsys
)

include_directories(
  ${fox_BINARY_DIR}/fsys
)

