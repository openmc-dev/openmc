module fox_m_fsys_realtypes

  implicit none
  private

  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)

  public :: sp
  public :: dp

end module fox_m_fsys_realtypes
