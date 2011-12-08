module source_header

  implicit none

!===============================================================================
! EXTSOURCE describes an external source of neutrons for a fixed-source problem
! or for the starting source in a criticality problem
!===============================================================================

  type ExtSource
     integer :: type                    ! type of source, e.g. 'box' or 'cell'
     real(8), allocatable :: values(:)  ! values for particular source type
  end type ExtSource

end module source_header
