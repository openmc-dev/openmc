module source_header

  implicit none

!===============================================================================
! EXTSOURCE describes an external source of neutrons for a fixed-source problem
! or for the starting source in a criticality problem
!===============================================================================

  type ExtSource
     integer :: type_space              ! spacial distributione, e.g. 'box' or 'point'
     integer :: type_angle              ! angle distribution, e.g. 'isotropic'
     integer :: type_energy             ! energy distribution, e.g. 'Watt'
     real(8), allocatable :: values(:)  ! values for particular source type
  end type ExtSource

end module source_header
