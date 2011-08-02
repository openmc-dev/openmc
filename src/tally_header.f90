module tally_header

  implicit none

!===============================================================================
! TALLYSCORE
!===============================================================================

  type TallyScore
     integer :: n_events
     real(8) :: val
     real(8) :: val_sq
  end type TallyScore

!===============================================================================
! TALLY
!===============================================================================

  type Tally
     integer :: uid
     integer :: type
     real(8) :: volume
     integer :: cell_type
     integer :: reaction_type
     integer :: material_type
     integer, allocatable :: reactions(:)
     integer, allocatable :: cells(:)
     integer, allocatable :: materials(:)
     integer, allocatable :: universes(:)
     real(8), allocatable :: energies(:)
     real(8) :: xyz_min(3)
     real(8) :: xyz_max(3)
     integer :: n_x
     integer :: n_y
     integer :: n_z
     type(TallyScore), allocatable :: score(:,:,:)
  end type Tally

end module tally_header
