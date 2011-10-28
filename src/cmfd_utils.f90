module cmfd_utils

  use global

  implicit none

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_input()

  end subroutine read_input

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(g,i,j,k,ng,nx,ny)

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: ng              ! max energy groups
    integer :: nx              ! maximum cells in x direction
    integer :: ny              ! maximum cells in y direction

    ! local variables
    integer :: nidx            ! index in matrix

    ! compute index
    nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

!===============================================================================
! PRINT_CMFD is a test routine to check if info from tally is being accessed 
!===============================================================================

  subroutine print_cmfd()

    write(7,*) 'hello world'

  end subroutine print_cmfd

end module cmfd_utils
