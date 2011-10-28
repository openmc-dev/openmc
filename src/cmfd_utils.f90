module cmfd_utils

  use global
  use mesh,          only: mesh_indices_to_bin
  use tally_header,  only: TallyObject, TallyScore

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

    integer :: bins(TALLY_TYPES)       ! bin for tally_types, for filters
    integer :: ijk(3)                  ! indices for mesh cell where tally is
    integer :: score_index             ! index in tally score to get value

    real(8) :: tally_val               ! value of tally being extracted    

    type(TallyObject), pointer :: t    ! pointer for a tally object
    type(StructuredMesh), pointer :: m ! pointer for mesh object

    ! associate pointers with objects
    t => tallies(1)
    m => meshes(t % mesh)

    ! set all bins to 1
    bins = 1

    ! get mesh indices, first we will first force to 1,1,1
    ijk = (/ 1, 1, 1 /)

    ! apply filters, here we will just try a mesh filter first
    bins(T_MESH) = mesh_indices_to_bin(m,ijk)

    ! calculate score index from bins
    score_index = sum((bins - 1) * t%stride) + 1

    ! get value from tally object
    tally_val = t%scores(score_index,1)%val

    ! write value to file
    write(7,*) "Tally value is:",tally_val

  end subroutine print_cmfd

end module cmfd_utils
