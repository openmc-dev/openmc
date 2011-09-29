module mesh

  use mesh_header

  implicit none

contains

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  subroutine get_mesh_bin(m, xyz, bin, in_mesh)

    type(StructuredMesh), pointer :: m
    real(8), intent(in)           :: xyz(:)
    integer, intent(out)          :: bin
    logical, intent(out)          :: in_mesh

    integer, allocatable :: ijk(:)

    ! Create indices array same size as xyz
    allocate(ijk(size(xyz)))
    
    ! Determine indices
    call get_mesh_indices(m, xyz, ijk, in_mesh)

    ! Convert indices to bin
    bin = mesh_indices_to_bin(m, ijk)

  end subroutine get_mesh_bin

!===============================================================================
! GET_MESH_INDICES determines the indices of a particle in a structured mesh
!===============================================================================

  subroutine get_mesh_indices(m, xyz, ijk, in_mesh)

    type(StructuredMesh), pointer :: m
    real(8), intent(in)           :: xyz(:)
    integer, intent(out)          :: ijk(:)
    logical, intent(out)          :: in_mesh

    ! Find particle in mesh
    ijk = ceiling((xyz - m % origin)/m % width)

    ! Determine if particle is in mesh
    if (any(ijk < 1) .or. any(ijk > m % dimension)) then
       in_mesh = .false.
    else
       in_mesh = .true.
    end if

  end subroutine get_mesh_indices

!===============================================================================
! MESH_INDICES_TO_BIN maps (i,j) or (i,j,k) indices to a single bin number for
! use in a TallyObject scores array
!===============================================================================

  function mesh_indices_to_bin(m, ijk) result(bin)

    type(StructuredMesh), pointer :: m
    integer, intent(in)           :: ijk(:)
    integer                       :: bin

    integer :: n_x
    integer :: n_y
    integer :: n_z

    n_x = m % dimension(1)
    n_y = m % dimension(2)

    if (m % n_dimension == 2) then
       bin = (ijk(1) - 1)*n_y + ijk(2)
    elseif (m % n_dimension == 3) then
       n_z = m % dimension(3)
       bin = (ijk(1) - 1)*n_y*n_z + (ijk(2) - 1)*n_z + ijk(3)
    end if

  end function mesh_indices_to_bin

!===============================================================================
! BIN_TO_MESH_INDICES maps a single mesh bin from a TallyObject scores array to
! (i,j) or (i,j,k) indices
!===============================================================================

  subroutine bin_to_mesh_indices(m, bin, ijk)

    type(StructuredMesh), pointer :: m
    integer, intent(in)           :: bin
    integer, intent(out)          :: ijk(:)

  end subroutine bin_to_mesh_indices

end module mesh
