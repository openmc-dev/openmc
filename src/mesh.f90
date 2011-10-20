module mesh

  use constants
  use mesh_header
  use particle_header, only: Particle

  implicit none

contains

!===============================================================================
! SURFACE_CROSSINGS determines which surfaces are crossed in a mesh after a
! collision and calls the appropriate subroutine to tally surface currents
!===============================================================================

  subroutine surface_crossings(m, p)

    type(StructuredMesh), pointer :: m
    type(Particle),       pointer :: p

    integer :: i             ! loop indices
    integer :: j             ! loop indices
    integer :: ijk0(3)       ! indices of starting coordinates
    integer :: ijk1(3)       ! indices of ending coordinates
    integer :: n_cross       ! number of surface crossings
    integer :: surface       ! surface/direction of crossing, e.g. IN_RIGHT
    real(8) :: uvw(3)        ! cosine of angle of particle
    real(8) :: xyz0(3)       ! starting/intermediate coordinates
    real(8) :: xyz1(3)       ! ending coordinates of particle
    real(8) :: xyz_cross(3)  ! coordinates of bounding surfaces
    real(8) :: d(3)          ! distance to each bounding surface
    real(8) :: distance      ! actual distance traveled
    logical :: start_in_mesh ! particle's starting xyz in mesh?
    logical :: end_in_mesh   ! particle's ending xyz in mesh?
    logical :: x_same        ! same starting/ending x index (i)
    logical :: y_same        ! same starting/ending y index (j)
    logical :: z_same        ! same starting/ending z index (k)

    ! Copy starting and ending location of particle
    xyz0 = p % last_xyz
    xyz1 = p % xyz

    ! Determine indices for starting and ending location
    call get_mesh_indices(m, xyz0, ijk0, start_in_mesh)
    write (*,'(A,1X,3I5)') 'Particle starts in:', ijk0

    call get_mesh_indices(m, xyz1, ijk1, end_in_mesh)
    write (*,'(A,1X,3I5)') 'Particle ends in:', ijk1

    ! Check to make sure start or end is in mesh
    if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) return

    ! Calculate number of surface crossings
    n_cross = sum(abs(ijk1 - ijk0))
    write (*,'(A,1X,I5)') 'Number of surface crossings:', n_cross
    if (n_cross == 0) return

    ! Copy particle's direction
    uvw = p % uvw

    ! ==========================================================================
    ! SPECIAL CASES WHERE TWO INDICES ARE THE SAME

    x_same = (ijk0(1) == ijk1(1))
    y_same = (ijk0(2) == ijk1(2))
    z_same = (ijk0(3) == ijk1(3))

    if (x_same .and. y_same) then
       ! Only z crossings
       print *, "Only z crossings"
       if (uvw(3) > 0) then
          do i = ijk0(3), ijk1(3) - 1
             ijk0(3) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_TOP"
          end do
       else
          do i = ijk0(3) - 1, ijk1(3), -1
             ijk0(3) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_TOP"
          end do
       end if
       return
    elseif (x_same .and. z_same) then
       ! Only y crossings
       print *, "Only y crossings"
       if (uvw(2) > 0) then
          do i = ijk0(2), ijk1(2) - 1
             ijk0(2) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_TOP"
          end do
       else
          do i = ijk0(2) - 1, ijk1(2), -1
             ijk0(2) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_TOP"
          end do
       end if
       return
    elseif (y_same .and. z_same) then
       ! Only x crossings
       print *, "Only x crossings"
       if (uvw(1) > 0) then
          do i = ijk0(1), ijk1(1) - 1
             ijk0(1) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_TOP"
          end do
       else
          do i = ijk0(1) - 1, ijk1(1), -1
             ijk0(1) = i
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_TOP"
          end do
       end if
       return
    end if

    ! ==========================================================================
    ! GENERIC CASE

    ! Bounding coordinates
    do i = 1, 3
       if (uvw(i) > 0) then
          xyz_cross(i) = m % origin(i) + ijk0(i) * m % width(i)
       else
          xyz_cross(i) = m % origin(i) + (ijk0(i) - 1) * m % width(i)
       end if
    end do

    do i = 1, n_cross
       ! Calculate distance to each bounding surface. We need to treat special
       ! case where the cosine of the angle is zero since this would result in a
       ! divide-by-zero.

       do j = 1, 3
          if (uvw(j) == 0) then
             d(j) = INFINITY
          else
             d(j) = (xyz_cross(j) - xyz0(j))/uvw(j)
          end if
       end do

       ! Determine the closest bounding surface of the mesh cell by calculating
       ! the minimum distance

       distance = minval(d)

       ! Now use the minimum distance and diretion of the particle to determine
       ! which surface was crossed

       if (distance == d(1)) then
          if (uvw(1) > 0) then
             ! Crossing into right mesh cell -- this is treated as outgoing
             ! current from (i,j,k)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_RIGHT"
             ijk0(1) = ijk0(1) + 1
             xyz_cross(1) = xyz_cross(1) + m % width(1)
          else
             ! Crossing into left mesh cell -- this is treated as incoming
             ! current in (i-1,j,k)
             ijk0(1) = ijk0(1) - 1
             xyz_cross(1) = xyz_cross(1) - m % width(1)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_RIGHT"

          end if
       elseif (distance == d(2)) then
          if (uvw(2) > 0) then
             ! Crossing into front mesh cell -- this is treated as outgoing
             ! current in (i,j,k)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_FRONT"
             ijk0(2) = ijk0(2) + 1
             xyz_cross(2) = xyz_cross(2) + m % width(2)
          else
             ! Crossing into back mesh cell -- this is treated as incoming
             ! current in (i,j-1,k)
             ijk0(2) = ijk0(2) - 1
             xyz_cross(2) = xyz_cross(2) - m % width(2)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_FRONT"
          end if
       else if (distance == d(3)) then
          if (uvw(3) > 0) then
             ! Crossing into top mesh cell -- this is treated as outgoing
             ! current in (i,j,k)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "OUT_TOP"
             ijk0(3) = ijk0(3) + 1
             xyz_cross(3) = xyz_cross(3) + m % width(3)
          else
             ! Crossing into bottom mesh cell -- this is treated as incoming
             ! current in (i,j,k-1)
             ijk0(3) = ijk0(3) - 1
             xyz_cross(3) = xyz_cross(3) - m % width(3)
             if (all(ijk0 >= 0) .and. all(ijk0 <= m % dimension)) print *, ijk0, "IN_TOP"
          end if
       end if

       ! Calculate new coordinates
       xyz0 = xyz0 + distance * uvw
    end do

  end subroutine surface_crossings

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  subroutine get_mesh_bin(m, xyz, bin, in_mesh)

    type(StructuredMesh), pointer :: m
    real(8), intent(in)           :: xyz(:)
    integer, intent(out)          :: bin
    logical, intent(out)          :: in_mesh

    integer              :: n
    integer, allocatable :: ijk(:)

    ! Get number of dimensions
    n = m % n_dimension

    ! Create indices array same size as xyz
    allocate(ijk(n))
    
    ! Determine indices
    call get_mesh_indices(m, xyz(1:n), ijk, in_mesh)

    ! Convert indices to bin
    bin = mesh_indices_to_bin(m, ijk)

    ! Release memory for ijk
    deallocate(ijk)

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

    integer :: n_y
    integer :: n_z

    if (m % n_dimension == 2) then
       n_y = m % dimension(2)

       ijk(1) = (bin - 1)/n_y + 1
       ijk(2) = mod(bin - 1, n_y) + 1
    else if (m % n_dimension == 3) then
       n_y = m % dimension(2)
       n_z = m % dimension(3)

       ijk(1) = (bin - 1)/(n_y*n_z) + 1
       ijk(2) = mod(bin - 1, n_y*n_z)/n_z + 1
       ijk(3) = mod(bin - 1, n_z) + 1
    end if

  end subroutine bin_to_mesh_indices

end module mesh
