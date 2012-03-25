module mesh

  use constants
  use global
  use mesh_header
  use particle_header, only: Particle

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  subroutine get_mesh_bin(m, xyz, bin)

    type(StructuredMesh), pointer :: m      ! mesh pointer
    real(8), intent(in)           :: xyz(:) ! coordinates
    integer, intent(out)          :: bin    ! tally bin

    integer :: n       ! size of mesh (2 or 3)
    integer :: ijk(3)  ! indices in mesh
    logical :: in_mesh ! was given coordinate in mesh at all?

    ! Get number of dimensions
    n = m % n_dimension

    ! Check for cases where particle is outside of mesh
    if (xyz(1) < m % lower_left(1)) then
       bin = NO_BIN_FOUND
       return
    elseif (xyz(1) > m % upper_right(1)) then
       bin = NO_BIN_FOUND
       return
    elseif (xyz(2) < m % lower_left(2)) then
       bin = NO_BIN_FOUND
       return
    elseif (xyz(2) > m % upper_right(2)) then
       bin = NO_BIN_FOUND
       return
    end if
    if (n > 2) then
       if (xyz(3) < m % lower_left(3)) then
          bin = NO_BIN_FOUND
          return
       elseif (xyz(3) > m % upper_right(3)) then
          bin = NO_BIN_FOUND
          return
       end if
    end if

    ! Determine indices
    call get_mesh_indices(m, xyz(1:n), ijk(1:n), in_mesh)

    ! Convert indices to bin
    if (in_mesh) then
       bin = mesh_indices_to_bin(m, ijk)
    else
       bin = NO_BIN_FOUND
    end if

  end subroutine get_mesh_bin

!===============================================================================
! GET_MESH_INDICES determines the indices of a particle in a structured mesh
!===============================================================================

  subroutine get_mesh_indices(m, xyz, ijk, in_mesh)

    type(StructuredMesh), pointer :: m
    real(8), intent(in)           :: xyz(:)  ! coordinates to check
    integer, intent(out)          :: ijk(:)  ! indices in mesh
    logical, intent(out)          :: in_mesh ! were given coords in mesh?

    ! Find particle in mesh
    ijk = ceiling((xyz - m % lower_left)/m % width)

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

    integer :: n_y ! number of mesh cells in y direction
    integer :: n_z ! number of mesh cells in z direction

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

    integer :: n_y ! number of mesh cells in y direction
    integer :: n_z ! number of mesh cells in z direction

    n_y = m % dimension(2)

    if (m % n_dimension == 2) then
       ijk(1) = (bin - 1)/n_y + 1
       ijk(2) = mod(bin - 1, n_y) + 1
    else if (m % n_dimension == 3) then
       n_z = m % dimension(3)
       ijk(1) = (bin - 1)/(n_y*n_z) + 1
       ijk(2) = mod(bin - 1, n_y*n_z)/n_z + 1
       ijk(3) = mod(bin - 1, n_z) + 1
    end if

  end subroutine bin_to_mesh_indices

!===============================================================================
! COUNT_FISSION_SITES determines the number of fission bank sites in each cell
! of a given mesh. This can be used for a variety of purposes (Shannon entropy,
! CMFD, uniform fission source weighting)
!===============================================================================

  subroutine count_fission_sites(m, count, total, sites_outside)

    type(StructuredMesh), pointer :: m             ! mesh to count sites
    real(8), intent(out)          :: count(:,:,:)  ! weight of sites in each cell
    real(8), intent(out)          :: total         ! total weight of sites
    logical, optional             :: sites_outside ! were there sites outside mesh?

    integer :: i       ! loop index for local fission sites
    integer :: ijk(3)  ! indices on mesh
    real(8) :: weight  ! accumulated weight of sites
    logical :: in_mesh ! was single site outside mesh?
    logical :: outside ! was any site outside mesh?
#ifdef MPI
    integer :: n       ! total size of count variable
#endif

    ! initialize variables
    count = ZERO
    weight = ZERO
    outside = .false.

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do i = 1, int(n_bank,4)
       ! determine scoring bin for entropy mesh
       call get_mesh_indices(m, fission_bank(i) % xyz, ijk, in_mesh)

       ! if outside mesh, skip particle
       if (.not. in_mesh) then
          outside = .true.
          cycle
       end if

       ! add weight
       weight = weight + ONE

       ! add to appropriate mesh box
       ! TODO: if tracking weight through bank, add weight instead
       count(ijk(1),ijk(2),ijk(3)) = count(ijk(1),ijk(2),ijk(3)) + 1
    end do FISSION_SITES

#ifdef MPI
    ! determine total number of mesh cells
    n = size(count,1) * size(count,2) * size(count,3)

    ! collect values from all processors
    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, count, n, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
    else
       call MPI_REDUCE(count, count, n, MPI_REAL8, MPI_SUM, 0, &
            MPI_COMM_WORLD, mpi_err)
    end if

    ! Check if there were sites outside the mesh for any processor
    if (present(sites_outside)) then
       call MPI_REDUCE(outside, sites_outside, 1, MPI_LOGICAL, MPI_LOR, 0, &
            MPI_COMM_WORLD, mpi_err)
    end if

    ! determine total weight of bank sites
    call MPI_REDUCE(weight, total, 1, MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, mpi_err)
#else
    total = weight
    sites_outside = outside
#endif

  end subroutine count_fission_sites

!===============================================================================
! MESH_INTERSECT
!===============================================================================

  function mesh_intersects(m, xyz0, xyz1) result(intersects)

    type(StructuredMesh), pointer :: m
    real(8), intent(in) :: xyz0(3)
    real(8), intent(in) :: xyz1(3)
    logical :: intersects

    real(8) :: x0, y0, z0 ! track start point
    real(8) :: x1, y1, z1 ! track end point
    real(8) :: xi, yi, zi ! track intersection point with mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)
    y0 = xyz0(2)
    z0 = xyz0(3)

    ! Copy coordinates of ending point
    x1 = xyz1(1)
    y1 = xyz1(2)
    z1 = xyz1(3)

    ! Check if line intersects bottom surface -- calculate the intersection
    ! point (y,z)
    xi = m % lower_left(1)
    yi = y0 + (xi - x0) * (y1 - y0) / (x1 - x0)
    zi = z0 + (xi - x0) * (z1 - z0) / (x1 - x0)
    if (yi >= y0 .and. yi < y1 .and. zi >= z0 .and. zi < z1) then
       intersects = .true.
       return
    end if
    
    ! Check if line intersects left surface -- calculate the intersection point
    ! (x,z)
    yi = m % lower_left(2)
    xi = x0 + (yi - y0) * (x1 - x0) / (y1 - y0)
    zi = z0 + (yi - y0) * (z1 - z0) / (y1 - y0)
    if (xi >= x0 .and. xi < x1 .and. zi >= z0 .and. zi < z1) then
       intersects = .true.
       return
    end if
    
    ! Check if line intersects front surface -- calculate the intersection point
    ! (x,y)
    zi = m % lower_left(3)
    xi = x0 + (zi - z0) * (x1 - x0) / (z1 - z0)
    yi = y0 + (zi - z0) * (y1 - y0) / (z1 - z0)
    if (xi >= x0 .and. xi < x1 .and. yi >= y0 .and. yi < y1) then
       intersects = .true.
       return
    end if
    
    ! Check if line intersects top surface -- calculate the intersection
    ! point (y,z)
    xi = m % upper_right(1)
    yi = y0 + (xi - x0) * (y1 - y0) / (x1 - x0)
    zi = z0 + (xi - x0) * (z1 - z0) / (x1 - x0)
    if (yi >= y0 .and. yi < y1 .and. zi >= z0 .and. zi < z1) then
       intersects = .true.
       return
    end if
    
    ! Check if line intersects right surface -- calculate the intersection point
    ! (x,z)
    yi = m % upper_right(2)
    xi = x0 + (yi - y0) * (x1 - x0) / (y1 - y0)
    zi = z0 + (yi - y0) * (z1 - z0) / (y1 - y0)
    if (xi >= x0 .and. xi < x1 .and. zi >= z0 .and. zi < z1) then
       intersects = .true.
       return
    end if
    
    ! Check if line intersects back surface -- calculate the intersection point
    ! (x,y)
    zi = m % upper_right(3)
    xi = x0 + (zi - z0) * (x1 - x0) / (z1 - z0)
    yi = y0 + (zi - z0) * (y1 - y0) / (z1 - z0)
    if (xi >= x0 .and. xi < x1 .and. yi >= y0 .and. yi < y1) then
       intersects = .true.
       return
    end if

  end function mesh_intersects

end module mesh
