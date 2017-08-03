module mesh

  use algorithm,  only: binary_search
  use constants
  use global
  use mesh_header
  use message_passing

  implicit none

contains

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  pure subroutine get_mesh_bin(m, xyz, bin)
    type(RegularMesh), intent(in) :: m      ! mesh pointer
    real(8), intent(in)           :: xyz(:) ! coordinates
    integer, intent(out)          :: bin    ! tally bin

    integer :: n       ! size of mesh
    integer :: d       ! mesh dimension index
    integer :: ijk(3)  ! indices in mesh
    logical :: in_mesh ! was given coordinate in mesh at all?

    ! Get number of dimensions
    n = m % n_dimension

    ! Loop over the dimensions of the mesh
    do d = 1, n

      ! Check for cases where particle is outside of mesh
      if (xyz(d) < m % lower_left(d)) then
        bin = NO_BIN_FOUND
        return
      elseif (xyz(d) > m % upper_right(d)) then
        bin = NO_BIN_FOUND
        return
      end if
    end do

    ! Determine indices
    call get_mesh_indices(m, xyz, ijk, in_mesh)

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

  pure subroutine get_mesh_indices(m, xyz, ijk, in_mesh)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in)           :: xyz(:)  ! coordinates to check
    integer, intent(out)          :: ijk(:)  ! indices in mesh
    logical, intent(out)          :: in_mesh ! were given coords in mesh?

    ! Find particle in mesh
    ijk(:m % n_dimension) = ceiling((xyz(:m % n_dimension) - m % lower_left)/m % width)

    ! Determine if particle is in mesh
    if (any(ijk(:m % n_dimension) < 1) .or. &
         any(ijk(:m % n_dimension) > m % dimension)) then
      in_mesh = .false.
    else
      in_mesh = .true.
    end if

  end subroutine get_mesh_indices

!===============================================================================
! MESH_INDICES_TO_BIN maps (i), (i,j), or (i,j,k) indices to a single bin number
! for use in a TallyObject results array
!===============================================================================

  pure function mesh_indices_to_bin(m, ijk) result(bin)
    type(RegularMesh), intent(in) :: m
    integer, intent(in)           :: ijk(:)
    integer                       :: bin

    if (m % n_dimension == 1) then
      bin = ijk(1)
    elseif (m % n_dimension == 2) then
      bin = (ijk(2) - 1) * m % dimension(1) + ijk(1)
    elseif (m % n_dimension == 3) then
      bin = ((ijk(3) - 1) * m % dimension(2) + (ijk(2) - 1)) &
           * m % dimension(1) + ijk(1)
    end if

  end function mesh_indices_to_bin

!===============================================================================
! BIN_TO_MESH_INDICES maps a single mesh bin from a TallyObject results array to
! (i), (i,j), or (i,j,k) indices
!===============================================================================

  pure subroutine bin_to_mesh_indices(m, bin, ijk)
    type(RegularMesh), intent(in) :: m
    integer, intent(in)           :: bin
    integer, intent(out)          :: ijk(:)

    if (m % n_dimension == 1) then
      ijk(1) = bin
    else if (m % n_dimension == 2) then
      ijk(1) = mod(bin - 1, m % dimension(1)) + 1
      ijk(2) = (bin - 1)/m % dimension(1) + 1
    else if (m % n_dimension == 3) then
      ijk(1) = mod(bin - 1, m % dimension(1)) + 1
      ijk(2) = mod(bin - 1, m % dimension(1) * m % dimension(2)) &
           / m % dimension(1) + 1
      ijk(3) = (bin - 1)/(m % dimension(1) * m % dimension(2)) + 1
    end if

  end subroutine bin_to_mesh_indices

!===============================================================================
! COUNT_BANK_SITES determines the number of fission bank sites in each cell of a
! given mesh as well as an optional energy group structure. This can be used for
! a variety of purposes (Shannon entropy, CMFD, uniform fission source
! weighting)
!===============================================================================

  subroutine count_bank_sites(m, bank_array, cnt, energies, size_bank, &
       sites_outside)

    type(RegularMesh), pointer :: m             ! mesh to count sites
    type(Bank), intent(in)     :: bank_array(:) ! fission or source bank
    real(8),    intent(out)    :: cnt(:,:,:,:)  ! weight of sites in each
    ! cell and energy group
    real(8), intent(in),    optional :: energies(:)   ! energy grid to search
    integer(8), intent(in), optional :: size_bank     ! # of bank sites (on each proc)
    logical, intent(inout), optional :: sites_outside ! were there sites outside mesh?
    real(8), allocatable :: cnt_(:,:,:,:)

    integer :: i        ! loop index for local fission sites
    integer :: n_sites  ! size of bank array
    integer :: ijk(3)   ! indices on mesh
    integer :: n        ! number of energy groups / size
    integer :: e_bin    ! energy_bin
#ifdef MPI
    integer :: mpi_err  ! MPI error code
#endif
    logical :: in_mesh  ! was single site outside mesh?
    logical :: outside  ! was any site outside mesh?

    ! initialize variables
    allocate(cnt_(size(cnt,1), size(cnt,2), size(cnt,3), size(cnt,4)))
    cnt_ = ZERO
    outside = .false.

    ! Set size of bank
    if (present(size_bank)) then
      n_sites = int(size_bank,4)
    else
      n_sites = size(bank_array)
    end if

    ! Determine number of energies in group structure
    if (present(energies)) then
      n = size(energies) - 1
    else
      n = 1
    end if

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do i = 1, n_sites
      ! determine scoring bin for entropy mesh
      call get_mesh_indices(m, bank_array(i) % xyz, ijk, in_mesh)

      ! if outside mesh, skip particle
      if (.not. in_mesh) then
        outside = .true.
        cycle
      end if

      ! determine energy bin
      if (present(energies)) then
        if (bank_array(i) % E < energies(1)) then
          e_bin = 1
        elseif (bank_array(i) % E > energies(n + 1)) then
          e_bin = n
        else
          e_bin = binary_search(energies, n + 1, bank_array(i) % E)
        end if
      else
        e_bin = 1
      end if

      ! add to appropriate mesh box
      cnt_(e_bin,ijk(1),ijk(2),ijk(3)) = cnt_(e_bin,ijk(1),ijk(2),ijk(3)) + &
           bank_array(i) % wgt
    end do FISSION_SITES

#ifdef MPI
    ! collect values from all processors
    n = size(cnt_)
    call MPI_REDUCE(cnt_, cnt, n, MPI_REAL8, MPI_SUM, 0, mpi_intracomm, mpi_err)

    ! Check if there were sites outside the mesh for any processor
    if (present(sites_outside)) then
      call MPI_REDUCE(outside, sites_outside, 1, MPI_LOGICAL, MPI_LOR, 0, &
           mpi_intracomm, mpi_err)
    end if
#else
    sites_outside = outside
    cnt = cnt_
#endif

  end subroutine count_bank_sites

!===============================================================================
! MESH_INTERSECTS determines if a line between xyz0 and xyz1 intersects the
! outer boundary of the given mesh. This is important for determining whether a
! track will score to a mesh tally.
!===============================================================================

  pure function mesh_intersects_1d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(1)
    real(8), intent(in) :: xyz1(1)
    logical :: intersects

    real(8) :: x0   ! track start point
    real(8) :: x1   ! track end point
    real(8) :: xm0 ! lower-left coordinates of mesh
    real(8) :: xm1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)

    ! Copy coordinates of ending point
    x1 = xyz1(1)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      intersects = .true.
      return
    end if

    ! Check if line intersects right surface
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      intersects = .true.
      return
    end if

  end function mesh_intersects_1d

  pure function mesh_intersects_2d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(2)
    real(8), intent(in) :: xyz1(2)
    logical :: intersects

    real(8) :: x0, y0   ! track start point
    real(8) :: x1, y1   ! track end point
    real(8) :: xi, yi   ! track intersection point with mesh
    real(8) :: xm0, ym0 ! lower-left coordinates of mesh
    real(8) :: xm1, ym1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)
    y0 = xyz0(2)

    ! Copy coordinates of ending point
    x1 = xyz1(1)
    y1 = xyz1(2)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)
    ym0 = m % lower_left(2)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)
    ym1 = m % upper_right(2)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface -- calculate the intersection point
    ! y
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects back surface -- calculate the intersection point
    ! x
    if ((y0 < ym0 .and. y1 > ym0) .or. (y0 > ym0 .and. y1 < ym0)) then
      xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects right surface -- calculate the intersection
    ! point y
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects front surface -- calculate the intersection point
    ! x
    if ((y0 < ym1 .and. y1 > ym1) .or. (y0 > ym1 .and. y1 < ym1)) then
      xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1) then
        intersects = .true.
        return
      end if
    end if

  end function mesh_intersects_2d

  pure function mesh_intersects_3d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(3)
    real(8), intent(in) :: xyz1(3)
    logical :: intersects

    real(8) :: x0, y0, z0    ! track start point
    real(8) :: x1, y1, z1    ! track end point
    real(8) :: xi, yi, zi    ! track intersection point with mesh
    real(8) :: xm0, ym0, zm0 ! lower-left coordinates of mesh
    real(8) :: xm1, ym1, zm1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)
    y0 = xyz0(2)
    z0 = xyz0(3)

    ! Copy coordinates of ending point
    x1 = xyz1(1)
    y1 = xyz1(2)
    z1 = xyz1(3)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)
    ym0 = m % lower_left(2)
    zm0 = m % lower_left(3)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)
    ym1 = m % upper_right(2)
    zm1 = m % upper_right(3)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface -- calculate the intersection point
    ! (y,z)
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0)
      zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects back surface -- calculate the intersection point
    ! (x,z)
    if ((y0 < ym0 .and. y1 > ym0) .or. (y0 > ym0 .and. y1 < ym0)) then
      xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0)
      zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects bottom surface -- calculate the intersection
    ! point (x,y)
    if ((z0 < zm0 .and. z1 > zm0) .or. (z0 > zm0 .and. z1 < zm0)) then
      xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0)
      yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0)
      if (xi >= xm0 .and. xi < xm1 .and. yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects right surface -- calculate the intersection point
    ! (y,z)
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0)
      zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects front surface -- calculate the intersection point
    ! (x,z)
    if ((y0 < ym1 .and. y1 > ym1) .or. (y0 > ym1 .and. y1 < ym1)) then
      xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0)
      zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects top surface -- calculate the intersection point
    ! (x,y)
    if ((z0 < zm1 .and. z1 > zm1) .or. (z0 > zm1 .and. z1 < zm1)) then
      xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0)
      yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0)
      if (xi >= xm0 .and. xi < xm1 .and. yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

  end function mesh_intersects_3d

end module mesh
