module mesh

  use algorithm,  only: binary_search
  use bank_header, only: bank
  use constants
  use mesh_header
  use message_passing

  implicit none

contains

!===============================================================================
! COUNT_BANK_SITES determines the number of fission bank sites in each cell of a
! given mesh as well as an optional energy group structure. This can be used for
! a variety of purposes (Shannon entropy, CMFD, uniform fission source
! weighting)
!===============================================================================

  subroutine count_bank_sites(m, bank_array, cnt, energies, size_bank, &
       sites_outside)

    type(RegularMesh), intent(in) :: m             ! mesh to count sites
    type(Bank), intent(in)     :: bank_array(:) ! fission or source bank
    real(8),    intent(out)    :: cnt(:,:)      ! weight of sites in each
    ! cell and energy group
    real(8), intent(in),    optional :: energies(:)   ! energy grid to search
    integer(8), intent(in), optional :: size_bank     ! # of bank sites (on each proc)
    logical, intent(inout), optional :: sites_outside ! were there sites outside mesh?
    real(8), allocatable :: cnt_(:,:)

    integer :: i        ! loop index for local fission sites
    integer :: n_sites  ! size of bank array
    integer :: n        ! number of energy groups / size
    integer :: mesh_bin ! mesh bin
    integer :: e_bin    ! energy bin
#ifdef MPI
    integer :: mpi_err  ! MPI error code
#endif
    logical :: outside  ! was any site outside mesh?

    ! initialize variables
    allocate(cnt_(size(cnt,1), size(cnt,2)))
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
      call m % get_bin(bank_array(i) % xyz, mesh_bin)

      ! if outside mesh, skip particle
      if (mesh_bin == NO_BIN_FOUND) then
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
      cnt_(e_bin, mesh_bin) = cnt_(e_bin, mesh_bin) + bank_array(i) % wgt
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

end module mesh
