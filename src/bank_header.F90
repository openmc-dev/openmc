module bank_header

  use, intrinsic :: ISO_C_BINDING

  use error, only: E_ALLOCATE, set_errmsg

  implicit none

!===============================================================================
! BANK is used for storing fission sites in eigenvalue calculations. Since all
! the state information of a neutron is not needed, this type allows sites to be
! stored with less memory
!===============================================================================

  type, bind(C) :: Bank
    real(C_DOUBLE) :: wgt           ! weight of bank site
    real(C_DOUBLE) :: xyz(3)        ! location of bank particle
    real(C_DOUBLE) :: uvw(3)        ! diretional cosines
    real(C_DOUBLE) :: E             ! energy / energy group if in MG mode.
    integer(C_INT) :: delayed_group ! delayed group
  end type Bank

  ! Source and fission bank
  type(Bank), allocatable, target :: source_bank(:)
  type(Bank), allocatable, target :: fission_bank(:)
#ifdef _OPENMP
  type(Bank), allocatable, target :: master_fission_bank(:)
#endif

  integer(8) :: n_bank       ! # of sites in fission bank

!$omp threadprivate(fission_bank, n_bank)

contains

!===============================================================================
! FREE_MEMORY_BANK deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_bank()

    ! Deallocate fission and source bank and entropy
!$omp parallel
    if (allocated(fission_bank)) deallocate(fission_bank)
!$omp end parallel
#ifdef _OPENMP
    if (allocated(master_fission_bank)) deallocate(master_fission_bank)
#endif
    if (allocated(source_bank)) deallocate(source_bank)

  end subroutine free_memory_bank

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_source_bank(ptr, n) result(err) bind(C)
    ! Return a pointer to the source bank
    type(C_PTR), intent(out) :: ptr
    integer(C_INT64_T), intent(out) :: n
    integer(C_INT) :: err

    if (.not. allocated(source_bank)) then
      err = E_ALLOCATE
      call set_errmsg("Source bank has not been allocated.")
    else
      err = 0
      ptr = C_LOC(source_bank)
      n = size(source_bank)
    end if
  end function openmc_source_bank

end module bank_header
