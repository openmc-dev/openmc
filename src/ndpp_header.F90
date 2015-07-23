module ndpp_header

  use constants, only: ZERO

  implicit none

!===============================================================================
! GROUPTRANSFER contains probability tables for the unresolved resonance range.
!===============================================================================

  type GroupTransfer
    real(8), allocatable :: outgoing(:,:) ! Outgoing transfer probabilities
                                          ! Dimension of (moments, gmin:gmax)
  end type GroupTransfer

!===============================================================================
! NDPP contains all the data produced by NDPP for a particular nuclide, s(a,b),
! or material.
!===============================================================================

  type Ndpp
    integer :: zaid    = 0         ! ZAID identifier = 1000*Z + A
    real(8) :: kT      = ZERO      ! Boltzmann constant * temperature (MeV)
    logical :: is_nuc  = .False.   ! Is our data a nuc or an sab?
    logical :: is_init = .False.   ! Is object initialized

    ! Elastic data is allocatable since it will be unique for each Ndpp object
    real(8), allocatable :: el_Ein(:)         ! Elastic Ein grid
    integer, allocatable :: el_Ein_srch(:)    ! Elastic Ein grid search bounds
    type(GroupTransfer), allocatable :: el(:) ! Elastic Data, Dim is # of Ein
    ! Inelastic data are pointers since it may not be unique and thus
    ! may point to data in another object.
    ! This is not a good design practice to have one object point to the member
    ! data of another object, but will do for now until the NDPP library format
    ! is modified such that all temperatures are written together, at which point
    ! the Ndpp object can contain all temperatures within one class.
    ! Doing so will require research about how to perform temperature
    ! interpolation and thus can be a bit in to the future.

    ! Inelastic Ein grid
    real(8), pointer :: inel_Ein(:) => null()
    ! Inelastic Ein grid search bounds
    integer, pointer :: inel_Ein_srch(:) => null()
    ! Inelastic Data, Dimension is # of Ein
    type(GroupTransfer), pointer :: inel(:) => null()
    ! Inelastic Data, Dimension is # of Ein
    type(GroupTransfer), pointer :: nuinel(:) => null()

    ! Chi data is temperature dependent, so it is allocatable
    real(8), allocatable :: chi_Ein(:)   ! Ein grid for all chi
    real(8), allocatable :: chi(:,:)     ! Data grid for ndpp chi data
                                         ! dimensions of chi: (g, Ein)
    real(8), allocatable :: chi_p(:,:)   ! Same for prompt only
    real(8), allocatable :: chi_d(:,:,:) ! Same, but additional dimension
                                         ! for precursor group

! Type-Bound procedures
contains
      procedure :: clear => ndpp_clear ! Clears NDPP structure
  end type Ndpp

  contains

!===============================================================================
! CLASS METHODS
!===============================================================================

!===============================================================================
! NDPP_CLEAR clears the object, of course.
!===============================================================================

  subroutine ndpp_clear(this)
    class(Ndpp), intent(inout) :: this ! Ndpp object to act on

    if (allocated(this % el)) then
      deallocate(this % el)
    end if
    if (associated(this % inel)) then
      deallocate(this % inel)
    end if
    if (associated(this % nuinel)) then
      deallocate(this % nuinel)
    end if
    if (allocated(this % el_Ein)) then
      deallocate(this % el_Ein)
    end if
    if (allocated(this % el_Ein_srch)) then
      deallocate(this % el_Ein_srch)
    end if
    if (associated(this % inel_Ein)) then
      deallocate(this % inel_Ein)
    end if
    if (associated(this % inel_Ein_srch)) then
      deallocate(this % inel_Ein_srch)
    end if
    if (allocated(this % chi)) then
      deallocate(this % chi)
    end if
    if (allocated(this % chi_p)) then
      deallocate(this % chi_p)
    end if
    if (allocated(this % chi_d)) then
      deallocate(this % chi_d)
    end if
    if (allocated(this % chi_Ein)) then
      deallocate(this % chi_Ein)
    end if

    this % zaid    = 0
    this % kT      = ZERO
    this % is_nuc  = .False.
    this % is_init = .False.

  end subroutine ndpp_clear

end module ndpp_header