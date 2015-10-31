module ace_header

  use constants,     only: MAX_FILE_LEN, ZERO
  use endf_header,   only: Tab1
  use list_header,   only: ListInt

  implicit none

!===============================================================================
! DISTANGLE contains data for a tabular secondary angle distribution whether it
! be tabular or 32 equiprobable cosine bins
!===============================================================================

  type DistAngle
    integer              :: n_energy    ! # of incoming energies
    real(8), allocatable :: energy(:)   ! incoming energy grid
    integer, allocatable :: type(:)     ! type of distribution
    integer, allocatable :: location(:) ! location of each table
    real(8), allocatable :: data(:)     ! angular distribution data

    ! Type-Bound procedures
    contains
      procedure :: clear => distangle_clear ! Deallocates DistAngle
  end type DistAngle

!===============================================================================
! DISTENERGY contains data for a secondary energy distribution for all
! scattering laws
!===============================================================================

  type DistEnergy
    integer    :: law                 ! secondary distribution law
    type(Tab1) :: p_valid             ! probability of law validity
    real(8), allocatable :: data(:)   ! energy distribution data

    ! For reactions that may have multiple energy distributions such as (n,2n),
    ! this pointer allows multiple laws to be stored
    type(DistEnergy), pointer :: next => null()

    ! Type-Bound procedures
    contains
      procedure :: clear => distenergy_clear ! Deallocates DistEnergy
  end type DistEnergy

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type Reaction
    integer :: MT                      ! ENDF MT value
    real(8) :: Q_value                 ! Reaction Q value
    integer :: multiplicity            ! Number of secondary particles released
    type(Tab1), pointer :: multiplicity_E => null() ! Energy-dependent neutron yield
    integer :: threshold               ! Energy grid index of threshold
    logical :: scatter_in_cm           ! scattering system in center-of-mass?
    logical :: multiplicity_with_E = .false. ! Flag to indicate E-dependent multiplicity
    real(8), allocatable :: sigma(:)   ! Cross section values
    logical :: has_angle_dist          ! Angle distribution present?
    logical :: has_energy_dist         ! Energy distribution present?
    type(DistAngle)           :: adist ! Secondary angular distribution
    type(DistEnergy), pointer :: edist => null() ! Secondary energy distribution

    ! Type-Bound procedures
    contains
      procedure :: clear => reaction_clear ! Deallocates Reaction
  end type Reaction

!===============================================================================
! URRDATA contains probability tables for the unresolved resonance range.
!===============================================================================

  type UrrData
    integer :: n_energy        ! # of incident neutron energies
    integer :: n_prob          ! # of probabilities
    integer :: interp          ! inteprolation (2=lin-lin, 5=log-log)
    integer :: inelastic_flag  ! inelastic competition flag
    integer :: absorption_flag ! other absorption flag
    logical :: multiply_smooth ! multiply by smooth cross section?
    real(8), allocatable :: energy(:)   ! incident energies
    real(8), allocatable :: prob(:,:,:) ! actual probabibility tables

    ! Type-Bound procedures
    contains
      procedure :: clear => urrdata_clear ! Deallocates UrrData
  end type UrrData

  contains

!===============================================================================
! DISTANGLE_CLEAR resets and deallocates data in Reaction.
!===============================================================================

    subroutine distangle_clear(this)

      class(DistAngle), intent(inout) :: this ! The DistAngle object to clear

      if (allocated(this % energy)) &
           deallocate(this % energy, this % type, this % location, this % data)

    end subroutine distangle_clear

!===============================================================================
! DISTENERGY_CLEAR resets and deallocates data in DistEnergy.
!===============================================================================

    recursive subroutine distenergy_clear(this)

      class(DistEnergy), intent(inout) :: this ! The DistEnergy object to clear

      ! Clear p_valid
      call this % p_valid % clear()

      if (allocated(this % data)) &
           deallocate(this % data)

      if (associated(this % next)) then
        ! recursively clear this item
        call this % next % clear()
        deallocate(this % next)
      end if

    end subroutine distenergy_clear

!===============================================================================
! REACTION_CLEAR resets and deallocates data in Reaction.
!===============================================================================

    subroutine reaction_clear(this)

      class(Reaction), intent(inout) :: this ! The Reaction object to clear

      if (allocated(this % sigma)) deallocate(this % sigma)

      if (associated(this % multiplicity_E)) deallocate(this % multiplicity_E)

      if (associated(this % edist)) then
        call this % edist % clear()
        deallocate(this % edist)
      end if

      call this % adist % clear()

    end subroutine reaction_clear

!===============================================================================
! URRDATA_CLEAR resets and deallocates data in Reaction.
!===============================================================================

    subroutine urrdata_clear(this)

      class(UrrData), intent(inout) :: this ! The UrrData object to clear

      if (allocated(this % energy)) &
           deallocate(this % energy, this % prob)

    end subroutine urrdata_clear


end module ace_header
