module ace_header

  use constants,   only: MAX_FILE_LEN, ZERO
  use dict_header, only: DictIntInt
  use endf_header, only: Tab1
  use stl_vector,  only: VectorInt

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
    type(Tab1), allocatable :: multiplicity_E ! Energy-dependent neutron yield
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
  end type UrrData

  contains

!===============================================================================
! DISTENERGY_CLEAR resets and deallocates data in DistEnergy.
!===============================================================================

    recursive subroutine distenergy_clear(this)

      class(DistEnergy), intent(inout) :: this ! The DistEnergy object to clear

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

      if (associated(this % edist)) then
        call this % edist % clear()
        deallocate(this % edist)
      end if

    end subroutine reaction_clear

end module ace_header
