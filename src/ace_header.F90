module ace_header

  use constants,   only: MAX_FILE_LEN, ZERO
  use dict_header, only: DictIntInt
  use endf_header, only: Tab1
  use secondary_header, only: SecondaryDistribution, AngleEnergyContainer
  use stl_vector,  only: VectorInt

  implicit none

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
    type(SecondaryDistribution) :: secondary

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
! REACTION_CLEAR resets and deallocates data in Reaction.
!===============================================================================

    subroutine reaction_clear(this)
      class(Reaction), intent(inout) :: this ! The Reaction object to clear

      if (associated(this % multiplicity_E)) deallocate(this % multiplicity_E)
    end subroutine reaction_clear

end module ace_header
