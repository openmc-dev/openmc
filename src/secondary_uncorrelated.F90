module secondary_uncorrelated

  use hdf5, only: HID_T

  use angle_distribution, only: AngleDistribution
  use angleenergy_header, only: AngleEnergy
  use constants, only: ONE, TWO, MAX_WORD_LEN
  use energy_distribution, only: EnergyDistribution, LevelInelastic, &
       ContinuousTabular, MaxwellEnergy, Evaporation, WattEnergy, DiscretePhoton
  use error, only: warning
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       object_exists
  use random_lcg, only: prn

!===============================================================================
! UNCORRELATEDANGLEENERGY represents an uncorrelated angle-energy
! distribution. This corresponds to when an energy distribution is given in ENDF
! File 5/6 and an angular distribution is given in ENDF File 4.
!===============================================================================

  type, extends(AngleEnergy) :: UncorrelatedAngleEnergy
    logical :: fission = .false.
    type(AngleDistribution) :: angle
    class(EnergyDistribution), allocatable :: energy
  contains
    procedure :: sample => uncorrelated_sample
    procedure :: from_hdf5 => uncorrelated_from_hdf5
  end type UncorrelatedAngleEnergy

contains

  subroutine uncorrelated_sample(this, E_in, E_out, mu)
    class(UncorrelatedAngleEnergy), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sampled scattering cosine

    ! Sample cosine of scattering angle
    if (this%fission) then
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! For fission, the angle is not used, so just assign a dummy value
      mu = ONE
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    elseif (allocated(this%angle%energy)) then
      mu = this%angle%sample(E_in)
    else
      ! no angle distribution given => assume isotropic for all energies
      mu = TWO*prn() - ONE
    end if

    ! Sample outgoing energy
    E_out = this%energy%sample(E_in)
  end subroutine uncorrelated_sample

  subroutine uncorrelated_from_hdf5(this, group_id)
    class(UncorrelatedAngleEnergy), intent(inout) :: this
    integer(HID_T),                 intent(in)    :: group_id

    integer(HID_T) :: energy_group
    integer(HID_T) :: angle_group
    character(MAX_WORD_LEN) :: type

    ! Check if angle group is present & read
    if (object_exists(group_id, 'angle')) then
      angle_group = open_group(group_id, 'angle')
      call this%angle%from_hdf5(angle_group)
      call close_group(angle_group)
    end if

    ! Check if energy group is present & read
    if (object_exists(group_id, 'energy')) then
      energy_group = open_group(group_id, 'energy')
      call read_attribute(type, energy_group, 'type')
      select case (type)
      case ('discrete_photon')
        allocate(DiscretePhoton :: this%energy)
      case ('level')
        allocate(LevelInelastic :: this%energy)
      case ('continuous')
        allocate(ContinuousTabular :: this%energy)
      case ('maxwell')
        allocate(MaxwellEnergy :: this%energy)
      case ('evaporation')
        allocate(Evaporation :: this%energy)
      case ('watt')
        allocate(WattEnergy :: this%energy)
      case default
        call warning("Energy distribution type '" // trim(type) &
             // "' not implemented.")
      end select

      if (allocated(this % energy)) then
        call this%energy%from_hdf5(energy_group)
      end if

      call close_group(energy_group)
    end if
  end subroutine uncorrelated_from_hdf5

end module secondary_uncorrelated
