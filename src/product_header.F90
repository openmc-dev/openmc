module product_header

  use hdf5, only: HID_T

  use angleenergy_header, only: AngleEnergyContainer
  use constants, only: ZERO, MAX_WORD_LEN, EMISSION_PROMPT, EMISSION_DELAYED, &
       EMISSION_TOTAL, NEUTRON, PHOTON
  use endf_header, only: Tabulated1D, Function1D, Polynomial
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       open_dataset, close_dataset, read_dataset
  use random_lcg, only: prn
  use secondary_correlated, only: CorrelatedAngleEnergy
  use secondary_kalbach, only: KalbachMann
  use secondary_nbody, only: NBodyPhaseSpace
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use string, only: to_str

!===============================================================================
! REACTIONPRODUCT stores a data for a reaction product including its yield and
! angle-energy distributions, each of which has a given probability of occurring
! for a given incoming energy. In general, most products only have one
! angle-energy distribution, but for some cases (e.g., (n,2n) in certain
! nuclides) multiple distinct distributions exist.
!===============================================================================

  type :: ReactionProduct
    integer :: particle
    integer :: emission_mode            ! prompt, delayed, or total emission
    real(8) :: decay_rate               ! Decay rate for delayed neutron precursors
    class(Function1D), pointer :: yield => null()  ! Energy-dependent neutron yield
    type(Tabulated1D), allocatable :: applicability(:)
    type(AngleEnergyContainer), allocatable :: distribution(:)
  contains
    procedure :: sample => reactionproduct_sample
    procedure :: from_hdf5 => reactionproduct_from_hdf5
  end type ReactionProduct

contains

  subroutine reactionproduct_sample(this, E_in, E_out, mu)
    class(ReactionProduct), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sampled scattering cosine

    integer :: i       ! loop counter
    integer :: n       ! number of angle-energy distributions
    real(8) :: prob    ! cumulative probability
    real(8) :: c       ! sampled cumulative probability

    n = size(this%applicability)
    if (n > 1) then
      prob = ZERO
      c = prn()
      do i = 1, n
        ! Determine probability that i-th energy distribution is sampled
        prob = prob + this % applicability(i) % evaluate(E_in)

        ! If i-th distribution is sampled, sample energy from the distribution
        if (c <= prob) then
          call this%distribution(i)%obj%sample(E_in, E_out, mu)
          exit
        end if
      end do
    else
      ! If only one distribution is present, go ahead and sample it
      call this%distribution(1)%obj%sample(E_in, E_out, mu)
    end if

  end subroutine reactionproduct_sample

  subroutine reactionproduct_from_hdf5(this, group_id)
    class(ReactionProduct), intent(inout) :: this
    integer(HID_T), intent(in)    :: group_id

    integer :: i
    integer :: n
    integer(HID_T) :: dgroup
    integer(HID_T) :: app
    integer(HID_T) :: yield
    character(MAX_WORD_LEN) :: temp

    ! Read particle type
    call read_attribute(temp, group_id, 'particle')
    select case (temp)
    case ('neutron')
      this % particle = NEUTRON
    case ('photon')
      this % particle = PHOTON
    end select

    ! Read emission mode and decay rate
    call read_attribute(temp, group_id, 'emission_mode')
    select case (temp)
    case ('prompt')
      this % emission_mode = EMISSION_PROMPT
    case ('delayed')
      this % emission_mode = EMISSION_DELAYED
    case ('total')
      this % emission_mode = EMISSION_TOTAL
    end select

    ! Read decay rate for delayed emission
    if (this % emission_mode == EMISSION_DELAYED) then
      call read_attribute(this % decay_rate, group_id, 'decay_rate')
    end if

    ! Read secondary particle yield
    yield = open_dataset(group_id, 'yield')
    call read_attribute(temp, yield, 'type')
    select case (temp)
    case ('Tabulated1D')
      allocate(Tabulated1D :: this % yield)
    case ('Polynomial')
      allocate(Polynomial :: this % yield)
    end select
    call this % yield % from_hdf5(yield)
    call close_dataset(yield)

    call read_attribute(n, group_id, 'n_distribution')
    allocate(this%applicability(n))
    allocate(this%distribution(n))

    do i = 1, n
      dgroup = open_group(group_id, trim('distribution_' // to_str(i - 1)))

      ! Read applicability
      if (n > 1) then
        app = open_dataset(dgroup, 'applicability')
        call this%applicability(i)%from_hdf5(app)
        call close_dataset(app)
      end if

      ! Read type of distribution and allocate accordingly
      call read_attribute(temp, dgroup, 'type')
      select case (temp)
      case ('uncorrelated')
        allocate(UncorrelatedAngleEnergy :: this%distribution(i)%obj)
      case ('correlated')
        allocate(CorrelatedAngleEnergy :: this%distribution(i)%obj)
      case ('nbody')
        allocate(NBodyPhaseSpace :: this%distribution(i)%obj)
      case ('kalbach-mann')
        allocate(KalbachMann :: this%distribution(i)%obj)
      end select

      ! Read distribution data
      call this%distribution(i)%obj%from_hdf5(dgroup)

      call close_group(dgroup)
    end do

  end subroutine reactionproduct_from_hdf5

end module product_header
