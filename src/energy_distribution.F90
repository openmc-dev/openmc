module energy_distribution

  use hdf5

  use algorithm,     only: binary_search
  use constants,     only: ZERO, ONE, HALF, TWO, PI, HISTOGRAM, LINEAR_LINEAR
  use endf_header,   only: Tabulated1D
  use hdf5_interface
  use math,          only: maxwell_spectrum, watt_spectrum
  use random_lcg,    only: prn

!===============================================================================
! ENERGYDISTRIBUTION (abstract) defines an energy distribution that is a
! function of the incident energy of a projectile. Each derived type must
! implement a sample() function that returns a sampled outgoing energy given an
! incoming energy
!===============================================================================

  type, abstract :: EnergyDistribution
  contains
    procedure(energy_distribution_sample_), deferred :: sample
    procedure(energy_distribution_from_hdf5_), deferred :: from_hdf5
  end type EnergyDistribution

  abstract interface
    function energy_distribution_sample_(this, E_in) result(E_out)
      import EnergyDistribution
      class(EnergyDistribution), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8) :: E_out
    end function energy_distribution_sample_

    subroutine energy_distribution_from_hdf5_(this, group_id)
      import EnergyDistribution
      import HID_T
      class(EnergyDistribution), intent(inout) :: this
      integer(HID_T), intent(in) :: group_id
    end subroutine energy_distribution_from_hdf5_
  end interface

  type :: EnergyDistributionContainer
    class(EnergyDistribution), allocatable :: obj
  end type EnergyDistributionContainer

!===============================================================================
!                             Derived classes
!===============================================================================

!===============================================================================
! TABULAREQUIPROBABLE represents an energy distribution with tabular
! equiprobable energy bins as given in ACE law 1. This is an older
! representation that has largely been replaced with ACE laws 4, 44, and 61.
!===============================================================================

  type, extends(EnergyDistribution) :: TabularEquiprobable
    integer :: n_region                      ! number of interpolation regions
    integer, allocatable :: breakpoints(:)   ! breakpoints of interpolation regions
    integer, allocatable :: interpolation(:) ! interpolation region codes
    real(8), allocatable :: energy_in(:)     ! incoming energies
    real(8), allocatable :: energy_out(:,:)  ! table of outgoing energies for
                                             ! each incoming energy
  contains
    procedure :: sample => equiprobable_sample
    procedure :: from_hdf5 => equiprobable_from_hdf5
  end type TabularEquiprobable

!===============================================================================
! DISCRETEPHOTON gives the energy distribution for a discrete photon (usually
! used for photon production from an incident-neutron reaction)
!===============================================================================

  type, extends(EnergyDistribution) :: DiscretePhoton
    integer :: primary_flag
    real(8) :: energy
    real(8) :: A
  contains
    procedure :: sample => discrete_photon_sample
    procedure :: from_hdf5 => discrete_photon_from_hdf5
  end type DiscretePhoton

!===============================================================================
! LEVELINELASTIC gives the energy distribution for level inelastic scattering by
! neutrons as in ENDF MT=51--90.
!===============================================================================

  type, extends(EnergyDistribution) :: LevelInelastic
    real(8) :: threshold
    real(8) :: mass_ratio
  contains
    procedure :: sample => level_inelastic_sample
    procedure :: from_hdf5 => level_inelastic_from_hdf5
  end type LevelInelastic

!===============================================================================
! CONTINUOUSTABULAR gives an energy distribution represented as a tabular
! distribution with histogram or linear-linear interpolation. This corresponds
! to ACE law 4, which NJOY produces for a number of ENDF energy distributions.
!===============================================================================

  type CTTable
    integer :: interpolation
    integer :: n_discrete
    real(8), allocatable :: e_out(:)
    real(8), allocatable :: p(:)
    real(8), allocatable :: c(:)
  end type CTTable

  type, extends(EnergyDistribution) :: ContinuousTabular
    integer :: n_region
    integer, allocatable :: breakpoints(:)
    integer, allocatable :: interpolation(:)
    real(8), allocatable :: energy(:)
    type(CTTable), allocatable :: distribution(:)
  contains
    procedure :: sample => continuous_sample
    procedure :: from_hdf5 => continuous_from_hdf5
  end type ContinuousTabular

!===============================================================================
! MAXWELLENERGY gives the energy distribution of neutrons emitted from a Maxwell
! fission spectrum. This corresponds to ACE law 7 and ENDF File 5, LF=7.
!===============================================================================

  type, extends(EnergyDistribution) :: MaxwellEnergy
    type(Tabulated1D) :: theta ! incoming-energy-dependent parameter
    real(8) :: u        ! restriction energy
  contains
    procedure :: sample => maxwellenergy_sample
    procedure :: from_hdf5 => maxwellenergy_from_hdf5
  end type MaxwellEnergy

!===============================================================================
! EVAPORATION represents an evaporation spectrum corresponding to ACE law 9 and
! ENDF File 5, LF=9.
!===============================================================================

  type, extends(EnergyDistribution) :: Evaporation
    type(Tabulated1D) :: theta
    real(8) :: u
  contains
    procedure :: sample => evaporation_sample
    procedure :: from_hdf5 => evaporation_from_hdf5
  end type Evaporation

!===============================================================================
! WATTENERGY gives the energy distribution of neutrons emitted from a Watt
! fission spectrum. This corresponds to ACE law 11 and ENDF File 5, LF=11.
!===============================================================================

  type, extends(EnergyDistribution) :: WattEnergy
    type(Tabulated1D) :: a
    type(Tabulated1D) :: b
    real(8) :: u
  contains
    procedure :: sample => watt_sample
    procedure :: from_hdf5 => watt_from_hdf5
  end type WattEnergy

contains

  function equiprobable_sample(this, E_in) result(E_out)
    class(TabularEquiprobable), intent(in) :: this
    real(8), intent(in) :: E_in  ! incoming energy
    real(8)             :: E_out ! sampled outgoing energy

    integer :: i, k, l        ! indices
    integer :: n_energy_in    ! number of incoming energies
    integer :: n_energy_out   ! number of outgoing energies
    real(8) :: r              ! interpolation factor on incoming energy
    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1
    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l

    ! Determine number of incoming/outgoing energies
    n_energy_in = size(this%energy_in)
    n_energy_out = size(this%energy_out, 1)

    ! Determine index on incoming energy grid and interpolation factor
    i = binary_search(this%energy_in, size(this%energy_in), E_in)
    r = (E_in - this%energy_in(i)) / &
         (this%energy_in(i+1) - this%energy_in(i))

    ! Sample outgoing energy bin
    k = 1 + int(n_energy_out * prn())

    ! Determine E_1 and E_K
    E_i_1 = this%energy_out(1, i)
    E_i_K = this%energy_out(n_energy_out, i)

    E_i1_1 = this%energy_out(1, i+1)
    E_i1_K = this%energy_out(n_energy_out, i+1)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

    ! Randomly select between the outgoing table for incoming energy E_i and
    ! E_(i+1)
    if (prn() < r) then
      l = i + 1
    else
      l = i
    end if

    ! Determine E_l_k and E_l_k+1
    E_l_k  = this%energy_out(k, l)
    E_l_k1 = this%energy_out(k+1, l)

    ! Determine E' (denoted here as E_out)
    E_out  = E_l_k + prn()*(E_l_k1 - E_l_k)

    ! Now interpolate between incident energy bins i and i + 1
    if (l == i) then
      E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
    else
      E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
    end if
  end function equiprobable_sample

  subroutine equiprobable_from_hdf5(this, group_id)
    class(TabularEquiprobable), intent(inout) :: this
    integer(HID_T),             intent(in)    :: group_id
  end subroutine equiprobable_from_hdf5

  function discrete_photon_sample(this, E_in) result(E_out)
    class(DiscretePhoton), intent(in) :: this
    real(8),               intent(in) :: E_in
    real(8)                           :: E_out

    if (this % primary_flag == 2) then
      E_out = this % energy + this % A/(this % A + 1)*E_in
    else
      E_out = this % energy
    end if
  end function discrete_photon_sample

  subroutine discrete_photon_from_hdf5(this, group_id)
    class(DiscretePhoton), intent(inout) :: this
    integer(HID_T),        intent(in)    :: group_id

    call read_attribute(this % primary_flag, group_id, 'primary_flag')
    call read_attribute(this % energy, group_id, 'energy')
    call read_attribute(this % A, group_id, 'atomic_weight_ratio')
  end subroutine discrete_photon_from_hdf5

  function level_inelastic_sample(this, E_in) result(E_out)
    class(LevelInelastic), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    E_out = this%mass_ratio*(E_in - this%threshold)
  end function level_inelastic_sample

  subroutine level_inelastic_from_hdf5(this, group_id)
    class(LevelInelastic), intent(inout) :: this
    integer(HID_T),        intent(in)    :: group_id

    call read_attribute(this%threshold, group_id, 'threshold')
    call read_attribute(this%mass_ratio, group_id, 'mass_ratio')
  end subroutine level_inelastic_from_hdf5

  function continuous_sample(this, E_in) result(E_out)
    class(ContinuousTabular), intent(in) :: this
    real(8), intent(in) :: E_in  ! incoming energy
    real(8)             :: E_out ! sampled outgoing energy

    integer :: i, k, l          ! indices
    integer :: n_energy_in      ! number of incoming energies
    integer :: n_energy_out     ! number of outgoing energies
    real(8) :: r                ! interpolation factor on incoming energy
    real(8) :: r1               ! random number on [0,1)
    real(8) :: frac             ! interpolation factor on outgoing energy
    real(8) :: E_i_1, E_i_K     ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K   ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K         ! endpoints interpolated between i and i+1
    real(8) :: E_l_k, E_l_k1    ! adjacent E on outgoing grid l
    real(8) :: p_l_k, p_l_k1    ! adjacent p on outgoing grid l
    real(8) :: c_k, c_k1        ! cumulative probability
    logical :: histogram_interp ! whether histogram interpolation is used

    ! Read number of interpolation regions and incoming energies
    if (this%n_region == 1) then
      histogram_interp = (this%interpolation(1) == 1)
    else
      histogram_interp = .false.
    end if

    ! Find energy bin and calculate interpolation factor -- if the energy is
    ! outside the range of the tabulated energies, choose the first or last bins
    n_energy_in = size(this%energy)
    if (E_in < this%energy(1)) then
      i = 1
      r = ZERO
    elseif (E_in > this%energy(n_energy_in)) then
      i = n_energy_in - 1
      r = ONE
    else
      i = binary_search(this%energy, n_energy_in, E_in)
      r = (E_in - this%energy(i)) / &
           (this%energy(i+1) - this%energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (histogram_interp) then
      l = i
    else
      if (r > prn()) then
        l = i + 1
      else
        l = i
      end if
    end if

    ! Interpolation for energy E1 and EK
    n_energy_out = size(this%distribution(i)%e_out)
    E_i_1 = this%distribution(i)%e_out(1)
    E_i_K = this%distribution(i)%e_out(n_energy_out)

    n_energy_out = size(this%distribution(i+1)%e_out)
    E_i1_1 = this%distribution(i+1)%e_out(1)
    E_i1_K = this%distribution(i+1)%e_out(n_energy_out)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

    ! Determine outgoing energy bin
    n_energy_out = size(this%distribution(l)%e_out)
    r1 = prn()
    c_k = this%distribution(l)%c(1)
    do k = 1, n_energy_out - 1
      c_k1 = this%distribution(l)%c(k+1)
      if (r1 < c_k1) exit
      c_k = c_k1
    end do

    ! Check to make sure k is <= NP - 1
    k = min(k, n_energy_out - 1)

    E_l_k = this%distribution(l)%e_out(k)
    p_l_k = this%distribution(l)%p(k)
    if (this%distribution(l)%interpolation == HISTOGRAM) then
      ! Histogram interpolation
      if (p_l_k > ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k
      end if

    elseif (this%distribution(l)%interpolation == LINEAR_LINEAR) then
      ! Linear-linear interpolation
      E_l_k1 = this%distribution(l)%e_out(k+1)
      p_l_k1 = this%distribution(l)%p(k+1)

      frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
      if (frac == ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k + (sqrt(max(ZERO, p_l_k*p_l_k + &
             TWO*frac*(r1 - c_k))) - p_l_k)/frac
      end if
    end if

    ! Now interpolate between incident energy bins i and i + 1
    if (.not. histogram_interp) then
      if (l == i) then
        E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
      else
        E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
      end if
    end if
  end function continuous_sample

  subroutine continuous_from_hdf5(this, group_id)
    class(ContinuousTabular), intent(inout) :: this
    integer(HID_T),           intent(in)    :: group_id

    integer :: i, j, k
    integer :: n
    integer :: n_energy
    integer(HID_T) :: dset_id
    integer(HSIZE_T) :: dims(1), dims2(2)
    integer, allocatable :: temp(:,:)
    integer, allocatable :: offsets(:)
    integer, allocatable :: interp(:)
    integer, allocatable :: n_discrete(:)
    real(8), allocatable :: eout(:,:)

    ! Open incoming energy dataset
    dset_id = open_dataset(group_id, 'energy')

    ! Get interpolation parameters
    call read_attribute(temp, dset_id, 'interpolation')
    allocate(this%breakpoints(size(temp, 1)))
    allocate(this%interpolation(size(temp, 1)))
    this%breakpoints(:) = temp(:, 1)
    this%interpolation(:) = temp(:, 2)
    this%n_region = size(this%breakpoints)

    ! Get incoming energies
    call get_shape(dset_id, dims)
    n_energy = int(dims(1), 4)
    allocate(this%energy(n_energy))
    allocate(this%distribution(n_energy))
    call read_dataset(this%energy, dset_id)
    call close_dataset(dset_id)

    ! Get outgoing energy distribution data
    dset_id = open_dataset(group_id, 'distribution')
    call read_attribute(offsets, dset_id, 'offsets')
    call read_attribute(interp, dset_id, 'interpolation')
    call read_attribute(n_discrete, dset_id, 'n_discrete_lines')
    call get_shape(dset_id, dims2)
    allocate(eout(dims2(1), dims2(2)))
    call read_dataset(eout, dset_id)
    call close_dataset(dset_id)

    do i = 1, n_energy
      ! Determine number of outgoing energies
      j = offsets(i)
      if (i < n_energy) then
        n = offsets(i+1) - j
      else
        n = size(eout, 1) - j
      end if

      associate (d => this % distribution(i))
        ! Assign interpolation scheme and number of discrete lines
        d % interpolation = interp(i)
        d % n_discrete = n_discrete(i)

        ! Allocate arrays for energies and PDF/CDF
        allocate(d % e_out(n))
        allocate(d % p(n))
        allocate(d % c(n))

        ! Copy data
        d % e_out(:) = eout(j+1:j+n, 1)
        d % p(:) = eout(j+1:j+n, 2)

        ! To get answers that match ACE data, for now we still use the tabulated
        ! CDF values that were passed through to the HDF5 library. At a later
        ! time, we can remove the CDF values from the HDF5 library and
        ! reconstruct them using the PDF
        if (.true.) then
          d % c(:) = eout(j+1:j+n, 3)
        else
          ! Calculate cumulative distribution function -- discrete portion
          do k = 1, n_discrete(i)
            if (k == 1) then
              d % c(k) = d % p(k)
            else
              d % c(k) = d % c(k-1) + d % p(k)
            end if
          end do

          ! Continuous portion
          do k = d % n_discrete + 1, n
            if (k == d % n_discrete + 1) then
              d % c(k) = sum(d % p(1:d % n_discrete))
            else
              if (d % interpolation == HISTOGRAM) then
                d % c(k) = d % c(k-1) + d % p(k-1) * &
                     (d % e_out(k) - d % e_out(k-1))
              elseif (d % interpolation == LINEAR_LINEAR) then
                d % c(k) = d % c(k-1) + HALF*(d % p(k-1) + d % p(k)) * &
                     (d % e_out(k) - d % e_out(k-1))
              end if
            end if
          end do

          ! Normalize density and distribution functions
          d % p(:) = d % p(:)/d % c(n)
          d % c(:) = d % c(:)/d % c(n)
        end if
      end associate
    end do
  end subroutine continuous_from_hdf5

  function maxwellenergy_sample(this, E_in) result(E_out)
    class(MaxwellEnergy), intent(in) :: this
    real(8), intent(in) :: E_in  ! incoming energy
    real(8)             :: E_out ! sampled outgoing energy

    real(8) :: theta ! Maxwell distribution parameter

    ! Get temperature corresponding to incoming energy
    theta = this % theta % evaluate(E_in)

    do
      ! Sample maxwell fission spectrum
      E_out = maxwell_spectrum(theta)

      ! Accept energy based on restriction energy
      if (E_out <= E_in - this%u) exit
    end do
  end function maxwellenergy_sample

  subroutine maxwellenergy_from_hdf5(this, group_id)
    class(MaxwellEnergy), intent(inout) :: this
    integer(HID_T),       intent(in)    :: group_id

    integer(HID_T) :: dset_id

    call read_attribute(this%u, group_id, 'u')
    dset_id = open_dataset(group_id, 'theta')
    call this%theta%from_hdf5(dset_id)
    call close_dataset(dset_id)
  end subroutine maxwellenergy_from_hdf5

  function evaporation_sample(this, E_in) result(E_out)
    class(Evaporation), intent(in) :: this
    real(8), intent(in) :: E_in  ! incoming energy
    real(8)             :: E_out ! sampled outgoing energy

    real(8) :: theta    ! evaporation spectrum parameter
    real(8) :: x, y, v

    ! Get temperature corresponding to incoming energy
    theta = this % theta % evaluate(E_in)

    y = (E_in - this%u)/theta
    v = 1 - exp(-y)

    ! Sample outgoing energy based on evaporation spectrum probability
    ! density function
    do
      x = -log((ONE - v*prn())*(ONE - v*prn()))
      if (x <= y) exit
    end do

    E_out = x*theta
  end function evaporation_sample

  subroutine evaporation_from_hdf5(this, group_id)
    class(Evaporation), intent(inout) :: this
    integer(HID_T),     intent(in)    :: group_id

    integer(HID_T) :: dset_id

    call read_attribute(this%u, group_id, 'u')
    dset_id = open_dataset(group_id, 'theta')
    call this%theta%from_hdf5(dset_id)
    call close_dataset(dset_id)
  end subroutine evaporation_from_hdf5

  function watt_sample(this, E_in) result(E_out)
    class(WattEnergy), intent(in) :: this
    real(8), intent(in) :: E_in  ! incoming energy
    real(8)             :: E_out ! sampled outgoing energy

    real(8) :: a, b ! Watt spectrum parameters

    ! Determine Watt parameter 'a' from tabulated function
    a = this % a % evaluate(E_in)

    ! Determine Watt parameter 'b' from tabulated function
    b = this % b % evaluate(E_in)

    do
      ! Sample energy-dependent Watt fission spectrum
      E_out = watt_spectrum(a, b)

      ! Accept energy based on restriction energy
      if (E_out <= E_in - this%u) exit
    end do
  end function watt_sample

  subroutine watt_from_hdf5(this, group_id)
    class(WattEnergy), intent(inout) :: this
    integer(HID_T),    intent(in)    :: group_id

    integer(HID_T) :: dset_id

    call read_attribute(this%u, group_id, 'u')

    dset_id = open_dataset(group_id, 'a')
    call this%a%from_hdf5(dset_id)
    call close_dataset(dset_id)

    dset_id = open_dataset(group_id, 'b')
    call this%b%from_hdf5(dset_id)
    call close_dataset(dset_id)
  end subroutine watt_from_hdf5

end module energy_distribution
