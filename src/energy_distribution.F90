module energy_distribution

  use constants,     only: ZERO, ONE, TWO, PI, HISTOGRAM, LINEAR_LINEAR
  use endf_header,   only: Tabulated1D
  use math,          only: maxwell_spectrum, watt_spectrum
  use random_lcg,    only: prn
  use search,        only: binary_search

!===============================================================================
! ENERGYDISTRIBUTION (abstract) defines an energy distribution that is a
! function of the incident energy of a projectile. Each derived type must
! implement a sample() function that returns a sampled outgoing energy given an
! incoming energy
!===============================================================================

  type, abstract :: EnergyDistribution
  contains
    procedure(energy_distribution_sample_), deferred :: sample
  end type EnergyDistribution

  abstract interface
    function energy_distribution_sample_(this, E_in) result(E_out)
      import EnergyDistribution
      class(EnergyDistribution), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8) :: E_out
    end function energy_distribution_sample_
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
  end type TabularEquiprobable

!===============================================================================
! LEVELINELASTIC gives the energy distribution for level inelastic scattering by
! neutrons as in ENDF MT=51--90.
!===============================================================================

  type, extends(EnergyDistribution) :: LevelInelastic
    real(8) :: threshold
    real(8) :: mass_ratio
  contains
    procedure :: sample => level_inelastic_sample
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


  function level_inelastic_sample(this, E_in) result(E_out)
    class(LevelInelastic), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    E_out = this%mass_ratio*(E_in - this%threshold)
  end function level_inelastic_sample


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

end module energy_distribution
