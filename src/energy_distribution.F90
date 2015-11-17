module energy_distribution

  use constants, only: ZERO, ONE, TWO, PI, HISTOGRAM, LINEAR_LINEAR
  use endf_header, only: Tab1
  use interpolation, only: interpolate_tab1
  use math, only: maxwell_spectrum, watt_spectrum
  use random_lcg, only: prn
  use search, only: binary_search

!===============================================================================
! ENERGYDISTRIBUTION (abstract) defines an energy distribution that is a
! function of the incident energy of a projectile
!===============================================================================

  type, abstract :: EnergyDistribution
  contains
    procedure(iSampleEnergy), deferred :: sample
  end type EnergyDistribution

  abstract interface
    function iSampleEnergy(this, E_in) result(E_out)
      import EnergyDistribution
      class(EnergyDistribution), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8) :: E_out
    end function iSampleEnergy
  end interface

  type :: EnergyDistributionContainer
    class(EnergyDistribution), allocatable :: obj
  end type EnergyDistributionContainer

!===============================================================================
! Derived classes
!===============================================================================

  type, extends(EnergyDistribution) :: TabularEquiprobable
    integer :: n_region
    integer, allocatable :: breakpoints(:)
    integer, allocatable :: interpolation(:)
    real(8), allocatable :: energy_in(:)
    real(8), allocatable :: energy_out(:,:)
  contains
    procedure :: sample => equiprobable_sample
  end type TabularEquiprobable

  type, extends(EnergyDistribution) :: LevelInelastic
    real(8) :: threshold
    real(8) :: mass_ratio
  contains
    procedure :: sample => level_inelastic_sample
  end type LevelInelastic

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
    real(8), allocatable :: energy_in(:)
    type(CTTable), allocatable :: energy_out(:)
  contains
    procedure :: sample => continuous_sample
  end type ContinuousTabular

  type, extends(EnergyDistribution) :: MaxwellEnergy
    type(Tab1) :: theta
    real(8) :: u
  contains
    procedure :: sample => maxwellenergy_sample
  end type MaxwellEnergy

  type, extends(EnergyDistribution) :: Evaporation
    type(Tab1) :: theta
    real(8) :: u
  contains
    procedure :: sample => evaporation_sample
  end type Evaporation

  type, extends(EnergyDistribution) :: WattEnergy
    type(Tab1) :: a
    type(Tab1) :: b
    real(8) :: u
  contains
    procedure :: sample => watt_sample
  end type WattEnergy

  type, extends(EnergyDistribution) :: NBodyPhaseSpace
    integer :: n_bodies
    real(8) :: mass_ratio
    real(8) :: A
    real(8) :: Q
  contains
    procedure :: sample => nbody_sample
  end type NBodyPhaseSpace

contains

  function equiprobable_sample(this, E_in) result(E_out)
    class(TabularEquiprobable), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    integer :: i, k, l
    integer :: n_energy_in
    integer :: n_energy_out
    real(8) :: r              ! interpolation factor on incoming energy
    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1
    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l

    ! Determine number of incoming/outgoing energies
    n_energy_in = size(this%energy_in)
    n_energy_out = size(this%energy_out, 1)

    ! determine index on incoming energy grid and interpolation factor
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
    real(8), intent(in) :: E_in
    real(8) :: E_out

    integer :: i, k, l
    integer :: n_energy_in
    integer :: n_energy_out
    real(8) :: r              ! interpolation factor on incoming energy
    real(8) :: r1             ! random number on [0,1)
    real(8) :: frac           ! interpolation factor on outgoing energy
    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1
    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l
    real(8) :: p_l_k, p_l_k1  ! adjacent p on outgoing grid l
    real(8) :: c_k, c_k1      ! cumulative probability
    logical :: histogram_interp

    ! read number of interpolation regions and incoming energies
    if (this%n_region == 1) then
      histogram_interp = (this%interpolation(1) == 1)
    else
      histogram_interp = .false.
    end if

    ! find energy bin and calculate interpolation factor -- if the energy is
    ! outside the range of the tabulated energies, choose the first or last bins
    n_energy_in = size(this%energy_in)
    if (E_in < this%energy_in(1)) then
      i = 1
      r = ZERO
    elseif (E_in > this%energy_in(n_energy_in)) then
      i = n_energy_in - 1
      r = ONE
    else
      i = binary_search(this%energy_in, n_energy_in, E_in)
      r = (E_in - this%energy_in(i)) / &
           (this%energy_in(i+1) - this%energy_in(i))
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

    ! interpolation for energy E1 and EK
    n_energy_out = size(this%energy_out(i)%e_out)
    E_i_1 = this%energy_out(i)%e_out(1)
    E_i_K = this%energy_out(i)%e_out(n_energy_out)

    n_energy_out = size(this%energy_out(i+1)%e_out)
    E_i1_1 = this%energy_out(i+1)%e_out(1)
    E_i1_K = this%energy_out(i+1)%e_out(n_energy_out)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

    ! determine outgoing energy bin
    n_energy_out = size(this%energy_out(l)%e_out)
    r1 = prn()
    c_k = this%energy_out(l)%c(1)
    do k = 1, n_energy_out - 1
      c_k1 = this%energy_out(l)%c(k+1)
      if (r1 < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
    k = min(k, n_energy_out - 1)

    E_l_k = this%energy_out(l)%e_out(k)
    p_l_k = this%energy_out(l)%p(k)
    if (this%energy_out(l)%interpolation == HISTOGRAM) then
      ! Histogram interpolation
      if (p_l_k > ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k
      end if

    elseif (this%energy_out(l)%interpolation == LINEAR_LINEAR) then
      ! Linear-linear interpolation
      E_l_k1 = this%energy_out(l)%e_out(k+1)
      p_l_k1 = this%energy_out(l)%p(k+1)

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
    real(8), intent(in) :: E_in
    real(8) :: E_out

    real(8) :: theta

    ! Get temperature corresponding to incoming energy
    theta = interpolate_tab1(this%theta, E_in)

    do
      ! sample maxwell fission spectrum
      E_out = maxwell_spectrum(theta)

      ! accept energy based on restriction energy
      if (E_out <= E_in - this%u) exit
    end do
  end function maxwellenergy_sample

  function evaporation_sample(this, E_in) result(E_out)
    class(Evaporation), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    real(8) :: theta
    real(8) :: x, y, v

    ! Get temperature corresponding to incoming energy
    theta = interpolate_tab1(this%theta, E_in)

    y = (E_in - this%U)/theta
    v = 1 - exp(-y)

    ! sample outgoing energy based on evaporation spectrum probability
    ! density function
    do
      x = -log((ONE - v*prn())*(ONE - v*prn()))
      if (x <= y) exit
    end do

    E_out = x*theta
  end function evaporation_sample

  function watt_sample(this, E_in) result(E_out)
    class(WattEnergy), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    real(8) :: a, b

    ! determine Watt parameter 'a' from tabulated function
    a = interpolate_tab1(this%a, E_in)

    ! determine Watt parameter 'b' from tabulated function
    b = interpolate_tab1(this%b, E_in)

    do
      ! Sample energy-dependent Watt fission spectrum
      E_out = watt_spectrum(a, b)

      ! accept energy based on restriction energy
      if (E_out <= E_in - this%u) exit
    end do
  end function watt_sample

  function nbody_sample(this, E_in) result(E_out)
    class(NBodyPhaseSpace), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8) :: E_out

    real(8) :: Ap
    real(8) :: E_max
    real(8) :: x, y, v
    real(8) :: r1, r2, r3, r4, r5, r6

    ! determine E_max parameter
    Ap = this%mass_ratio
    E_max = (Ap - ONE)/Ap * (this%A/(this%A + ONE)*E_in + this%Q)

    ! x is essentially a Maxwellian distribution
    x = maxwell_spectrum(ONE)

    select case (this%n_bodies)
    case (3)
      y = maxwell_spectrum(ONE)
    case (4)
      r1 = prn()
      r2 = prn()
      r3 = prn()
      y = -log(r1*r2*r3)
    case (5)
      r1 = prn()
      r2 = prn()
      r3 = prn()
      r4 = prn()
      r5 = prn()
      r6 = prn()
      y = -log(r1*r2*r3*r4) - log(r5) * cos(PI/TWO*r6)**2
    end select

    ! now determine v and E_out
    v = x/(x+y)
    E_out = E_max * v
  end function nbody_sample

end module energy_distribution
