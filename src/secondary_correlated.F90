module secondary_correlated

  use angleenergy_header, only: AngleEnergy
  use constants, only: ZERO, ONE, TWO, HISTOGRAM, LINEAR_LINEAR
  use distribution_univariate, only: DistributionContainer
  use random_lcg, only: prn
  use search, only: binary_search

!===============================================================================
! CORRELATEDANGLEENERGY represents a correlated angle-energy distribution. This
! corresponds to ACE law 61 and ENDF File 6, LAW=1, LANG/=2.
!===============================================================================

  type AngleEnergyTable
    integer :: interpolation
    integer :: n_discrete
    real(8), allocatable :: e_out(:)
    real(8), allocatable :: p(:)
    real(8), allocatable :: c(:)
    type(DistributionContainer), allocatable :: angle(:)
  end type AngleEnergyTable

  type, extends(AngleEnergy) :: CorrelatedAngleEnergy
    integer :: n_region                      ! number of interpolation regions
    integer, allocatable :: breakpoints(:)   ! breakpoints of interpolation regions
    integer, allocatable :: interpolation(:) ! interpolation region codes
    real(8), allocatable :: energy(:)        ! incoming energies
    type(AngleEnergyTable), allocatable :: distribution(:) ! outgoing E/mu distributions
  contains
    procedure :: sample => correlated_sample
  end type CorrelatedAngleEnergy

contains

  subroutine correlated_sample(this, E_in, E_out, mu)
    class(CorrelatedAngleEnergy), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sapmled scattering cosine

    integer :: i, k, l        ! indices
    integer :: n_energy_in    ! number of incoming energies
    integer :: n_energy_out   ! number of outgoing energies
    real(8) :: r              ! interpolation factor on incoming energy
    real(8) :: r1             ! random number on [0,1)
    real(8) :: frac           ! interpolation factor on outgoing energy
    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1
    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l
    real(8) :: p_l_k, p_l_k1  ! adjacent p on outgoing grid l
    real(8) :: c_k, c_k1      ! cumulative probability

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Before the secondary distribution refactor, an isotropic polar cosine was
    ! always sampled but then overwritten with the polar cosine sampled from the
    ! correlated distribution. To preserve the random number stream, we keep
    ! this dummy sampling here but can remove it later (will change answers)
    mu = TWO*prn() - ONE
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! find energy bin and calculate interpolation factor -- if the energy is
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
    if (r > prn()) then
      l = i + 1
    else
      l = i
    end if

    ! interpolation for energy E1 and EK
    n_energy_out = size(this%distribution(i)%e_out)
    E_i_1 = this%distribution(i)%e_out(1)
    E_i_K = this%distribution(i)%e_out(n_energy_out)

    n_energy_out = size(this%distribution(i+1)%e_out)
    E_i1_1 = this%distribution(i+1)%e_out(1)
    E_i1_K = this%distribution(i+1)%e_out(n_energy_out)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

    ! determine outgoing energy bin
    n_energy_out = size(this%distribution(l)%e_out)
    r1 = prn()
    c_k = this%distribution(l)%c(1)
    do k = 1, n_energy_out - 1
      c_k1 = this%distribution(l)%c(k+1)
      if (r1 < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
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
    if (l == i) then
      E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
    else
      E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
    end if

    ! Find correlated angular distribution for closest outgoing energy bin
    if (r1 - c_k < c_k1 - r1) then
      mu = this%distribution(l)%angle(k)%obj%sample()
    else
      mu = this%distribution(l)%angle(k + 1)%obj%sample()
    end if
  end subroutine correlated_sample

end module secondary_correlated
