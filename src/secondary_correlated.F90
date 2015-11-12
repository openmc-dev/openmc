module secondary_correlated

  use constants, only: ZERO, ONE, TWO, HISTOGRAM, LINEAR_LINEAR
  use distribution_univariate, only: Tabular
  use error, only: fatal_error
  use secondary_header, only: SecondaryDistribution
  use random_lcg, only: prn
  use search, only: binary_search

  type AngleEnergyTable
    type(Tabular) :: energy
    type(Tabular), allocatable :: angle(:)
  end type AngleEnergyTable

  type, extends(SecondaryDistribution) :: CorrelatedAngleEnergy
    integer :: n_region
    integer, allocatable :: breakpoints(:)
    integer, allocatable :: interpolation(:)
    real(8), allocatable :: energy_in(:)
    type(AngleEnergyTable), allocatable :: table(:)
  contains
    procedure :: sample => correlated_sample
  end type CorrelatedAngleEnergy

contains

  subroutine correlated_sample(this, E_in, E_out, mu)
    class(CorrelatedAngleEnergy), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8), intent(out) :: E_out
    real(8), intent(out) :: mu

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

    ! TODO: Write error during initialization
    if (this%n_region > 1) then
      call fatal_error("Multiple interpolation regions not supported while &
           &attempting to sample Kalbach-Mann distribution.")
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
    if (r > prn()) then
      l = i + 1
    else
      l = i
    end if

    ! interpolation for energy E1 and EK
    n_energy_out = size(this%table(i)%energy%x)
    E_i_1 = this%table(i)%energy%x(1)
    E_i_K = this%table(i)%energy%x(n_energy_out)

    n_energy_out = size(this%table(i+1)%energy%x)
    E_i1_1 = this%table(i+1)%energy%x(1)
    E_i1_K = this%table(i+1)%energy%x(n_energy_out)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

!!$    ! TODO: Write error at initizliation
!!$    if (this%table(l)%n_discrete > 0) then
!!$      ! discrete lines present
!!$      call fatal_error("Discrete lines in continuous tabular distributed not &
!!$           &yet supported")
!!$    end if

    ! determine outgoing energy bin
    n_energy_out = size(this%table(l)%energy%x)
    r1 = prn()
    c_k = this%table(l)%energy%c(1)
    do k = 1, n_energy_out - 1
      c_k1 = this%table(l)%energy%c(k+1)
      if (r1 < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
    k = min(k, n_energy_out - 1)

    E_l_k = this%table(l)%energy%x(k)
    p_l_k = this%table(l)%energy%p(k)
    if (this%table(l)%energy%interpolation == HISTOGRAM) then
      ! Histogram interpolation
      if (p_l_k > ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k
      end if

    elseif (this%table(l)%energy%interpolation == LINEAR_LINEAR) then
      ! Linear-linear interpolation
      E_l_k1 = this%table(l)%energy%x(k+1)
      p_l_k1 = this%table(l)%energy%p(k+1)

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
      mu = this%table(l)%angle(k)%sample()
    else
      mu = this%table(l)%angle(k + 1)%sample()
    end if
  end subroutine correlated_sample

end module secondary_correlated
