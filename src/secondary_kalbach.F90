module secondary_kalbach

  use constants, only: ZERO, ONE, TWO, HISTOGRAM, LINEAR_LINEAR
  use error, only: fatal_error
  use secondary_header, only: SecondaryDistribution
  use random_lcg, only: prn
  use search, only: binary_search

  type KalbachMannTable
    integer :: n_discrete
    integer :: interpolation
    real(8), allocatable :: e_out(:)
    real(8), allocatable :: p(:)
    real(8), allocatable :: c(:)
    real(8), allocatable :: r(:)
    real(8), allocatable :: a(:)
  end type KalbachMannTable

  type, extends(SecondaryDistribution) :: KalbachMann
    integer :: n_region
    integer, allocatable :: breakpoints(:)
    integer, allocatable :: interpolation(:)
    real(8), allocatable :: energy_in(:)
    type(KalbachMannTable), allocatable :: table(:)
  contains
    procedure :: sample => kalbachmann_sample
  end type KalbachMann

contains

  subroutine kalbachmann_sample(this, E_in, E_out, mu)
    class(KalbachMann), intent(in) :: this
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
    real(8) :: km_r, km_a     ! Kalbach-Mann parameters
    real(8) :: T

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
    n_energy_out = size(this%table(i)%e_out)
    E_i_1 = this%table(i)%e_out(1)
    E_i_K = this%table(i)%e_out(n_energy_out)

    n_energy_out = size(this%table(i+1)%e_out)
    E_i1_1 = this%table(i+1)%e_out(1)
    E_i1_K = this%table(i+1)%e_out(n_energy_out)

    E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
    E_K = E_i_K + r*(E_i1_K - E_i_K)

    ! TODO: Write error at initizliation
    if (this%table(l)%n_discrete > 0) then
      ! discrete lines present
      call fatal_error("Discrete lines in continuous tabular distributed not &
           &yet supported")
    end if

    ! determine outgoing energy bin
    n_energy_out = size(this%table(l)%e_out)
    r1 = prn()
    c_k = this%table(l)%c(1)
    do k = 1, n_energy_out - 1
      c_k1 = this%table(l)%c(k+1)
      if (r1 < c_k1) exit
      c_k = c_k1
    end do

    ! check to make sure k is <= NP - 1
    k = min(k, n_energy_out - 1)

    E_l_k = this%table(l)%e_out(k)
    p_l_k = this%table(l)%p(k)
    if (this%table(l)%interpolation == HISTOGRAM) then
      ! Histogram interpolation
      if (p_l_k > ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k
      end if

      ! Determine Kalbach-Mann parameters
      km_r = this%table(l)%r(k)
      km_a = this%table(l)%a(k)

    elseif (this%table(l)%interpolation == LINEAR_LINEAR) then
      ! Linear-linear interpolation
      E_l_k1 = this%table(l)%e_out(k+1)
      p_l_k1 = this%table(l)%p(k+1)

      frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
      if (frac == ZERO) then
        E_out = E_l_k + (r1 - c_k)/p_l_k
      else
        E_out = E_l_k + (sqrt(max(ZERO, p_l_k*p_l_k + &
             TWO*frac*(r1 - c_k))) - p_l_k)/frac
      end if

      ! Determine Kalbach-Mann parameters
      km_r = this%table(l)%r(k) + (E_out - E_l_k)/(E_l_k1 - E_l_k) * &
           (this%table(l)%r(k+1) - this%table(l)%r(k))
      km_a = this%table(l)%a(k) + (E_out - E_l_k)/(E_l_k1 - E_l_k) * &
           (this%table(l)%a(k+1) - this%table(l)%a(k))
    end if

    ! Now interpolate between incident energy bins i and i + 1
    if (l == i) then
      E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
    else
      E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
    end if

    ! Sampled correlated angle from Kalbach-Mann parameters
    if (prn() > km_r) then
      T = (TWO*prn() - ONE) * sinh(km_a)
      mu = log(T + sqrt(T*T + ONE))/km_a
    else
      r1 = prn()
      mu = log(r1*exp(km_a) + (ONE - r1)*exp(-km_a))/km_a
    end if
  end subroutine kalbachmann_sample

end module secondary_kalbach
