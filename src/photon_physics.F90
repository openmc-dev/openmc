module photon_physics

  use algorithm,       only: binary_search
  use constants
  use particle_header
  use photon_header,   only: BremsstrahlungData, ttb_e_grid, ttb
  use random_lcg,      only: prn
  use settings

contains

!===============================================================================
! THICK_TARGET_BREMSSTRAHLUNG
!===============================================================================

  subroutine thick_target_bremsstrahlung(p, E_lost) bind(C)
    type(Particle), intent(inout) :: p
    real(C_DOUBLE), intent(out) :: E_lost

    integer :: i, j
    integer :: i_e, i_w
    integer :: n
    integer :: n_e
    real(8) :: a
    real(8) :: f
    real(8) :: e, e_l, e_r
    real(8) :: y, y_l, y_r
    real(8) :: w, w_l, w_r
    real(8) :: p_l, p_r
    real(8) :: c, c_l, c_max
    type(BremsstrahlungData), pointer :: mat

    if (p % material == MATERIAL_VOID) return

    if (p % E < energy_cutoff(PHOTON)) return

    ! Get bremsstrahlung data for this material and particle type
    if (p % type == POSITRON) then
      mat => ttb(p % material) % positron
    else
      mat => ttb(p % material) % electron
    end if

    e = log(p % E)
    n_e = size(ttb_e_grid)

    ! Find the lower bounding index of the incident electron energy
    j = binary_search(ttb_e_grid, n_e, e)
    if (j == n_e) j = j - 1

    ! Get the interpolation bounds
    e_l = ttb_e_grid(j)
    e_r = ttb_e_grid(j+1)
    y_l = mat % yield(j)
    y_r = mat % yield(j+1)

    ! Calculate the interpolation weight w_j+1 of the bremsstrahlung energy PDF
    ! interpolated in log energy, which can be interpreted as the probability
    ! of index j+1
    f = (e - e_l)/(e_r - e_l)

    ! Get the photon number yield for the given energy using linear
    ! interpolation on a log-log scale
    y = exp(y_l + (y_r - y_l)*f)

    ! Sample number of secondary bremsstrahlung photons
    n = int(y + prn())

    E_lost = ZERO
    if (n == 0) return

    ! Sample index of the tabulated PDF in the energy grid, j or j+1
    if (prn() <= f .or. j == 1) then
      i_e = j + 1

      ! Interpolate the maximum value of the CDF at the incoming particle
      ! energy on a log-log scale
      p_l = mat % pdf(i_e-1, i_e)
      p_r = mat % pdf(i_e, i_e)
      c_l = mat % cdf(i_e-1, i_e)
      a = log(p_r/p_l)/(e_r - e_l) + ONE
      c_max = c_l + exp(e_l)*p_l/a*(exp(a*(e - e_l)) - ONE)
    else
      i_e = j

      ! Maximum value of the CDF
      c_max = mat % cdf(i_e, i_e)
    end if

    ! Sample the energies of the emitted photons
    do i = 1, n
      ! Generate a random number r and determine the index i for which
      ! cdf(i) <= r*cdf,max <= cdf(i+1)
      c = prn()*c_max
      i_w = binary_search(mat % cdf(:i_e,i_e), i_e, c)

      ! Sample the photon energy
      w_l = ttb_e_grid(i_w)
      w_r = ttb_e_grid(i_w+1)
      p_l = mat % pdf(i_w, i_e)
      p_r = mat % pdf(i_w+1, i_e)
      c_l = mat % cdf(i_w, i_e)
      a = log(p_r/p_l)/(w_r - w_l) + ONE
      w = exp(w_l)*(a*(c - c_l)/(exp(w_l)*p_l) + ONE)**(ONE/a)

      if (w > energy_cutoff(PHOTON)) then
        ! Create secondary photon
        call particle_create_secondary(p, p % coord(1) % uvw, w, PHOTON, &
             run_ce=.true._C_BOOL)
        E_lost = E_lost + w
      end if
    end do

  end subroutine thick_target_bremsstrahlung

end module photon_physics
