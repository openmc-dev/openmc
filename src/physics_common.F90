module physics_common

  use constants
  use particle_header,        only: Particle
  use random_lcg,             only: prn
  use settings,               only: weight_cutoff, weight_survive

  implicit none

contains

!===============================================================================
! RUSSIAN_ROULETTE
!===============================================================================

  subroutine russian_roulette(p)

    type(Particle), intent(inout) :: p

    if (p % wgt < weight_cutoff) then
      if (prn() < p % wgt / weight_survive) then
        p % wgt = weight_survive
        p % last_wgt = p % wgt
      else
        p % wgt = ZERO
        p % last_wgt = ZERO
        p % alive = .false.
      end if
    end if

  end subroutine russian_roulette

end module physics_common
