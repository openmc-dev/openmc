module physics_common

  use, intrinsic :: ISO_C_BINDING

  use particle_header,        only: Particle

  implicit none

!===============================================================================
! RUSSIAN_ROULETTE FROM C
!===============================================================================

  interface
    subroutine russian_roulette(p) bind(C)
      import Particle
      type(Particle), intent(inout) :: p
    end subroutine russian_roulette
  end interface

end module physics_common
