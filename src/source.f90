module source

  use global
  use mcnp_random, only: rang, RN_init_particle

  implicit none

contains

!=====================================================================
! INIT_SOURCE initializes particles in the source bank
!=====================================================================

  subroutine init_source()

    integer :: i,j
    real(8) :: r(3)

    ! Allocate fission and source banks
    allocate( source_bank(n_particles)    )
    allocate( fission_bank(3*n_particles) )

    ! Initialize first cycle source bank
    do i = 1, n_particles
       call RN_init_particle(int(i,8))
       r = (/ (rang(), j = 1,3) /)
       source_bank(i)%uid = i
       source_bank(i)%xyz = r
    end do

  end subroutine init_source

end module source
