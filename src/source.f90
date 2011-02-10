module source

  use global
  use mcnp_random, only: rang, RN_init_particle

  implicit none

  integer :: source_index

contains

!=====================================================================
! INIT_SOURCE initializes particles in the source bank
!=====================================================================

  subroutine init_source()

    integer :: i,j
    real(8) :: r(3)
    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    real(8) :: p_min(3), p_max(3)

    ! Allocate fission and source banks
    allocate( source_bank(n_particles)    )
    allocate( fission_bank(3*n_particles) )

    ! Check external source type
    if ( external_source%type == SRC_BOX ) then
       p_min = external_source%values(1:3)
       p_max = external_source%values(4:6)
    end if

    ! Initialize first cycle source bank
    do i = 1, n_particles
       call RN_init_particle(int(i,8))

       ! position
       r = (/ (rang(), j = 1,3) /)
       source_bank(i)%uid = i
       source_bank(i)%xyz = p_min + r*(p_max - p_min)

       ! angle
       phi = 2.0_8*PI*rang()
       mu = 2.0_8*rang() - 1.0_8
       source_bank(i)%uvw(1) = mu
       source_bank(i)%uvw(2) = sqrt(1. - mu**2) * cos(phi)
       source_bank(i)%uvw(3) = sqrt(1. - mu**2) * sin(phi)

       ! cell
       source_bank(i)%cell = 0

       ! set particle to be alive
       source_bank(i)%alive = .true.
    end do

    ! Reset source index
    source_index = 0

  end subroutine init_source

!=====================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!=====================================================================

  function get_source_particle()

    type(Neutron), pointer :: get_source_particle

    ! increment index
    source_index = source_index + 1

    ! if at end of bank, return null pointer
    if ( source_index > size(source_bank) ) then
       get_source_particle => null()
       return
    end if

    ! point to next source particle
    get_source_particle => source_bank(source_index)

  end function get_source_particle

end module source
