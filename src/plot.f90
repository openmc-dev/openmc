module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, dist_to_boundary, cross_surface, &
                             cross_lattice
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use particle_header, only: Particle

  implicit none

contains

!===============================================================================
! RUN_PLOT
!===============================================================================

  subroutine run_plot()

  end subroutine run_plot

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport_no_collision(p)

    type(Particle), pointer :: p

    integer        :: surf           ! surface which particle is on
    integer        :: last_cell      ! most recent cell particle was in
    real(8)        :: distance       ! distance particle travels
    logical        :: found_cell     ! found cell which particle is in?
    logical        :: in_lattice     ! is surface crossing in lattice?
    character(MAX_LINE_LEN) :: msg   ! output/error message
    type(Universe), pointer :: univ

    if (p % cell == 0) then
       univ => universes(BASE_UNIVERSE)
       call find_cell(univ, p, found_cell)

       ! if particle couldn't be located, print error
       if (.not. found_cell) then
          write(msg, '(A,3ES11.3)') & 
               "Could not locate cell for particle at: ", p % xyz
          call fatal_error(msg)
       end if
    end if

    ! find energy index, interpolation factor
    do while (p % alive)

       ! Find the distance to the nearest boundary
       call dist_to_boundary(p, distance, surf, in_lattice)

       ! Advance particle
       p%xyz = p%xyz + distance * p%uvw
       p%xyz_local = p%xyz_local + distance * p%uvw

       last_cell = p % cell
       p % cell = 0
       if (in_lattice) then
          p % surface = 0
          call cross_lattice(p)
       else
          p % surface = surf
          call cross_surface(p, last_cell)
       end if
       
    end do

  end subroutine transport_no_collision

end module plot
