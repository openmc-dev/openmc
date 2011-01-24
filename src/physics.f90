module physics

  use global
  use geometry,    only: find_cell, dist_to_boundary, cross_boundary
  use types,       only: Neutron
  use mcnp_random, only: rang
  use output,      only: error

  implicit none

contains

!=====================================================================
! TRANSPORT encompasses the main logic for moving a particle through
! geometry.
!=====================================================================

  subroutine transport( neut )

    type(Neutron), pointer, intent(inout) :: neut

    character(250) :: msg
    integer :: i, surf
    logical :: alive
 
    real(8) :: d_to_boundary
    real(8) :: d_to_collision
    real(8) :: distance

    if ( neut%cell == 0 ) then
       call find_cell( neut )
    end if

    do while ( neut%alive )

       ! Determine distance neutron moves
       call dist_to_boundary( neut, d_to_boundary, surf )
       d_to_collision = -log(rang()) / 1.0
       distance = min( d_to_boundary, d_to_collision )

       ! Advance neutron
       neut%xyz = neut%xyz + distance*neut%uvw

       ! Add pathlength tallies

       if ( d_to_collision > d_to_boundary ) then
          neut%surface = surf
          neut%cell = 0
          call cross_boundary( neut )
       else
          ! collision
          msg = "Collision not implemented yet!"
          call error( msg )
       end if
       
    end do

  end subroutine transport

end module physics
