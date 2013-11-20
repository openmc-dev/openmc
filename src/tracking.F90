module tracking

  use cross_section,   only: calculate_xs
  use error,           only: fatal_error, warning
  use geometry,        only: find_cell, distance_to_boundary, &
                             cross_lattice, check_cell_overlap, &
                             handle_lost_particle
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use mesh
  use output,          only: write_message
  use particle_header, only: LocalCoord, Particle, deallocate_coord
  use physics,         only: collision
  use random_lcg,      only: prn
  use string,          only: to_str
  use tally,           only: score_analog_tally, score_tracklength_tally, &
                             score_surface_current
                             

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), intent(inout) :: p

    integer :: surface_crossed ! surface which particle is on
    integer :: lattice_crossed ! lattice boundary which particle crossed
    integer :: last_cell       ! most recent cell particle was in
    integer :: n_event         ! number of collisions/crossings
    real(8) :: d_boundary      ! distance to nearest boundary
    real(8) :: d_collision     ! sampled distance to collision
    real(8) :: distance        ! distance particle travels
    logical :: found_cell      ! found cell which particle is in?
    type(LocalCoord), pointer, save :: coord => null()
!$omp threadprivate(coord)

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      message = "Simulating Particle " // trim(to_str(p % id))
      call write_message()
    end if

    ! If the cell hasn't been determined based on the particle's location,
    ! initiate a search for the current cell
    if (p % coord % cell == NONE) then
      call find_cell(p, found_cell)

      ! Particle couldn't be located
      if (.not. found_cell) then
        message = "Could not locate particle " // trim(to_str(p % id))
        call fatal_error()
      end if

      ! set birth cell attribute
      p % cell_born = p % coord % cell
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
!$omp critical
    total_weight = total_weight + p % wgt
!$omp end critical

    ! Force calculation of cross-sections by setting last energy to zero 
    micro_xs % last_E = ZERO

    do while (p % alive)

      if (check_overlaps) call check_cell_overlap(p)

      ! Calculate microscopic and macroscopic cross sections -- note: if the
      ! material is the same as the last material and the energy of the
      ! particle hasn't changed, we don't need to lookup cross sections again.

      if (p % material /= p % last_material) call calculate_xs(p)

      ! Find the distance to the nearest boundary
      call distance_to_boundary(p, d_boundary, surface_crossed, lattice_crossed)

      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY
      else
        d_collision = -log(prn()) / material_xs % total
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)

      ! Advance particle
      coord => p % coord0
      do while (associated(coord))
        coord % xyz = coord % xyz + distance * coord % uvw
        coord => coord % next
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) &
           call score_tracklength_tally(p, distance)

      ! Score track-length estimate of k-eff
!$omp critical
      global_tallies(K_TRACKLENGTH) % value = &
           global_tallies(K_TRACKLENGTH) % value + p % wgt * distance * &
           material_xs % nu_fission
!$omp end critical

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        last_cell = p % coord % cell
        p % coord % cell = NONE
        if (lattice_crossed /= NONE) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(p, lattice_crossed)
          p % event = EVENT_LATTICE
        else
          ! Particle crosses surface
          p % surface = surface_crossed
          call cross_surface(p, last_cell)
          p % event = EVENT_SURFACE
        end if
      else
        ! ====================================================================
        ! PARTICLE HAS COLLISION

        ! Score collision estimate of keff
!$omp critical
        global_tallies(K_COLLISION) % value = &
             global_tallies(K_COLLISION) % value + p % wgt * &
             material_xs % nu_fission / material_xs % total
!$omp end critical

        ! score surface current tallies -- this has to be done before the collision
        ! since the direction of the particle will change and we need to use the
        ! pre-collision direction to figure out what mesh surfaces were crossed

        if (active_current_tallies % size() > 0) call score_surface_current(p)

        ! Clear surface component
        p % surface = NONE

        call collision(p)

        ! Score collision estimator tallies -- this is done after a collision
        ! has occurred rather than before because we need information on the
        ! outgoing energy for any tallies with an outgoing energy filter

        if (active_analog_tallies % size() > 0) call score_analog_tally(p)

        ! Reset banked weight during collision
        p % n_bank   = 0
        p % wgt_bank = ZERO

        ! Reset fission logical
        p % fission = .false.

        ! Save coordinates for tallying purposes
        p % last_xyz = p % coord0 % xyz

        ! Set last material to none since cross sections will need to be
        ! re-evaluated
        p % last_material = NONE

        ! Set all uvws to base level -- right now, after a collision, only the
        ! base level uvws are changed
        coord => p % coord0
        do while(associated(coord % next))
          if (coord % next % rotated) then
            ! If next level is rotated, apply rotation matrix
            coord % next % uvw = matmul(cells(coord % cell) % &
                 rotation, coord % uvw)
          else
            ! Otherwise, copy this level's direction
            coord % next % uvw = coord % uvw
          end if

          ! Advance coordinate level
          coord => coord % next
        end do
      end if

      ! If particle has too many events, display warning and kill it
      n_event = n_event + 1
      if (n_event == MAX_EVENTS) then
        message = "Particle " // trim(to_str(p%id)) // " underwent maximum &
             &number of events."
        call warning()
        p % alive = .false.
      end if

    end do

  end subroutine transport

!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

  subroutine cross_surface(p, last_cell)

    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: last_cell  ! last cell particle was in

    real(8) :: x         ! x-x0 for sphere
    real(8) :: y         ! y-y0 for sphere
    real(8) :: z         ! z-z0 for sphere
    real(8) :: R         ! radius of sphere
    real(8) :: u         ! x-component of direction
    real(8) :: v         ! y-component of direction
    real(8) :: w         ! z-component of direction
    real(8) :: n1        ! x-component of surface normal
    real(8) :: n2        ! y-component of surface normal
    real(8) :: n3        ! z-component of surface normal
    real(8) :: dot_prod  ! dot product of direction and normal
    real(8) :: norm      ! "norm" of surface normal
    integer :: i_surface ! index in surfaces
    logical :: found     ! particle found in universe?
    type(Surface), pointer, save :: surf => null()
!$omp threadprivate(surf)

    i_surface = abs(p % surface)
    surf => surfaces(i_surface)
    if (verbosity >= 10 .or. trace) then
      message = "    Crossing surface " // trim(to_str(surf % id))
      call write_message()
    end if

    if (surf % bc == BC_VACUUM .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE LEAKS OUT OF PROBLEM

      ! Kill particle
      p % alive = .false.

      ! Score any surface current tallies -- note that the particle is moved
      ! forward slightly so that if the mesh boundary is on the surface, it is
      ! still processed

      if (active_current_tallies % size() > 0) then
        ! TODO: Find a better solution to score surface currents than
        ! physically moving the particle forward slightly

        p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
        call score_surface_current(p)
      end if

      ! Score to global leakage tally
      if (tallies_on) then
!$omp critical
        global_tallies(LEAKAGE) % value = &
           global_tallies(LEAKAGE) % value + p % wgt
!$omp end critical
      end if

      ! Display message
      if (verbosity >= 10 .or. trace) then
        message = "    Leaked out of surface " // trim(to_str(surf % id))
        call write_message()
      end if
      return

    elseif (surf % bc == BC_REFLECT .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE REFLECTS FROM SURFACE

      ! Do not handle reflective boundary conditions on lower universes
      if (.not. associated(p % coord, p % coord0)) then
        message = "Cannot reflect particle " // trim(to_str(p % id)) // &
             " off surface in a lower universe."
        call handle_lost_particle(p)
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing in coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        p % coord0 % xyz = p % coord0 % xyz - TINY_BIT * p % coord0 % uvw
        call score_surface_current(p)
        p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
      end if

      ! Copy particle's direction cosines
      u = p % coord0 % uvw(1)
      v = p % coord0 % uvw(2)
      w = p % coord0 % uvw(3)

      select case (surf%type)
      case (SURF_PX)
        u = -u

      case (SURF_PY)
        v = -v

      case (SURF_PZ)
        w = -w

      case (SURF_PLANE)
        ! Find surface coefficients and norm of vector normal to surface
        n1 = surf % coeffs(1)
        n2 = surf % coeffs(2)
        n3 = surf % coeffs(3)
        norm = n1*n1 + n2*n2 + n3*n3
        dot_prod = u*n1 + v*n2 + w*n3

        ! Reflect direction according to normal
        u = u - 2*dot_prod*n1/norm
        v = v - 2*dot_prod*n2/norm
        w = w - 2*dot_prod*n3/norm

      case (SURF_CYL_X)
        ! Find y-y0, z-z0 and dot product of direction and surface normal
        y = p % coord0 % xyz(2) - surf % coeffs(1)
        z = p % coord0 % xyz(3) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = v*y + w*z

        ! Reflect direction according to normal
        v = v - 2*dot_prod*y/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CYL_Y)
        ! Find x-x0, z-z0 and dot product of direction and surface normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        z = p % coord0 % xyz(3) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = u*x + w*z

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CYL_Z)
        ! Find x-x0, y-y0 and dot product of direction and surface normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = u*x + v*y

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        v = v - 2*dot_prod*y/(R*R)

      case (SURF_SPHERE)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = u*x + v*y + w*z

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        v = v - 2*dot_prod*y/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CONE_X)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (v*y + w*z - R*u*x)/((R + ONE)*R*x*x)

        ! Reflect direction according to normal
        u = u + 2*dot_prod*R*x
        v = v - 2*dot_prod*y
        w = w - 2*dot_prod*z

      case (SURF_CONE_Y)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (u*x + w*z - R*v*y)/((R + ONE)*R*y*y)

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x
        v = v + 2*dot_prod*R*y
        w = w - 2*dot_prod*z

      case (SURF_CONE_Z)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (u*x + v*y - R*w*z)/((R + ONE)*R*z*z)

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x
        v = v - 2*dot_prod*y
        w = w + 2*dot_prod*R*z

      case default
        message = "Reflection not supported for surface " // &
             trim(to_str(surf % id))
        call fatal_error()
      end select

      ! Set new particle direction
      norm = sqrt(u*u + v*v + w*w)
      p % coord0 % uvw = [u, v, w] / norm

      ! Reassign particle's cell and surface
      p % coord0 % cell = last_cell
      p % surface = -p % surface

      ! If a reflective surface is coincident with a lattice or universe
      ! boundary, it is necessary to redetermine the particle's coordinates in
      ! the lower universes.

      if (associated(p % coord0 % next)) then
        call deallocate_coord(p % coord0 % next)
        call find_cell(p, found)
        if (.not. found) then
          message = "Couldn't find particle after reflecting from surface."
          call handle_lost_particle(p)
          return
        end if
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        message = "    Reflected from surface " // trim(to_str(surf%id))
        call write_message()
      end if
      return
    end if

    ! ==========================================================================
    ! SEARCH NEIGHBOR LISTS FOR NEXT CELL

    if (p % surface > 0 .and. allocated(surf % neighbor_pos)) then
      ! If coming from negative side of surface, search all the neighboring
      ! cells on the positive side

      call find_cell(p, found, surf % neighbor_pos)
      if (found) return

    elseif (p % surface < 0  .and. allocated(surf % neighbor_neg)) then
      ! If coming from positive side of surface, search all the neighboring
      ! cells on the negative side

      call find_cell(p, found, surf % neighbor_neg)
      if (found) return

    end if

    ! ==========================================================================
    ! COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

    ! Remove lower coordinate levels and assignment of surface
    p % surface = NONE
    p % coord => p % coord0
    call deallocate_coord(p % coord % next)
    call find_cell(p, found)

    if (run_mode /= MODE_PLOTTING .and. (.not. found)) then
      ! If a cell is still not found, there are two possible causes: 1) there is
      ! a void in the model, and 2) the particle hit a surface at a tangent. If
      ! the particle is really traveling tangent to a surface, if we move it
      ! forward a tiny bit it should fix the problem.

      p % coord => p % coord0
      call deallocate_coord(p % coord % next)
      p % coord % xyz = p % coord % xyz + TINY_BIT * p % coord % uvw
      call find_cell(p, found)

      ! Couldn't find next cell anywhere! This probably means there is an actual
      ! undefined region in the geometry.

      if (.not. found) then
        message = "After particle " // trim(to_str(p % id)) // " crossed surface " &
             // trim(to_str(surfaces(i_surface) % id)) // " it could not be &
             &located in any cell and it did not leak."
        call handle_lost_particle(p)
        return
      end if
    end if
       
  end subroutine cross_surface

end module tracking
