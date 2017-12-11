module tracking

  use constants
  use cross_section,      only: calculate_xs
  use error,              only: fatal_error, warning, write_message
  use geometry_header,    only: cells
  use geometry,           only: find_cell, distance_to_boundary, cross_lattice, &
                                check_cell_overlap
  use message_passing
  use mgxs_header
  use nuclide_header
  use particle_header,    only: LocalCoord, Particle
  use physics,            only: collision
  use physics_mg,         only: collision_mg
  use random_lcg,         only: prn
  use settings
  use simulation_header
  use string,             only: to_str
  use surface_header
  use tally_header
  use tally,              only: score_analog_tally, score_tracklength_tally, &
                                score_collision_tally, score_surface_current, &
                                score_track_derivative, score_surface_tally, &
                                score_collision_derivative, zero_flux_derivs
  use track_output,       only: initialize_particle_track, write_particle_track, &
                                add_particle_track, finalize_particle_track

  implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), intent(inout) :: p

    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! sampled distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      call write_message("Simulating Particle " // trim(to_str(p % id)))
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
!$omp atomic
    total_weight = total_weight + p % wgt

    ! Force calculation of cross-sections by setting last energy to zero
    if (run_CE) then
      micro_xs % last_E = ZERO
    end if

    ! Prepare to write out particle track.
    if (p % write_track) then
      call initialize_particle_track()
    endif

    ! Every particle starts with no accumulated flux derivative.
    if (active_tallies % size() > 0) call zero_flux_derivs()

    EVENT_LOOP: do
      ! Store pre-collision particle properties
      p % last_wgt = p % wgt
      p % last_E   = p % E
      p % last_uvw = p % coord(1) % uvw
      p % last_xyz = p % coord(1) % xyz

      ! If the cell hasn't been determined based on the particle's location,
      ! initiate a search for the current cell. This generally happens at the
      ! beginning of the history and again for any secondary particles
      if (p % coord(p % n_coord) % cell == NONE) then
        call find_cell(p, found_cell)
        if (.not. found_cell) then
          call fatal_error("Could not locate particle " // trim(to_str(p % id)))
        end if

        ! set birth cell attribute
        if (p % cell_born == NONE) p % cell_born = p % coord(p % n_coord) % cell
      end if

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

      if (check_overlaps) call check_cell_overlap(p)

      ! Calculate microscopic and macroscopic cross sections
      if (run_CE) then
        ! If the material is the same as the last material and the temperature
        ! hasn't changed, we don't need to lookup cross sections again.
        if (p % material /= p % last_material .or. &
             p % sqrtkT /= p % last_sqrtkT) call calculate_xs(p)
      else
        ! Since the MGXS can be angle dependent, this needs to be done
        ! After every collision for the MGXS mode
        if (p % material /= MATERIAL_VOID) then
          ! Update the temperature index
          call macro_xs(p % material) % obj % find_temperature(p % sqrtkT)
          ! Get the data
          call macro_xs(p % material) % obj % calculate_xs(p % g, &
               p % coord(p % n_coord) % uvw, material_xs)
        else
          material_xs % total      = ZERO
          material_xs % absorption = ZERO
          material_xs % nu_fission = ZERO
        end if

        ! Finally, update the particle group while we have already checked for
        ! if multi-group
        p % last_g = p % g
      end if

      ! Find the distance to the nearest boundary
      call distance_to_boundary(p, d_boundary, surface_crossed, &
           lattice_translation, next_level)

      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY
      else
        d_collision = -log(prn()) / material_xs % total
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)

      ! Advance particle
      do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) then
        call score_tracklength_tally(p, distance)
      end if

      ! Score track-length estimate of k-eff
      if (run_mode == MODE_EIGENVALUE) then
        global_tally_tracklength = global_tally_tracklength + p % wgt * &
             distance * material_xs % nu_fission
      end if

      ! Score flux derivative accumulators for differential tallies.
      if (active_tallies % size() > 0) call score_track_derivative(p, distance)

      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        if (next_level > 0) p % n_coord = next_level

        ! Saving previous cell data
        do j = 1, p % n_coord
          p % last_cell(j) = p % coord(j) % cell
        end do
        p % last_n_coord = p % n_coord

        p % coord(p % n_coord) % cell = NONE
        if (any(lattice_translation /= 0)) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(p, lattice_translation)
          p % event = EVENT_LATTICE
        else
          ! Particle crosses surface
          p % surface = surface_crossed

          call cross_surface(p)
          p % event = EVENT_SURFACE
        end if
        ! Score cell to cell partial currents
        if(active_surface_tallies % size() > 0) call score_surface_tally(p)
      else
        ! ====================================================================
        ! PARTICLE HAS COLLISION

        ! Score collision estimate of keff
        if (run_mode == MODE_EIGENVALUE) then
          global_tally_collision = global_tally_collision + p % wgt * &
               material_xs % nu_fission / material_xs % total
        end if

        ! score surface current tallies -- this has to be done before the collision
        ! since the direction of the particle will change and we need to use the
        ! pre-collision direction to figure out what mesh surfaces were crossed

        if (active_current_tallies % size() > 0) call score_surface_current(p)

        ! Clear surface component
        p % surface = NONE

        if (run_CE) then
          call collision(p)
        else
          call collision_mg(p)
        end if

        ! Score collision estimator tallies -- this is done after a collision
        ! has occurred rather than before because we need information on the
        ! outgoing energy for any tallies with an outgoing energy filter
        if (active_collision_tallies % size() > 0) call score_collision_tally(p)
        if (active_analog_tallies % size() > 0) call score_analog_tally(p)

        ! Reset banked weight during collision
        p % n_bank   = 0
        p % wgt_bank = ZERO
        p % n_delayed_bank(:) = 0

        ! Reset fission logical
        p % fission = .false.

        ! Save coordinates for tallying purposes
        p % last_xyz_current = p % coord(1) % xyz

        ! Set last material to none since cross sections will need to be
        ! re-evaluated
        p % last_material = NONE

        ! Set all uvws to base level -- right now, after a collision, only the
        ! base level uvws are changed
        do j = 1, p % n_coord - 1
          if (p % coord(j + 1) % rotated) then
            ! If next level is rotated, apply rotation matrix
            p % coord(j + 1) % uvw = matmul(cells(p % coord(j) % cell) % &
                 rotation_matrix, p % coord(j) % uvw)
          else
            ! Otherwise, copy this level's direction
            p % coord(j + 1) % uvw = p % coord(j) % uvw
          end if
        end do

        ! Score flux derivative accumulators for differential tallies.
        if (active_tallies % size() > 0) call score_collision_derivative(p)
      end if

      ! If particle has too many events, display warning and kill it
      n_event = n_event + 1
      if (n_event == MAX_EVENTS) then
        if (master) call warning("Particle " // trim(to_str(p%id)) &
             &// " underwent maximum number of events.")
        p % alive = .false.
      end if

      ! Check for secondary particles if this particle is dead
      if (.not. p % alive) then
        if (p % n_secondary > 0) then
          call p % initialize_from_source(p % secondary_bank(p % n_secondary), &
                                          run_CE, energy_bin_avg)
          p % n_secondary = p % n_secondary - 1
          n_event = 0

          ! Enter new particle in particle track file
          if (p % write_track) call add_particle_track()
        else
          exit EVENT_LOOP
        end if
      end if
    end do EVENT_LOOP

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

  subroutine cross_surface(p)
    type(Particle), intent(inout) :: p

    real(8) :: u          ! x-component of direction
    real(8) :: v          ! y-component of direction
    real(8) :: w          ! z-component of direction
    real(8) :: norm       ! "norm" of surface normal
    real(8) :: d          ! distance between point and plane
    real(8) :: xyz(3)     ! Saved global coordinate
    integer :: i_surface  ! index in surfaces
    logical :: rotational ! if rotational periodic BC applied
    logical :: found      ! particle found in universe?
    class(Surface), pointer :: surf

    i_surface = abs(p % surface)
    surf => surfaces(i_surface)%obj
    if (verbosity >= 10 .or. trace) then
      call write_message("    Crossing surface " // trim(to_str(surf % id)))
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

        p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
      end if

      ! Score to global leakage tally
      global_tally_leakage = global_tally_leakage + p % wgt

      ! Display message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Leaked out of surface " &
             &// trim(to_str(surf % id)))
      end if
      return

    elseif (surf % bc == BC_REFLECT .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE REFLECTS FROM SURFACE

      ! Do not handle reflective boundary conditions on lower universes
      if (p % n_coord /= 1) then
        call p % mark_as_lost("Cannot reflect particle " &
             // trim(to_str(p % id)) // " off surface in a lower universe.")
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing is coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        xyz = p % coord(1) % xyz
        p % coord(1) % xyz = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
        p % coord(1) % xyz = xyz
      end if

      ! Reflect particle off surface
      call surf%reflect(p%coord(1)%xyz, p%coord(1)%uvw)

      ! Make sure new particle direction is normalized
      u = p%coord(1)%uvw(1)
      v = p%coord(1)%uvw(2)
      w = p%coord(1)%uvw(3)
      norm = sqrt(u*u + v*v + w*w)
      p%coord(1)%uvw(:) = [u, v, w] / norm

      ! Reassign particle's cell and surface
      p % coord(1) % cell = p % last_cell(p % last_n_coord)
      p % surface = -p % surface

      ! If a reflective surface is coincident with a lattice or universe
      ! boundary, it is necessary to redetermine the particle's coordinates in
      ! the lower universes.

      p % n_coord = 1
      call find_cell(p, found)
      if (.not. found) then
        call p % mark_as_lost("Couldn't find particle after reflecting&
             & from surface " // trim(to_str(surf % id)) // ".")
        return
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz_current = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Reflected from surface " &
             &// trim(to_str(surf%id)))
      end if
      return
    elseif (surf % bc == BC_PERIODIC .and. run_mode /= MODE_PLOTTING) then
      ! =======================================================================
      ! PERIODIC BOUNDARY

      ! Do not handle periodic boundary conditions on lower universes
      if (p % n_coord /= 1) then
        call p % mark_as_lost("Cannot transfer particle " &
             // trim(to_str(p % id)) // " across surface in a lower universe.&
             & Boundary conditions must be applied to universe 0.")
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing is coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        xyz = p % coord(1) % xyz
        p % coord(1) % xyz = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
        p % coord(1) % xyz = xyz
      end if

      rotational = .false.
      select type (surf)
      type is (SurfaceXPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceXPlane)
          p % coord(1) % xyz(1) = opposite % x0
        type is (SurfaceYPlane)
          rotational = .true.

          ! Rotate direction
          u = p % coord(1) % uvw(1)
          v = p % coord(1) % uvw(2)
          p % coord(1) % uvw(1) = v
          p % coord(1) % uvw(2) = -u

          ! Rotate position
          p % coord(1) % xyz(1) = surf % x0 + p % coord(1) % xyz(2) - opposite % y0
          p % coord(1) % xyz(2) = opposite % y0
        end select

      type is (SurfaceYPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceYPlane)
          p % coord(1) % xyz(2) = opposite % y0
        type is (SurfaceXPlane)
          rotational = .true.

          ! Rotate direction
          u = p % coord(1) % uvw(1)
          v = p % coord(1) % uvw(2)
          p % coord(1) % uvw(1) = -v
          p % coord(1) % uvw(2) = u

          ! Rotate position
          p % coord(1) % xyz(2) = surf % y0 + p % coord(1) % xyz(1) - opposite % x0
          p % coord(1) % xyz(1) = opposite % x0
        end select

      type is (SurfaceZPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceZPlane)
          p % coord(1) % xyz(3) = opposite % z0
        end select

      type is (SurfacePlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfacePlane)
          ! Get surface normal for opposite plane
          xyz(:) = opposite % normal(p % coord(1) % xyz)

          ! Determine distance to plane
          norm = xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)
          d = opposite % evaluate(p % coord(1) % xyz) / norm

          ! Move particle along normal vector based on distance
          p % coord(1) % xyz(:) = p % coord(1) % xyz(:) - d*xyz
        end select
      end select

      ! Reassign particle's surface
      if (rotational) then
        p % surface = surf % i_periodic
      else
        p % surface = sign(surf % i_periodic, p % surface)
      end if

      ! Figure out what cell particle is in now
      p % n_coord = 1
      call find_cell(p, found)
      if (.not. found) then
        call p % mark_as_lost("Couldn't find particle after hitting &
             &periodic boundary on surface " // trim(to_str(surf % id)) // ".")
        return
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz_current = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Hit periodic boundary on surface " &
             // trim(to_str(surf%id)))
      end if
      return
    end if

    ! ==========================================================================
    ! SEARCH NEIGHBOR LISTS FOR NEXT CELL

    if (p % surface > 0 .and. allocated(surf%neighbor_pos)) then
      ! If coming from negative side of surface, search all the neighboring
      ! cells on the positive side

      call find_cell(p, found, surf%neighbor_pos)
      if (found) return

    elseif (p % surface < 0  .and. allocated(surf%neighbor_neg)) then
      ! If coming from positive side of surface, search all the neighboring
      ! cells on the negative side

      call find_cell(p, found, surf%neighbor_neg)
      if (found) return

    end if

    ! ==========================================================================
    ! COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

    ! Remove lower coordinate levels and assignment of surface
    p % surface = NONE
    p % n_coord = 1
    call find_cell(p, found)

    if (run_mode /= MODE_PLOTTING .and. (.not. found)) then
      ! If a cell is still not found, there are two possible causes: 1) there is
      ! a void in the model, and 2) the particle hit a surface at a tangent. If
      ! the particle is really traveling tangent to a surface, if we move it
      ! forward a tiny bit it should fix the problem.

      p % n_coord = 1
      p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
      call find_cell(p, found)

      ! Couldn't find next cell anywhere! This probably means there is an actual
      ! undefined region in the geometry.

      if (.not. found) then
        call p % mark_as_lost("After particle " // trim(to_str(p % id)) &
             // " crossed surface " // trim(to_str(surf % id)) &
             // " it could not be located in any cell and it did not leak.")
        return
      end if
    end if

  end subroutine cross_surface

end module tracking
