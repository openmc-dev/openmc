module tracking

  use constants,       only: MODE_EIGENVALUE
  use cross_section,   only: calculate_xs
  use dd_header,       only: dd_type
  use dd_tracking,    only: recalc_initial_xs, &
                            check_domain_boundary_crossing, &
                            cross_domain_boundary
  use error,           only: fatal_error, warning
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, check_cell_overlap, &
                             distance_to_mesh_surface, calibrate_coord
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use material_header
  use mesh,            only: get_mesh_bin
  use output,          only: write_message
  use particle_header, only: LocalCoord, Particle
  use physics,         only: collision
  use random_lcg,      only: prn, prn_seed
  use string,          only: to_str
  use tally,           only: score_analog_tally, score_tracklength_tally, &
                             score_collision_tally, score_surface_current
  use track_output,    only: initialize_particle_track, write_particle_track, &
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
    integer :: last_cell              ! most recent cell particle was in
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! sampled distance to collision
    real(8) :: d_dd_mesh              ! distance to boundary on the DD mesh
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    logical :: dd_boundary_crossed    ! domain decomposition boundary crossed
    real(8) :: xyz_temp(3)            ! temporary storage of particle position

    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      call write_message("Simulating Particle " // trim(to_str(p % id)))
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! Add paricle's starting weight to count for normalizing tallies later
    if (p % new_particle) then
!$omp atomic
      total_weight = total_weight + p % wgt
    end if

    ! Force calculation of cross-sections by setting last energy to zero
    micro_xs % last_E = ZERO
    micro_xs % last_index_sab = NONE

    ! Prepare to write out particle track.
    if (p % write_track) then
      call initialize_particle_track()
    endif

    ! Make sure we start with the proper XS and coordinates for DD runs
    if (dd_run .and. current_stage > 1) then
      call recalc_initial_xs(p)

      ! Move particle slightly forward to recover particle coordinates
      ! (The lower coordinates will be recalculted by calibrate_coord function)
      xyz_temp = p % coord(1) % xyz
      p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
      call find_cell(p, found_cell)
      if (.not. found_cell) then
        call fatal_error("Could not locate particle after crossing domain.")
      end if
      p % coord(1) % xyz = xyz_temp
    end if

    EVENT_LOOP: do
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

      else
        ! calibrate coordinates to guarantee numerical precision
        ! especially for domain decom
        if (dd_run) call calibrate_coord(p)

      end if

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

      if (check_overlaps) call check_cell_overlap(p)

      ! Calculate microscopic and macroscopic cross sections -- note: if the
      ! material is the same as the last material and the energy of the
      ! particle hasn't changed, we don't need to lookup cross sections again.
      if (p % material /= p % last_material ) then

        if (dd_run) then
          ! NOTE: calculate_xs does not re-calculate XS on a per-nuclide
          ! basis if the energy hasn't changed. For DD reproducibility we'd
          ! need to send a prn seed for each nuclide along with particles as
          ! they cross domain boundaries.  This would just about double the
          ! amount of info that goes with particles across boundaries, so we
          ! disable it for DD runs.  If we care about XS calculation speed more
          ! than network communication load, we might consider communicating
          ! these seeds with particles.  Or, we could comment out these two
          ! lines and lose random number reproducibility.
          micro_xs % last_E = ZERO
          micro_xs % last_index_sab = NONE
        end if
        call calculate_xs(p)
      end if

      ! determine whether the particle is now transferred to right domain
      if (dd_run .and. p % stored_distance > ZERO) then
        xyz_temp = p % coord(1) % xyz
        p % coord(1) % xyz = xyz_temp + p % fly_dd_distance * p % coord(1) % uvw
        call distance_to_mesh_surface(p, domain_decomp % mesh, &
             distance, d_dd_mesh, meshbin=domain_decomp % meshbin)
        p % coord(1) % xyz = xyz_temp
        p % fly_dd_distance = p % fly_dd_distance + d_dd_mesh
        if (p % stored_distance - p % fly_dd_distance > FP_COINCIDENT) then
          ! Prepare particle for communication and stop tracking it
          call cross_domain_boundary(p, domain_decomp, p % stored_distance, &
               p % fly_dd_distance)
          exit
        end if
      end if

      ! Find the distance to the nearest boundary
      call distance_to_boundary(p, d_boundary, surface_crossed, &
           lattice_translation, next_level)

      ! Sample a distance to collision
      if (dd_run .and. p % stored_distance > ZERO) then
        if (abs(p % stored_distance - d_boundary) < 100*FP_COINCIDENT) then
          d_boundary  = p % stored_distance
          d_collision = INFINITY
        else
          d_collision = p % stored_distance
        end if
          p % stored_distance = ZERO
          p % fly_dd_distance = ZERO
      else
        if (material_xs % total == ZERO) then
          d_collision = INFINITY
        else
          d_collision = -log(prn()) / material_xs % total
        end if
      end if

      ! Select smaller of the two distances
      distance = min(d_boundary, d_collision)

      ! Check domain mesh boundary
      if (dd_run) then
        call distance_to_mesh_surface(p, domain_decomp % mesh, &
             distance, d_dd_mesh, meshbin=domain_decomp % meshbin)
      end if

      ! Check for domain boundary crossing
      if (dd_run .and. d_dd_mesh <= distance) then

        ! Determine if we should communicate the particle
        call check_domain_boundary_crossing(d_dd_mesh, d_collision, d_boundary,&
             lattice_translation, surface_crossed, dd_boundary_crossed)

        if (dd_boundary_crossed) then
          ! =================================================================
          ! PARTICLE CROSSES DOMAIN BOUNDARY

          ! Prepare particle for communication and stop tracking it
          call cross_domain_boundary(p, domain_decomp, distance, d_dd_mesh)
          exit

        end if

      end if

      ! Advance particle
      do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
      end do

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) &
           call score_tracklength_tally(p, distance)

      ! Score track-length estimate of k-eff
      if (run_mode == MODE_EIGENVALUE) then
        global_tally_tracklength = global_tally_tracklength + p % wgt * &
             distance * material_xs % nu_fission
      end if

      ! Check for collisions and surface crossings
      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        if (next_level > 0) p % n_coord = next_level
        last_cell = p % coord(p % n_coord) % cell
        p % coord(p % n_coord) % cell = NONE
        if (any(lattice_translation /= 0)) then
          ! Particle crosses lattice boundary
          p % surface = NONE
          call cross_lattice(p, lattice_translation)
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

        call collision(p)

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
        p % last_xyz = p % coord(1) % xyz

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
          call p % initialize_from_source(p % secondary_bank(p % n_secondary))
          p % n_secondary = p % n_secondary - 1

          ! Enter new particle in particle track file
          if (p % write_track) call add_particle_track()
        else
          exit EVENT_LOOP
        end if
      end if
    end do EVENT_LOOP

    ! Add number of events to DD object
    if (dd_run .and. domain_decomp % count_interactions) then
      domain_decomp % interaction_count = domain_decomp % interaction_count + n_event
    end if

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

end module tracking
