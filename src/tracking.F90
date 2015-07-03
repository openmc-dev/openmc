module tracking

  use constants
  use cross_section,     only: calculate_xs
  use dd_header,         only: dd_type
  use dd_tracking,       only: check_domain_boundary_crossing, recalc_initial_xs
  use error,             only: fatal_error, warning
  use geometry,          only: find_cell, distance_to_boundary, cross_surface, &
                               cross_lattice, check_cell_overlap, &
                               distance_to_mesh_surface
  use geometry_header,   only: Universe, BASE_UNIVERSE
  use global
  use material_header
  use mesh,              only: get_mesh_bin
  use output,            only: write_message
  use particle_header,   only: LocalCoord, Particle
  use physics,           only: collision
  use random_lcg,        only: prn, prn_seed
  use string,            only: to_str
  use tally,             only: score_analog_tally, score_tracklength_tally, &
                               score_surface_current
  use track_output,      only: initialize_particle_track, &
                               write_particle_track, finalize_particle_track
  implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), intent(inout) :: p

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
    type(LocalCoord), pointer :: coord

    ! DD debugging vars
    integer :: meshbin
    integer(8) :: starting_seed
    integer(8) :: debug1 = 5464171501723750348_8
    integer(8) :: debug2 = 331860069508688933_8
    integer(8) :: debug3 = 0_8
    integer(8) :: debug4 = 0_8

    starting_seed = prn_seed(1)
    !print *, 'starting', starting_seed
    if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
        print *,'starting particle', starting_seed, p % id, p % new_particle
    ! Display message if high verbosity or trace is on
    if (verbosity >= 9 .or. trace) then
      call write_message("Simulating Particle " // trim(to_str(p % id)))
    end if

    ! If the cell hasn't been determined based on the particle's location,
    ! initiate a search for the current cell
    if (p % coord % cell == NONE) then
      call find_cell(p, found_cell)

      ! Particle couldn't be located
      if (.not. found_cell) then
        call fatal_error("Could not locate particle " // trim(to_str(p % id)))
      end if

      ! set birth cell attribute
      p % cell_born = p % coord % cell
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

    ! Make sure we start with the proper XS for DD runs
    if (dd_run) call recalc_initial_xs(p)

    do while (p % alive)
      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
          print *, prn_seed(1), p % coord0 % xyz
!      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
!        then
!          call get_mesh_bin(domain_decomp % mesh, p % coord0 % xyz, meshbin)
!          print *, 'meshbin',meshbin
!        end if

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

      if (check_overlaps) call check_cell_overlap(p)

!      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
!          print *, prn_seed(1), p % material /= p % last_material, p % inst /= p % last_inst, p % coord0 % xyz(2)

      ! Calculate microscopic and macroscopic cross sections -- note: if the
      ! material is the same as the last material and the energy of the
      ! particle hasn't changed, we don't need to lookup cross sections again.
      if (p % material /= p % last_material &
            & .or. p % inst /= p % last_inst) then

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

      ! Find the distance to the nearest boundary
!      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
!          print *, prn_seed(1), 'd2b',p % coord0 % xyz
      call distance_to_boundary(p, d_boundary, surface_crossed, &
           &lattice_translation)

      ! Sample a distance to collision
      if (material_xs % total == ZERO) then
        d_collision = INFINITY
      else
        if (p % stored_distance > ZERO) then
          d_collision = p % stored_distance
          p % stored_distance = ZERO
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
        distance = min(distance, d_dd_mesh)
      end if

      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
          print *, prn_seed(1), 'distances', d_boundary, d_dd_mesh, d_collision
!      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
!          print *, prn_seed(1), 'direction', p % coord0 % uvw

      ! Advance particle
      coord => p % coord0
      do while (associated(coord))
        coord % xyz = coord % xyz + distance * coord % uvw
        coord => coord % next
      end do
!      if (starting_seed == debug1 .or. starting_seed == debug2 .or. starting_seed == debug3 .or. starting_seed == debug4) &
!          print *, prn_seed(1), 'advanced', p % coord0 % xyz

      ! Score track-length tallies
      if (active_tracklength_tallies % size() > 0) &
           call score_tracklength_tally(p, distance)

      ! Score track-length estimate of k-eff
!$omp atomic
      global_tallies(K_TRACKLENGTH) % value = &
           global_tallies(K_TRACKLENGTH) % value + p % wgt * distance * &
           material_xs % nu_fission

      ! Check for domain boundary crossing
      if (dd_run) then
        call check_domain_boundary_crossing(domain_decomp, p, distance, &
            d_dd_mesh, d_collision, d_boundary, lattice_translation, &
            surface_crossed, dd_boundary_crossed)
        ! If the particle crosses a domain boundary, stop tracking it
        if (dd_boundary_crossed) exit
      end if

      ! Check for collisions and surface crossings
      if (d_collision > d_boundary) then
        ! ====================================================================
        ! PARTICLE CROSSES SURFACE

        last_cell = p % coord % cell
        p % coord % cell = NONE
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
!$omp atomic
        global_tallies(K_COLLISION) % value = &
             global_tallies(K_COLLISION) % value + p % wgt * &
             material_xs % nu_fission / material_xs % total

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
                 rotation_matrix, coord % uvw)
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
        if (master) call warning("Particle " // trim(to_str(p%id)) &
             &// " underwent maximum number of events.")
        p % alive = .false.
      end if

    end do

    ! Add number of events to DD object
    if (dd_run .and. domain_decomp % count_interactions) then
      domain_decomp % interaction_count = domain_decomp % interaction_count + n_event
    end if

!    if (current_batch == 1) then
!      if (.not. p % alive) print *, prn_seed(1), starting_seed, p % coord0 % xyz(1), rank, p % id
!    end if

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

end module tracking
