module tracking

  use cross_section,   only: calculate_xs
  use error,           only: fatal_error, warning
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, check_cell_overlap
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use output,          only: write_message
  use particle_header, only: LocalCoord, Particle
  use physics,         only: collision
  use random_lcg,      only: prn
  use string,          only: to_str
  use tally,           only: score_analog_tally, score_tracklength_tally, &
                             score_surface_current
  use track_output,    only: initialize_particle_track, write_particle_track, &
                             finalize_particle_track

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

    ! Prepare to write out particle track.
    if (p % write_track) then
      call initialize_particle_track()
    endif

    do while (p % alive)

      ! Write particle track.
      if (p % write_track) call write_particle_track(p)

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

    ! Finish particle track output.
    if (p % write_track) then
      call write_particle_track(p)
      call finalize_particle_track(p)
    endif

  end subroutine transport

end module tracking
