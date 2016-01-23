module dd_tracking

  use constants
  use cross_section,     only: calculate_xs
  use dd_header,         only: DomainDecomType
  use error,             only: fatal_error
  use global,            only: domain_decomp, surfaces, rank, verbosity, trace
  use particle_header,   only: Particle
  use mesh,              only: get_mesh_bin
  use output,            only: write_message
  use random_lcg,        only: prn_seed
  use random_lcg_header, only: N_STREAMS
  use string,            only: to_str

  implicit none
  public

contains

!===============================================================================
! CHECK_DOMAIN_BOUNDARY_CROSSING checks if a particle would cross a domain
! boundary and should be communicated to the process that tracks another domain.
! Care must be taken to account for reflective boundary conditions and lattice
! translations, which might be coincident with domain domain boundaries.
!===============================================================================

  subroutine check_domain_boundary_crossing(d_dd_mesh, &
       d_collision, d_boundary, lattice_translation, surface_crossed, &
       boundary_crossed)

    real(8), intent(in)  :: d_dd_mesh              ! distance to DD boundary
    real(8), intent(in)  :: d_collision            ! distance to collision
    real(8), intent(in)  :: d_boundary             ! distance to cell boundary
    integer, intent(in)  :: lattice_translation(3) ! lattice translation vector
    integer, intent(in)  :: surface_crossed        ! surface which particle is on
    logical, intent(out) :: boundary_crossed

    boundary_crossed = .false.

    ! Check for coincidence with a boundary condition - in this case we
    ! don't need to communicate the particle.  Here we rely on lattice
    ! boundaries NOT being selected by distance_to_boundary when they are
    ! coincident with boundary condition surfaces.
    if (any(lattice_translation /= 0)) then ! communicate if lattice trans
      ! Exit the particle tracking loop without killing the particle
      boundary_crossed = .true.
    elseif (surface_crossed /= None) then
      if(.not. (surfaces(abs(surface_crossed))%obj % bc /= BC_TRANSMIT .and. &
               d_collision > d_boundary .and. &
               abs(d_dd_mesh - d_boundary) < FP_COINCIDENT)) then
        ! Exit the particle tracking loop without killing the particle
        boundary_crossed = .true.
      end if
    end if

  end subroutine check_domain_boundary_crossing

!===============================================================================
! CROSS_DOMAIN_BOUNDARY determines which domain a particle will scatter to,
! stores the position and direction of the particle, and sets the outscatter
! flag for sending this particle later.  Note that the particle is not moved and
! tallies are not recorded. After transfering the particle to the right domain
! (which maybe NOT a neighbor!), the particle will be recovered and tracked.
! We need to use the same distance to the next point for reproducibility. The
! accumulated distance to domain boundary is stored to determine whether the
! particle enters into right domain.
!===============================================================================

  subroutine cross_domain_boundary(p, dd, tracking_dist, flying_dist)

    type(Particle), intent(inout) :: p
    type(DomainDecomType), intent(inout)  :: dd
    real(8), intent(in)           :: tracking_dist ! distance p needs to travel
    real(8), intent(in)           :: flying_dist   ! distance p traveled already

    real(8) :: xyz(3)
    integer :: to_meshbin  ! domain meshbin the particle is traveling to
    integer :: to_bin      ! local relative bin the particle is traveling to

    ! Calculate current point and calculate the bin in the DD mesh
    xyz = p % coord(1) % xyz + (TINY_BIT + flying_dist)* p % coord(1) % uvw
    call get_mesh_bin(dd % mesh, xyz, to_meshbin)

    ! Check for particle leaking out of domain mesh - this is a user input error
    if (to_meshbin == NO_BIN_FOUND) then
      if (.not. dd % allow_truncation)  then
        call fatal_error("For non-truncated DD mode, particle " // &
                  trim(to_str(p % id)) // " leaked out of DD mesh at (" // &
                  trim(to_str(p % coord(1) % xyz(1))) // ", " // &
                  trim(to_str(p % coord(1) % xyz(2))) // ", " // &
                  trim(to_str(p % coord(1) % xyz(3))) // ") on rank " // &
                  trim(to_str(rank)) // ". Does the DD mesh " // &
                  "completely envelope the defined geometry?")
      end if
      return
    end if

    ! Check for a bad determination of a change - this would be a bug
    if (to_meshbin == dd % meshbin) then
      call fatal_error("Can't determine which domain to send particle " // &
           "on rank " // trim(to_str(rank)) // ". Particle "// &
                  trim(to_str(p % id)) // " at (" // &
                  trim(to_str(p % coord(1) % xyz(1))) // ", " // &
                  trim(to_str(p % coord(1) % xyz(2))) // ", " // &
                  trim(to_str(p % coord(1) % xyz(3))) // ") on rank " // &
                  trim(to_str(rank)) // ". (to_meshbin, dd % meshbin): (" // &
                  trim(to_str(to_meshbin)) // ", " // &
                  trim(to_str(dd % meshbin)) // ").  prn_seed: " // &
                  trim(to_str(prn_seed(1))))
    end if

    if (verbosity >= 10 .or. trace) then
      call write_message('Scatter from domain ' // &
                trim(to_str(dd % meshbin)) // &
                ' to ' // trim(to_str(to_meshbin)) // ' pid = ' // &
                trim(to_str(p % id)))
    end if

    ! Check if we're scattering further than a direct neighbor - bug if so
    if (.not. dd % bins_dict % has_key(to_meshbin)) then
      call fatal_error("Not transferring to direct neighbor! From domain " // &
                trim(to_str(dd % meshbin)) // ' to ' // &
                trim(to_str(to_meshbin)) // ' pid = ' // trim(to_str(p % id)))
    end if

    ! Convert destination domain meshbin to relative local bin
    to_bin = dd % bins_dict % get_key(to_meshbin)

    ! Note where this particle will be transmitted (after all particles run)
    p % outscatter_destination = to_bin

    ! Increment count of how many are going to the local neighbor
    dd % n_scatters_local(to_bin) = dd % n_scatters_local(to_bin) + 1

    ! Save the transport info needed to restart the particle in the new domain
    p % stored_xyz      = p % coord(1) % xyz
    p % stored_uvw      = p % coord(1) % uvw
    p % stored_distance = tracking_dist
    p % fly_dd_distance = flying_dist
    p % prn_seed        = prn_seed

  end subroutine cross_domain_boundary


!===============================================================================
! RECALC_INITIAL_XS recalculates the inital cross sections for a particle using
! a stored random number seed. For DD runs, if we normally wouldn't have to
! recalculate the cross section after a scatter then we need to make sure that
! we recalculate it before starting transport with the same random number seed
! so we get the same thing as we would have gotten if we tracked the particle to
! completion without transporting it across domains.  This is needed entirely
! because URR ptables use a random number from the stream.
!===============================================================================

  subroutine recalc_initial_xs(p)

    type(Particle), intent(inout) :: p

    integer(8) :: tmp_seed(N_STREAMS) ! Temporary variable to hold prn_seed

    if (p % material /= NONE) then  ! .and. p % material == p % last_material
      tmp_seed = prn_seed
      prn_seed = p % xs_seed
      call calculate_xs(p)
      prn_seed = tmp_seed
    end if

  end subroutine recalc_initial_xs

end module dd_tracking
