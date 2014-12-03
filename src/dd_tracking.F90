module dd_tracking

  use constants
  use dd_header,        only: dd_type
  use error,            only: fatal_error
  use global,           only: domain_decomp, rank, verbosity, trace
  use particle_header,  only: Particle
  use mesh,             only: get_mesh_bin
  use output,           only: write_message
  use random_lcg,       only: prn_seed
  use string,           only: to_str
  
  implicit none
  private
  public :: cross_domain_boundary

contains

!===============================================================================
! CROSS_DOMAIN_BOUNDARY determines which domain a particle will scatter to,
! stores the position and direction of the particle, and sets the outscatter
! flag for sending this particle later
!===============================================================================

  subroutine cross_domain_boundary(p, dd, dist)

    type(Particle), intent(inout) :: p
    type(dd_type), intent(inout)  :: dd
    real(8), intent(in)           :: dist ! distance p still needs to travel
  
    real(8) :: xyz(3)
    integer :: to_meshbin  ! domain meshbin the particle is traveling to
    integer :: to_bin      ! local relative bin the particle is traveling to
  
    ! Advance particle a little and recalculate the bin in the DD mesh
    xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
    call get_mesh_bin(dd % mesh, xyz, to_meshbin)
    
    ! Check for particle leaking out of domain mesh - this is a user input error
    if (to_meshbin == NO_BIN_FOUND) then
      if (.not. dd % allow_truncation)  then
        call fatal_error("For non-truncated DD mode, particle " // &
                  trim(to_str(p % id)) // " leaked out of DD mesh at (" // &
                  trim(to_str(p % coord0 % xyz(1))) // ", " // &
                  trim(to_str(p % coord0 % xyz(2))) // ", " // &
                  trim(to_str(p % coord0 % xyz(3))) // ") on rank " // &
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
                  trim(to_str(p % coord0 % xyz(1))) // ", " // &
                  trim(to_str(p % coord0 % xyz(2))) // ", " // &
                  trim(to_str(p % coord0 % xyz(3))) // ") on rank " // &
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
    p % stored_xyz      = p % coord0 % xyz
    p % stored_uvw      = p % coord0 % uvw
    p % stored_distance = dist
    p % prn_seed        = prn_seed
    
  end subroutine cross_domain_boundary

end module dd_tracking
