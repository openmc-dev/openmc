module dd_testing_setup

  use constants
  use dd_header,        only: dd_type
  use error,            only: warning
  use global,           only: n_procs, rank, message, n_particles
  use string,           only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none
  public

contains

!===============================================================================
! CHECK_PROCS makes sure we have the right number of processors
!===============================================================================

  function check_procs(n_procs) result(skip)
  
    integer, intent(in) :: n_procs

    logical :: skip
    skip = .false.
  
#ifdef MPI
    if (.not. n_procs == 5) then
      message = "Skipping test: must be run with MPI " // &
                trim(to_str(n_procs)) // " procs"
      call warning()
      skip = .true.
    end if
#else
    message = "Skipping test: requires MPI"
    call warning()
    skip = .true.
#endif

  end function check_procs

!===============================================================================
! DD_SIMPLE_FOUR_DOMAINS hardcodes in a simple 4-domain xml input
!===============================================================================

  subroutine dd_simple_four_domains(dd)
  
    type(dd_type), intent(inout) :: dd
  
    ! Mimic everything done in input_xml for DD setup
    allocate(dd % mesh)
    allocate(dd % mesh % dimension(3))
    allocate(dd % mesh % lower_left(3))
    allocate(dd % mesh % upper_right(3))
    allocate(dd % mesh % width(3))
    dd % mesh % n_dimension = 3

    dd % mesh % dimension = (/ 2, 2, 1/)
    dd % n_domains = 4
    dd % mesh % lower_left = (/ -2, -2, -1/)
    dd % mesh % upper_right = (/ 2, 2, 1/)
    dd % mesh % width = (/ 2, 2, 2/)
  
    allocate(dd % domain_load_dist(4))
    dd % domain_load_dist = (/ 2.0, 1.0, 1.0, 1.0/)

    ! Set n_particles so buffers are properly allocated
    n_particles = 97 ! we will have 55 transferring

  end subroutine dd_simple_four_domains

!===============================================================================
! DD_SIMPLE_FOUR_DOMAIN_SCATTERS hardcodes in a few particles to transfer, based
! on dd_simple_four_domains
!===============================================================================

  subroutine dd_simple_four_domain_scatters(dd)
  
    type(dd_type), intent(inout) :: dd
    integer :: bin

    ! Initialize particle buffer
    dd % particle_buffer % outscatter_destination = NO_OUTSCATTER

    dd % n_scatters_local = 0
    select case(rank)
      case (0)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  7 !  7 from 1 -> 2
        dd % particle_buffer(2) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(9) % outscatter_destination = bin
        dd % particle_buffer(11) % outscatter_destination = bin
        dd % particle_buffer(12) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  5 !  5 from 1 -> 3
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(10) % outscatter_destination = bin
        dd % particle_buffer(14) % outscatter_destination = bin
      case (1)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  2 !  2 from 1 -> 2
        dd % particle_buffer(2) % outscatter_destination = bin
        dd % particle_buffer(13) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  9 !  9 from 1 -> 3
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(9) % outscatter_destination = bin
        dd % particle_buffer(12) % outscatter_destination = bin
        dd % particle_buffer(14) % outscatter_destination = bin
      case (2)
        bin = dd % bins_dict % get_key(1)
        dd % n_scatters_local(bin) = 12 ! 12 from 2 -> 1
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(6) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(13) % outscatter_destination = bin
        dd % particle_buffer(15) % outscatter_destination = bin
        dd % particle_buffer(17) % outscatter_destination = bin
        dd % particle_buffer(18) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(21) % outscatter_destination = bin
        dd % particle_buffer(22) % outscatter_destination = bin
        dd % particle_buffer(24) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(4)
        dd % n_scatters_local(bin) =  3 !  3 from 2 -> 4
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(16) % outscatter_destination = bin
        dd % particle_buffer(20) % outscatter_destination = bin
      case (3)
        bin = dd % bins_dict % get_key(1)
        dd % n_scatters_local(bin) =  7 !  7 from 3 -> 1
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(26) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(4)
        dd % n_scatters_local(bin) =  1 !  1 from 3 -> 4
        dd % particle_buffer(14) % outscatter_destination = bin
      case (4)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  0 !  0 from 4 -> 2
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  9 !  9 from 4 -> 3
        dd % particle_buffer(14) % outscatter_destination = bin
        dd % particle_buffer(15) % outscatter_destination = bin
        dd % particle_buffer(16) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(21) % outscatter_destination = bin
        dd % particle_buffer(23) % outscatter_destination = bin
        dd % particle_buffer(24) % outscatter_destination = bin
        dd % particle_buffer(25) % outscatter_destination = bin
        dd % particle_buffer(27) % outscatter_destination = bin
    end select

  end subroutine dd_simple_four_domain_scatters

end module dd_testing_setup
