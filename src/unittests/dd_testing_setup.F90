module dd_testing_setup

  use constants
  use dd_header,        only: dd_type
  use global,           only: master, n_procs, rank, n_particles, n_tallies, &
                              tallies, n_cells, n_user_tallies, micro_xs, &
                              material_xs, total_weight
  use output,           only: write_message
  use particle_header,  only: Particle
  use string,           only: to_str
  use tally,            only: setup_active_usertallies
  use tally,            only: score_analog_tally
  use tally_initialize, only: add_tallies, configure_tallies

#ifdef MPI
  use mpi
#endif

  implicit none
  public

contains

!===============================================================================
! CHECK_PROCS makes sure we have the right number of processors
!===============================================================================

  function check_procs(procs) result(skip)
  
    integer, intent(in) :: procs

    logical :: skip
    skip = .false.
  
#ifdef MPI
    if (.not. n_procs == procs) then
      if (master) call write_message("Skipping test: must be run with " // &
          "MPI and " // trim(to_str(procs)) // " procs")
      skip = .true.
    end if
#else
    call warning("Skipping test: requires MPI")
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

!===============================================================================
! DD_SIMPLE_FOUR_DOMAIN_TALLIES initializes a fake set of tallies designed to
! test on-the-fly memory loading with the simple four-domain test case
! We'll have two scores and a cell filter with 4 cells.  Tests can simulate
! scoring events to this tally where the cells overlap domains in any way.
!===============================================================================

  subroutine dd_simple_four_domain_tallies(p)

    type(Particle), intent(inout) :: p

    ! This hardcodes in what would be done in input_xml and during
    ! initialization, and would need to be modified if anything changes there

    n_tallies = 0 ! This gets set in add_tallies

    n_cells = 4

    n_user_tallies = 1
    call add_tallies("user", n_user_tallies)

    tallies(1) % id = 1
    tallies(1) % label = 'DD test tally for simple four-domain test case'
    tallies(1) %  on_the_fly_allocation = .true.
    tallies(1) %  type = TALLY_VOLUME
    tallies(1) % estimator = ESTIMATOR_ANALOG
    tallies(1) % n_filters = 1
    allocate(tallies(1) % filters(1))
    tallies(1) % filters(1) % type = FILTER_CELL
    tallies(1) % filters(1) % n_bins = 4
    allocate(tallies(1) % filters(1) % int_bins(4))
    tallies(1) % filters(1) % int_bins = (/ 1, 2, 3, 4/)
    allocate(tallies(1) % nuclide_bins(1))
    tallies(1) % nuclide_bins(1) = -1
    tallies(1) % n_nuclide_bins = 1
    allocate(tallies(1) % score_bins(2))
    allocate(tallies(1) % moment_order(2))
    tallies(1) % moment_order = 0
    tallies(1) % score_bins(1) = SCORE_FLUX
    tallies(1) % score_bins(2) = SCORE_FISSION
    tallies(1) % n_score_bins = 2
    tallies(1) % n_user_score_bins = 2

    call configure_tallies()
    call setup_active_usertallies()

    ! Set up a particle to score to this tally
    call p % initialize()

    p % event_nuclide = 1

    allocate(micro_xs(1))
    micro_xs % absorption = 1.0_8
    micro_xs % fission = 1.0_8

    material_xs % total = 1.0_8

  end subroutine dd_simple_four_domain_tallies

!===============================================================================
! DD_SCORE_TO_FOUR_DOMAIN_TALLIES
!===============================================================================

  subroutine dd_score_to_four_domain_tallies(p)

    type(Particle), intent(inout) :: p

    ! Simulate scoring to tally bins
    select case(rank)
      case(0)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        total_weight = 7.0
      case(1)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 1
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        total_weight = 7.0
      case(2)
        p % coord % cell = 3
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        p % coord % cell = 3
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        total_weight = 7.0
      case(3)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        total_weight = 2.0
      case(4)
        p % coord % cell = 4
        call score_analog_tally(p)
        p % coord % cell = 3
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
        total_weight = 5.0
    end select

  end subroutine dd_score_to_four_domain_tallies

end module dd_testing_setup
