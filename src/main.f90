program main

  use constants
  use cross_section,   only: read_xs, read_xsdata, material_total_xs
  use energy_grid,     only: unionized_grid, original_indices
  use fileio,          only: read_input, read_command_line, read_count, &
                             normalize_ao, build_universe
  use geometry,        only: neighbor_lists
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use logging,         only: create_log
  use mcnp_random,     only: RN_init_problem, RN_init_particle
  use mpi_routines,    only: setup_mpi, synchronize_bank
  use output,          only: title, echo_input, message, print_summary, &
                             print_particle, header, print_runtime
  use particle_header, only: Particle
  use physics,         only: transport
  use score,           only: calculate_keff
  use source,          only: init_source, get_source_particle
  use string,          only: int_to_str
  use timing,          only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

  implicit none

  character(max_word_len) :: filename
  character(max_line_len) :: msg
  type(Universe), pointer :: univ

  ! Start timers
  call timer_start(time_total)
  call timer_start(time_init)

  ! Setup MPI
  call setup_mpi()

  ! Read command line arguments
  call read_command_line()
  if (master) call create_log()

  ! Print the OpenMC title and version/date/time information
  if (master) call title()
  ! Initialize random number generator. The first argument corresponds to which
  ! random number generator to use- in this case one of the L'Ecuyer 63-bit
  ! RNGs.

  ! Print initialization header block
  if (master) call header("INITIALIZATION", 1)

  ! initialize random number generator
  call RN_init_problem(3, 0_8, 0_8, 0_8, 0)

  ! Set default values for settings
  call set_defaults()

  ! Read input file -- make a first pass through the file to count cells,
  ! surfaces, etc in order to allocate arrays, then do a second pass to actually
  ! read values
  call read_count(path_input)
  call read_input(path_input)

  ! determine at which level universes are and link cells to parenting cells
  univ => universes(BASE_UNIVERSE)
  call build_universe(univ, 0, 0)

  ! After reading input and basic geometry setup is complete, build lists of
  ! neighboring cells for efficient tracking
  call neighbor_lists()

  ! Read cross section summary file to determine what files contain
  ! cross-sections
  call read_xsdata(path_xsdata)

  ! With the AWRs from the xsdata, change all material specifications so that
  ! they contain atom percents summing to 1
  call normalize_ao()

  ! Read ACE-format cross sections
  call read_xs()

  ! Construct unionized energy grid from cross-sections
  call unionized_grid()
  call original_indices()

  ! calculate total material cross-sections for sampling path lenghts
  call material_total_xs()

  ! create source particles
  call init_source()

  ! stop timer for initialization
  call timer_stop(time_init)
  if (master) then
     call echo_input()
     call print_summary()
  end if

  ! start problem
  call run_problem()

  ! show timing statistics
  call timer_stop(time_total)
  if (master) call print_runtime()

  ! deallocate arrays
  call free_memory()
  
contains

!===============================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are performed over
! the cycles and histories.
!===============================================================================

  subroutine run_problem()

    integer    :: i, j
    integer    :: ierr
    integer    :: i_cycle    ! cycle index
    integer(8) :: i_particle ! history index
    integer    :: total_bank ! total number of particles banked
    real(8)    :: t0
    real(8)    :: t1
    type(Particle), pointer :: p => null()
    character(max_line_len) :: msg

    if (master) call header("BEGIN SIMULATION", 1)

    tallies_on = .false.

    ! ==========================================================================
    ! LOOP OVER CYCLES
    CYCLE_LOOP: do i_cycle = 1, n_cycles

       ! Start timer for computation
       call timer_start(time_compute)

       msg = "Simulating cycle " // trim(int_to_str(i_cycle)) // "..."
       call message(msg, 8)
       
       ! Set all tallies to zero
       n_bank = 0

       ! =======================================================================
       ! LOOP OVER HISTORIES
       HISTORY_LOOP: do

          ! grab source particle from bank
          p => get_source_particle()
          if ( .not. associated(p) ) then
             ! no particles left in source bank
             exit HISTORY_LOOP
          end if

          ! set random number seed
          i_particle = (i_cycle-1)*n_particles + p % uid
          call RN_init_particle(i_particle)

          ! transport particle
          call transport(p)

       end do HISTORY_LOOP

       ! Accumulate time for computation
       call timer_stop(time_compute)

       ! =======================================================================
       ! WRAP UP FISSION BANK AND COMPUTE TALLIES, KEFF, ETC

       ! Start timer for inter-cycle synchronization
       call timer_start(time_intercycle)

       ! Distribute fission bank across processors evenly
       call synchronize_bank(i_cycle)

       ! Collect results and statistics
       call calculate_keff(i_cycle)

       ! print cycle information

       ! Turn tallies on once inactive cycles are complete
       if (i_cycle == n_inactive) tallies_on = .true.

       ! Stop timer for inter-cycle synchronization
       call timer_stop(time_intercycle)

    end do CYCLE_LOOP

  end subroutine run_problem

end program main


