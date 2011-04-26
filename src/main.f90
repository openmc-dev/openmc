program main

  use global
  use fileio,        only: read_input, read_command_line, read_count, &
       &                   normalize_ao, build_universe
  use output,        only: title, echo_input, message, warning, error, &
       &                   print_summary, print_particle
  use geometry,      only: neighbor_lists
  use mcnp_random,   only: RN_init_problem, RN_init_particle
  use source,        only: init_source, get_source_particle
  use physics,       only: transport
  use cross_section, only: read_xsdata, material_total_xs
  use ace,           only: read_xs
  use energy_grid,   only: unionized_grid, original_indices
  use mpi_routines,  only: setup_mpi, synchronize_bank, t_sync
  use score,         only: calculate_keff

#ifdef MPI
  use mpi
#endif

  implicit none

  character(16) :: filename
  character(250) :: msg
  type(Universe), pointer :: univ

  ! Setup MPI
  call setup_mpi()

  ! Print the OpenMC title and version/date/time information
  if (master) call title()

  ! Initialize random number generator. The first argument corresponds
  ! to which random number generator to use- in this case one of the
  ! L'Ecuyer 63-bit RNGs.
  call RN_init_problem(3, 0_8, 0_8, 0_8, 0)

  ! Set default values for settings
  call set_defaults()

  ! Read command line arguments
  call read_command_line()

  ! Read input file -- make a first pass through the file to count
  ! cells, surfaces, etc in order to allocate arrays, then do a second
  ! pass to actually read values
  call read_count(path_input)
  call read_input(path_input)

  ! determine at which level universes are and link cells to parenting
  ! cells
  univ => universes(BASE_UNIVERSE)
  call build_universe(univ, 0, 0)

  ! After reading input and basic geometry setup is complete, build
  ! lists of neighboring cells for efficient tracking
  call neighbor_lists()

  ! Read cross section summary file to determine what files contain
  ! cross-sections
  call read_xsdata(path_xsdata)

  ! With the AWRs from the xsdata, change all material specifications
  ! so that they contain atom percents summing to 1
  call normalize_ao()

  ! Read ACE-format cross sections
  call read_xs()

  ! Construct unionized energy grid from cross-sections
  call unionized_grid()
  call original_indices()

  ! calculate total material cross-sections for sampling path lenghts
  call material_total_xs()

  if (master) then
     call echo_input()
     call print_summary()
  end if

  ! create source particles
  call init_source()

  ! start problem
  call run_problem()

  ! deallocate arrays
  call free_memory()
  
contains

!=====================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are
! performed over the cycles and histories.
!=====================================================================

  subroutine run_problem()

    integer :: i, j
    integer :: ierr
    integer :: i_cycle    ! cycle index
    integer :: i_particle ! history index
    integer :: total_bank ! total number of particles banked
    real(8) :: t0
    real(8) :: t1
    type(Particle), pointer :: p => null()
    character(250) :: msg

    msg = "Running problem..."
    call message(msg, 6)

    tallies_on = .false.

#ifdef MPI
    t0 = MPI_WTIME()
#endif

    CYCLE_LOOP: do i_cycle = 1, n_cycles

       msg = "Simulating cycle " // trim(int_to_str(i_cycle)) // "..."
       call message(msg, 8)
       
       ! Set all tallies to zero
       n_bank = 0

       HISTORY_LOOP: do

          ! grab source particle from bank
          p => get_source_particle()
          if ( .not. associated(p) ) then
             ! no particles left in source bank
             exit HISTORY_LOOP
          end if

          ! set random number seed
          i_particle = (i_cycle-1)*n_particles + p % uid
          call RN_init_particle(int(i_particle,8))

          ! transport particle
          call transport(p)

       end do HISTORY_LOOP

       call synchronize_bank(i_cycle)

       ! Collect results and statistics
       call calculate_keff(i_cycle)

       ! print cycle information

       ! Turn tallies on once inactive cycles are complete
       if (i_cycle == n_inactive) tallies_on = .true.

    end do CYCLE_LOOP

#ifdef MPI
    ! print run time
    t1 = MPI_WTIME()
    if (master) then
       print *, "Time elapsed   = " // real_to_str(t1 - t0)
       print *, "Init time      = " // real_to_str(t_sync(1))
       print *, "Sample time    = " // real_to_str(t_sync(2))
       print *, "Send/recv time = " // real_to_str(t_sync(3))
       print *, "Rebuild time   = " // real_to_str(t_sync(4))
    end if
#endif

  end subroutine run_problem

end program main


