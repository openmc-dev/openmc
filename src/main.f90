program main

  use global
  use fileio,        only: read_input, read_command_line, read_count, &
       &                   normalize_ao, build_universe
  use output,        only: title, echo_input, message, warning, error, &
       &                   print_summary
  use geometry,      only: sense, cell_contains, neighbor_lists
  use mcnp_random,   only: RN_init_problem, rang, RN_init_particle
  use source,        only: init_source, get_source_particle
  use physics,       only: transport
  use cross_section, only: read_xsdata, material_total_xs
  use ace,           only: read_xs
  use energy_grid,   only: unionized_grid, original_indices

  implicit none

  character(16) :: filename
  character(250) :: msg
  type(Universe), pointer :: univ

  ! Print the OpenMC title and version/date/time information
  call title()
  verbosity = 10

  ! Initialize random number generator
  call RN_init_problem( 3, 0_8, 0_8, 0_8, 0 )

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

  call echo_input()
  call print_summary()

  ! create source particles
  call init_source()

  ! start problem
  surfaces(3)%bc = BC_VACUUM
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
    integer :: i_cycle
    integer :: i_particle
    type(Particle), pointer :: p => null()

    CYCLE_LOOP: do i_cycle = 1, n_cycles
       
       ! Set all tallies to zero

       HISTORY_LOOP: do j = 1, n_particles

          ! grab source particle from bank
          p => get_source_particle()
          if ( .not. associated(p) ) then
             ! no particles left in source bank
             exit HISTORY_LOOP
          end if

          ! set random number seed
          i_particle = (i-1)*n_particles + j
          call RN_init_particle(int(i_particle,8))

          ! transport particle
          call transport(p)

       end do HISTORY_LOOP

       ! Collect results and statistics

       ! print cycle information

    end do CYCLE_LOOP

    ! Collect all tallies and print

    ! print run time

  end subroutine run_problem

end program main


