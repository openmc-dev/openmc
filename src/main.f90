program main

  use global
  use fileio,      only: read_input, read_command_line, read_count
  use output,      only: title, echo_input, message, warning, error
  use geometry,    only: sense, cell_contains
  use mcnp_random, only: RN_init_problem, rang, RN_init_particle
  use source,      only: init_source

  implicit none

  character(16) :: filename
  character(250) :: msg

  real(8) :: point(3)
  integer :: s(4)

  ! Print the OpenMC title and version/date/time information
  call title()

  ! Initialize random number generator
  call RN_init_problem( 3, 0_8, 0_8, 0_8, 0 )

  ! Read command line arguments
  call read_command_line()

  ! Read input file -- make a first pass through the file to count
  ! cells, surfaces, etc in order to allocate arrays, then do a second
  ! pass to actually read values
  call read_count(inputfile)
  call read_input(inputfile)
  call echo_input()

  ! create source particles
  call init_source()

  ! start problem
  ! call run_problem()


contains

!=====================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are
! performed over the cycles and histories.
!=====================================================================

  subroutine run_problem()

    integer :: i, j
    integer :: i_cycle
    integer :: i_particle
    type(Neutron), pointer :: particle => null()

    CYCLE_LOOP: do i_cycle = 1, n_cycles
       
       ! Set all tallies to zero

       HISTORY_LOOP: do j = 1, n_particles

          ! grab source particle from bank
          ! particle => get_source_particle()

          ! set random number seed
          i_particle = (i-1)*n_particles + j
          call RN_init_particle(int(i_particle,8))

          ! transport particle
          ! call transport(particle)

       end do HISTORY_LOOP

       ! Collect results and statistics

       ! print cycle information

    end do CYCLE_LOOP

    ! Collect all tallies and print

    ! print run time

    call free_memory()
    
  end subroutine run_problem

end program main


