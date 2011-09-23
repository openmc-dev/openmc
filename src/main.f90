program main

  use constants
  use global
  use initialize,      only: initialize_run
  use mcnp_random,     only: RN_init_particle
  use mpi_routines,    only: synchronize_bank
  use output,          only: message, header, print_runtime
  use particle_header, only: Particle
  use physics,         only: transport
  use tally,           only: calculate_keff
  use source,          only: get_source_particle
  use string,          only: int_to_str
  use tally,           only: synchronize_tallies
  use timing,          only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

  implicit none

  ! start timer for total run time
  call timer_start(time_total)

  ! set up problem
  call initialize_run()

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

    integer                 :: i_cycle     ! cycle index
    integer(8)              :: i_particle  ! history index
    character(MAX_LINE_LEN) :: msg         ! output/error message
    type(Particle), pointer :: p => null()

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

       ! Collect tallies
       call synchronize_tallies()

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


