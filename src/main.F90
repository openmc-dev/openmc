program main

  use constants
  use global
  use finalize,   only: finalize_run
  use initialize, only: initialize_run
  use intercycle, only: shannon_entropy, calculate_keff, synchronize_bank, &
                        count_source_for_ufs
  use output,     only: write_message, header
  use plotter,    only: run_plot
  use physics,    only: transport
  use source,     only: get_source_particle
  use string,     only: to_str
  use tally,      only: synchronize_tallies
  use timing,     only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

  implicit none

  ! set up problem
  call initialize_run()

  ! start problem
  if (plotting) then
     call run_plot()
  else
     call run_problem()
  end if
     
  ! finalize run
  call finalize_run()
  
contains

!===============================================================================
! RUN_PROBLEM encompasses all the main logic where iterations are performed over
! the cycles and histories.
!===============================================================================

  subroutine run_problem()

    integer(8) :: i  ! index over histories in single cycle

    if (master) call header("BEGIN SIMULATION", level=1)

    tallies_on = .false.
    call timer_start(time_inactive)

    ! Allocate particle
    allocate(p)

    ! Display column titles
    if (entropy_on) then
       message = " Batch   k(batch)   Entropy         Average k"
       call write_message(1)
       message = " =====   ========   =======    ==================="
       call write_message(1)
    else
       message = " Batch   k(batch)          Average k"
       call write_message(1)
       message = " =====   ========     ==================="
       call write_message(1)
    end if

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_batches

       message = "Simulating batch " // trim(to_str(current_batch)) // "..."
       call write_message(8)

       ! =======================================================================
       ! LOOP OVER GENERATIONS
       GENERATION_LOOP: do current_gen = 1, gen_per_batch

          ! Set all tallies to zero
          n_bank = 0

          ! Count source sites if using uniform fission source weighting
          if (ufs) call count_source_for_ufs()

          ! ====================================================================
          ! LOOP OVER HISTORIES

          ! Start timer for transport
          call timer_start(time_transport)

          HISTORY_LOOP: do i = 1, work

             ! grab source particle from bank
             call get_source_particle(i)

             ! transport particle
             call transport()

          end do HISTORY_LOOP

          ! Accumulate time for transport
          call timer_stop(time_transport)

          ! ====================================================================
          ! WRAP UP FISSION BANK AND COMPUTE TALLIES, KEFF, ETC

          ! Start timer for inter-cycle synchronization
          call timer_start(time_intercycle)

          ! Distribute fission bank across processors evenly
          call synchronize_bank()

          ! Stop timer for inter-cycle synchronization
          call timer_stop(time_intercycle)
          
       end do GENERATION_LOOP

       ! Collect tallies
       if (tallies_on) then
          call timer_start(time_ic_tallies)
          call synchronize_tallies()
          call timer_stop(time_ic_tallies)
       end if

       ! Calculate shannon entropy
       if (entropy_on) call shannon_entropy()

       ! Collect results and statistics
       call calculate_keff()

       ! Turn tallies on once inactive cycles are complete
       if (current_batch == n_inactive) then
          tallies_on = .true.
          call timer_stop(time_inactive)
          call timer_start(time_active)
       end if

    end do BATCH_LOOP

    call timer_stop(time_active)

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_problem

end program main
