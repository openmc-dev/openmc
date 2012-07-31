module criticality

  use constants,   only: ZERO
  use global
  use intercycle,  only: shannon_entropy, calculate_keff, synchronize_bank, &
                         count_source_for_ufs
  use output,      only: write_message, header, print_columns
  use physics,     only: transport
  use source,      only: get_source_particle, write_source_binary
  use state_point, only: create_state_point
  use string,      only: to_str
  use tally,       only: synchronize_tallies
  use timing,      only: timer_start, timer_stop

contains

!===============================================================================
! RUN_CRITICALITY encompasses all the main logic where iterations are performed
! over the cycles and histories.
!===============================================================================

  subroutine run_criticality()

    integer(8) :: i  ! index over histories in single cycle

    if (master) call header("CRITICALITY TRANSPORT SIMULATION", level=1)

    tallies_on = .false.
    call timer_start(time_inactive)

    ! Allocate particle
    allocate(p)

    ! Display column titles
    call print_columns()

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_batches

       ! In a restart run, skip any batches that have already been simulated
       if (restart_run .and. current_batch <= restart_batch) then
          if (current_batch == n_inactive) tallies_on = .true.
          if (current_batch > n_inactive) n_realizations = n_realizations + 1
          cycle BATCH_LOOP
       end if

       call initialize_batch()

       ! =======================================================================
       ! LOOP OVER GENERATIONS
       GENERATION_LOOP: do current_gen = 1, gen_per_batch

          call initialize_generation()

          ! Start timer for transport
          call timer_start(time_transport)

          ! ====================================================================
          ! LOOP OVER PARTICLES
          PARTICLE_LOOP: do i = 1, work

             ! grab source particle from bank
             call get_source_particle(i)

             ! transport particle
             call transport()

          end do PARTICLE_LOOP

          ! Accumulate time for transport
          call timer_stop(time_transport)

          ! Distribute fission bank across processors evenly
          call timer_start(time_intercycle)
          call synchronize_bank()
          call timer_stop(time_intercycle)
          
       end do GENERATION_LOOP

       call finalize_batch()

    end do BATCH_LOOP

    call timer_stop(time_active)

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_criticality

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

       message = "Simulating batch " // trim(to_str(current_batch)) // "..."
       call write_message(8)

       ! Reset total starting particle weight used for normalizing tallies
       total_weight = ZERO

  end subroutine initialize_batch

!===============================================================================
! INITIALIZE_GENERATION
!===============================================================================

  subroutine initialize_generation()

    ! Reset number of fission bank sites
    n_bank = 0

    ! Count source sites if using uniform fission source weighting
    if (ufs) call count_source_for_ufs()

  end subroutine initialize_generation

!===============================================================================
! FINALIZE_BATCH handles synchronization and accumulation of tallies,
! calculation of Shannon entropy, getting single-batch estimate of keff, and
! turning on tallies when appropriate
!===============================================================================

  subroutine finalize_batch()

    integer :: i ! loop index for state point batches

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

    ! Write out state point if it's been specified for this batch
    do i = 1, n_state_points
       if (current_batch == statepoint_batch(i)) then
          ! Create state point file
          if (master) call create_state_point()

          ! Create binary source file
          call write_source_binary('source.' // &
               trim(to_str(current_batch)) // '.binary')
          exit
       end if
    end do

    ! Turn tallies on once inactive cycles are complete
    if (current_batch == n_inactive) then
       tallies_on = .true.
       call timer_stop(time_inactive)
       call timer_start(time_active)
    end if

  end subroutine finalize_batch

end module criticality
