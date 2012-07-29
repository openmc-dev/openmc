module fixed_source

  use constants,  only: ZERO
  use global
  use output,     only: write_message, header
  use physics,    only: transport
  use source,     only: get_source_particle
  use string,     only: to_str
  use tally,      only: synchronize_tallies
  use timing,     only: timer_start, timer_stop

contains

  subroutine run_fixedsource()

    integer(8) :: i  ! index over histories in single cycle

    if (master) call header("BEGIN SIMULATION", level=1)

    tallies_on = .true.
    call timer_start(time_inactive)

    ! Allocate particle
    allocate(p)

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_batches

       call initialize_batch()

       ! Start timer for transport
       call timer_start(time_transport)

       ! =======================================================================
       ! LOOP OVER PARTICLES
       PARTICLE_LOOP: do i = 1, work

          ! grab source particle from bank
          call get_source_particle(i)

          ! transport particle
          call transport()

       end do PARTICLE_LOOP

       ! Accumulate time for transport
       call timer_stop(time_transport)

       call timer_start(time_ic_tallies)
       call synchronize_tallies()
       call timer_stop(time_ic_tallies)

    end do BATCH_LOOP

    call timer_stop(time_active)

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_fixedsource

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

       message = "Simulating batch " // trim(to_str(current_batch)) // "..."
       call write_message()

       ! Reset total starting particle weight used for normalizing tallies
       total_weight = ZERO

  end subroutine initialize_batch

end module fixed_source
