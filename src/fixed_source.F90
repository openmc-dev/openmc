module fixed_source

#ifdef MPI
  use message_passing
#endif

  use constants,       only: ZERO, MAX_LINE_LEN
  use global
  use output,          only: write_message, header
  use particle_header, only: Particle
  use random_lcg,      only: set_particle_seed
  use source,          only: initialize_source, get_source_particle
  use state_point,     only: write_state_point
  use string,          only: to_str
  use tally,           only: synchronize_tallies, setup_active_usertallies
  use trigger,         only: check_triggers
  use tracking,        only: transport

  implicit none

contains

  subroutine run_fixedsource()

    type(Particle) :: p
    integer(8)     :: i_work ! index over histories in single cycle

    if (.not. restart_run) call initialize_source()

    if (master) call header("FIXED SOURCE TRANSPORT SIMULATION", level=1)

    ! Turn timer and tallies on
    tallies_on = .true.
!$omp parallel
    call setup_active_usertallies()
!$omp end parallel
    call time_active % start()

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_max_batches

      ! In a restart run, skip any batches that have already been simulated
      if (restart_run .and. current_batch <= restart_batch) then
        if (current_batch > n_inactive) n_realizations = n_realizations + 1
        cycle BATCH_LOOP
      end if

      call initialize_batch()
      overall_gen = current_batch

      ! Start timer for transport
      call time_transport % start()

      ! =======================================================================
      ! LOOP OVER PARTICLES
!$omp parallel do schedule(static) firstprivate(p)
      PARTICLE_LOOP: do i_work = 1, work
        current_work = i_work

        ! grab source particle from bank
        call get_source_particle(p, current_work)

        ! transport particle
        call transport(p)

      end do PARTICLE_LOOP
!$omp end parallel do

      ! Accumulate time for transport
      call time_transport % stop()

      call finalize_batch()

      if (satisfy_triggers) exit BATCH_LOOP

    end do BATCH_LOOP

    call time_active % stop()

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_fixedsource

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    call write_message("Simulating batch " // trim(to_str(current_batch)) &
         &// "...", 1)

    ! Reset total starting particle weight used for normalizing tallies
    total_weight = ZERO

  end subroutine initialize_batch

!===============================================================================
! FINALIZE_BATCH
!===============================================================================

  subroutine finalize_batch()

! Update global tallies with the omp private accumulation variables
!$omp parallel
!$omp critical
    global_tallies(LEAKAGE) % value = &
         global_tallies(LEAKAGE) % value + global_tally_leakage
!$omp end critical

    ! reset private tallies
    global_tally_leakage = ZERO
!$omp end parallel

    ! Collect and accumulate tallies
    call time_tallies % start()
    call synchronize_tallies()
    call time_tallies % stop()

    ! Check_triggers
    if (master) call check_triggers()
#ifdef MPI
    call MPI_BCAST(satisfy_triggers, 1, MPI_LOGICAL, 0, &
         MPI_COMM_WORLD, mpi_err)
#endif
    if (satisfy_triggers .or. &
         (trigger_on .and. current_batch == n_max_batches)) then
      call statepoint_batch % add(current_batch)
    end if

    ! Write out state point if it's been specified for this batch
    if (statepoint_batch % contains(current_batch)) then
      call write_state_point()
    end if

  end subroutine finalize_batch

end module fixed_source
