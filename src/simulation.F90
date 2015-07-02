module simulation

#ifdef MPI
  use mpi
#endif

  use cmfd_execute, only: cmfd_init_batch, execute_cmfd
  use constants,    only: ZERO
  use eigenvalue,   only: count_source_for_ufs, calculate_average_keff, &
                          calculate_combined_keff, calculate_generation_keff, &
                          shannon_entropy, synchronize_bank, keff_generation
#ifdef _OPENMP
  use eigenvalue,   only: join_bank_from_threads
#endif
  use global
  use output,       only: write_message, header, print_columns, &
                          print_batch_keff, print_generation
  use particle_header, only: Particle
  use source,       only: get_source_particle, initialize_source
  use state_point,  only: write_state_point, write_source_point
  use string,       only: to_str
  use tally,        only: synchronize_tallies, setup_active_usertallies, &
                          reset_result
  use trigger,      only: check_triggers
  use tracking,     only: transport

  implicit none
  private
  public :: run_simulation

contains

!===============================================================================
! RUN_EIGENVALUE encompasses all the main logic where iterations are performed
! over the batches, generations, and histories in a k-eigenvalue calculation.
!===============================================================================

  subroutine run_simulation()

    type(Particle) :: p
    integer(8)     :: i_work

    if (.not. restart_run) call initialize_source()

    ! Display header
    if (master) then
      if (run_mode == MODE_FIXEDSOURCE) then
        call header("FIXED SOURCE TRANSPORT SIMULATION", level=1)
      elseif (run_mode == MODE_EIGENVALUE) then
        call header("K EIGENVALUE SIMULATION", level=1)
        call print_columns()
      end if
    end if

    ! Turn on inactive timer
    call time_inactive % start()

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_max_batches

      call initialize_batch()

      ! Handle restart runs
      if (restart_run .and. current_batch <= restart_batch) then
        call replay_batch_history()
        cycle BATCH_LOOP
      end if

      ! =======================================================================
      ! LOOP OVER GENERATIONS
      GENERATION_LOOP: do current_gen = 1, gen_per_batch

        call initialize_generation()

        ! Start timer for transport
        call time_transport % start()

        ! ====================================================================
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

        call finalize_generation()

      end do GENERATION_LOOP

      call finalize_batch()

      if (satisfy_triggers) exit BATCH_LOOP

    end do BATCH_LOOP

    call time_active % stop()

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

    ! Clear particle
    call p % clear()

  end subroutine run_simulation

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    if (run_mode == MODE_FIXEDSOURCE) then
      call write_message("Simulating batch " // trim(to_str(current_batch)) &
           // "...", 1)
    end if

    ! Reset total starting particle weight used for normalizing tallies
    total_weight = ZERO

    if (current_batch == n_inactive + 1) then
      ! Switch from inactive batch timer to active batch timer
      call time_inactive % stop()
      call time_active % start()

      ! Enable active batches (and tallies_on if it hasn't been enabled)
      active_batches = .true.
      tallies_on = .true.

      ! Add user tallies to active tallies list
!$omp parallel
      call setup_active_usertallies()
!$omp end parallel
    end if

    ! check CMFD initialize batch
    if (run_mode == MODE_EIGENVALUE) then
      if (cmfd_run) call cmfd_init_batch()
    end if

  end subroutine initialize_batch

!===============================================================================
! INITIALIZE_GENERATION
!===============================================================================

  subroutine initialize_generation()

    ! set overall generation number
    overall_gen = gen_per_batch*(current_batch - 1) + current_gen

    if (run_mode == MODE_EIGENVALUE) then
      ! Reset number of fission bank sites
      n_bank = 0

      ! Count source sites if using uniform fission source weighting
      if (ufs) call count_source_for_ufs()

      ! Store current value of tracklength k
      keff_generation = global_tallies(K_TRACKLENGTH) % value
    end if

  end subroutine initialize_generation

!===============================================================================
! FINALIZE_GENERATION
!===============================================================================

  subroutine finalize_generation()

    ! Update global tallies with the omp private accumulation variables
!$omp parallel
!$omp critical
    if (run_mode == MODE_EIGENVALUE) then
      global_tallies(K_COLLISION) % value = &
           global_tallies(K_COLLISION) % value + global_tally_collision
      global_tallies(K_ABSORPTION) % value = &
           global_tallies(K_ABSORPTION) % value + global_tally_absorption
      global_tallies(K_TRACKLENGTH) % value = &
           global_tallies(K_TRACKLENGTH) % value + global_tally_tracklength
    end if
    global_tallies(LEAKAGE) % value = &
         global_tallies(LEAKAGE) % value + global_tally_leakage
!$omp end critical

    ! reset private tallies
    if (run_mode == MODE_EIGENVALUE) then
      global_tally_collision = 0
      global_tally_absorption = 0
      global_tally_tracklength = 0
    end if
    global_tally_leakage = 0
!$omp end parallel

    if (run_mode == MODE_EIGENVALUE) then
#ifdef _OPENMP
      ! Join the fission bank from each thread into one global fission bank
      call join_bank_from_threads()
#endif

      ! Distribute fission bank across processors evenly
      call time_bank % start()
      call synchronize_bank()
      call time_bank % stop()

      ! Calculate shannon entropy
      if (entropy_on) call shannon_entropy()

      ! Collect results and statistics
      call calculate_generation_keff()
      call calculate_average_keff()

      ! Write generation output
      if (master .and. current_gen /= gen_per_batch) call print_generation()
    end if

  end subroutine finalize_generation

!===============================================================================
! FINALIZE_BATCH handles synchronization and accumulation of tallies,
! calculation of Shannon entropy, getting single-batch estimate of keff, and
! turning on tallies when appropriate
!===============================================================================

  subroutine finalize_batch()

    ! Collect tallies
    call time_tallies % start()
    call synchronize_tallies()
    call time_tallies % stop()

    ! Reset global tally results
    if (.not. active_batches) then
      call reset_result(global_tallies)
      n_realizations = 0
    end if

    if (run_mode == MODE_EIGENVALUE) then
      ! Perform CMFD calculation if on
      if (cmfd_on) call execute_cmfd()

      ! Display output
      if (master) call print_batch_keff()

      ! Calculate combined estimate of k-effective
      if (master) call calculate_combined_keff()
    end if

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

    ! Write out source point if it's been specified for this batch
    if ((sourcepoint_batch % contains(current_batch) .or. source_latest) .and. &
         source_write) then
      call write_source_point()
    end if

    if (master .and. current_batch == n_max_batches .and. &
         run_mode == MODE_EIGENVALUE) then
      ! Make sure combined estimate of k-effective is calculated at the last
      ! batch in case no state point is written
      call calculate_combined_keff()
    end if

  end subroutine finalize_batch

!===============================================================================
! REPLAY_BATCH_HISTORY displays keff and entropy for each generation within a
! batch using data read from a state point file
!===============================================================================

  subroutine replay_batch_history

    ! Write message at beginning
    if (current_batch == 1) then
      call write_message("Replaying history from state point...", 1)
    end if

    if (run_mode == MODE_EIGENVALUE) then
      do current_gen = 1, gen_per_batch
        overall_gen = overall_gen + 1
        call calculate_average_keff()

        ! print out batch keff
        if (current_gen < gen_per_batch) then
          if (master) call print_generation()
        else
          if (master) call print_batch_keff()
        end if
      end do
    end if

    ! Write message at end
    if (current_batch == restart_batch) then
      call write_message("Resuming simulation...", 1)
    end if

  end subroutine replay_batch_history

end module simulation
