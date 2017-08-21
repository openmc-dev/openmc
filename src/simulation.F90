module simulation

  use, intrinsic :: ISO_C_BINDING

#ifdef _OPENMP
  use omp_lib
#endif

  use cmfd_execute,    only: cmfd_init_batch, execute_cmfd
  use constants,       only: ZERO
  use eigenvalue,      only: count_source_for_ufs, calculate_average_keff, &
                             calculate_generation_keff, shannon_entropy, &
                             synchronize_bank, keff_generation, k_sum
#ifdef _OPENMP
  use eigenvalue,      only: join_bank_from_threads
#endif
  use global
  use message_passing
  use output,          only: write_message, header, print_columns, &
                             print_batch_keff, print_generation, print_runtime, &
                             print_results, print_overlap_check, write_tallies
  use particle_header, only: Particle
  use random_lcg,      only: set_particle_seed
  use source,          only: initialize_source, sample_external_source
  use state_point,     only: write_state_point, write_source_point, load_state_point
  use string,          only: to_str
  use tally,           only: accumulate_tallies, setup_active_tallies, &
                             init_tally_routines
  use tally_header,    only: configure_tallies
  use trigger,         only: check_triggers
  use tracking,        only: transport

  implicit none
  private
  public :: openmc_run

contains

!===============================================================================
! OPENMC_RUN encompasses all the main logic where iterations are performed
! over the batches, generations, and histories in a fixed source or k-eigenvalue
! calculation.
!===============================================================================

  subroutine openmc_run() bind(C)

    type(Particle) :: p
    integer(8)     :: i_work

    call initialize_simulation()

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
!$omp parallel do schedule(static) firstprivate(p) copyin(tally_derivs)
        PARTICLE_LOOP: do i_work = 1, work
          current_work = i_work

          ! grab source particle from bank
          call initialize_history(p, current_work)

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

    call finalize_simulation()

    ! Clear particle
    call p % clear()

  end subroutine openmc_run

!===============================================================================
! INITIALIZE_HISTORY
!===============================================================================

  subroutine initialize_history(p, index_source)

    type(Particle), intent(inout) :: p
    integer(8),     intent(in)    :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    integer :: i

    ! set defaults
    call p % initialize_from_source(source_bank(index_source), run_CE, &
         energy_bin_avg)

    ! set identifier for particle
    p % id = work_index(rank) + index_source

    ! set random number seed
    particle_seed = (total_gen + overall_generation() - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! set particle trace
    trace = .false.
    if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
         p % id == trace_particle) trace = .true.

    ! Set particle track.
    p % write_track = .false.
    if (write_all_tracks) then
      p % write_track = .true.
    else if (allocated(track_identifiers)) then
      do i=1, size(track_identifiers(1,:))
        if (current_batch == track_identifiers(1,i) .and. &
             &current_gen == track_identifiers(2,i) .and. &
             &p % id == track_identifiers(3,i)) then
          p % write_track = .true.
          exit
        end if
      end do
    end if

  end subroutine initialize_history

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    integer :: i

    if (run_mode == MODE_FIXEDSOURCE) then
      call write_message("Simulating batch " // trim(to_str(current_batch)) &
           // "...", 6)
    end if

    ! Reset total starting particle weight used for normalizing tallies
    total_weight = ZERO

    if (current_batch == n_inactive + 1) then
      ! Switch from inactive batch timer to active batch timer
      call time_inactive % stop()
      call time_active % start()

      do i = 1, n_tallies
        tallies(i) % obj % active = .true.
      end do
    end if

    ! check CMFD initialize batch
    if (run_mode == MODE_EIGENVALUE) then
      if (cmfd_run) call cmfd_init_batch()
    end if

    ! Add user tallies to active tallies list
    call setup_active_tallies()

  end subroutine initialize_batch

!===============================================================================
! INITIALIZE_GENERATION
!===============================================================================

  subroutine initialize_generation()

    if (run_mode == MODE_EIGENVALUE) then
      ! Reset number of fission bank sites
      n_bank = 0

      ! Count source sites if using uniform fission source weighting
      if (ufs) call count_source_for_ufs()

      ! Store current value of tracklength k
      keff_generation = global_tallies(RESULT_VALUE, K_TRACKLENGTH)
    end if

  end subroutine initialize_generation

!===============================================================================
! FINALIZE_GENERATION
!===============================================================================

  subroutine finalize_generation()

    integer(8) :: i

    ! Update global tallies with the omp private accumulation variables
!$omp parallel
!$omp critical
    if (run_mode == MODE_EIGENVALUE) then
      global_tallies(RESULT_VALUE, K_COLLISION) = &
           global_tallies(RESULT_VALUE, K_COLLISION) + global_tally_collision
      global_tallies(RESULT_VALUE, K_ABSORPTION) = &
           global_tallies(RESULT_VALUE, K_ABSORPTION) + global_tally_absorption
      global_tallies(RESULT_VALUE, K_TRACKLENGTH) = &
           global_tallies(RESULT_VALUE, K_TRACKLENGTH) + global_tally_tracklength
    end if
    global_tallies(RESULT_VALUE, LEAKAGE) = &
         global_tallies(RESULT_VALUE, LEAKAGE) + global_tally_leakage
!$omp end critical

    ! reset private tallies
    if (run_mode == MODE_EIGENVALUE) then
      global_tally_collision = ZERO
      global_tally_absorption = ZERO
      global_tally_tracklength = ZERO
    end if
    global_tally_leakage = ZERO
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
      if (master .and. verbosity >= 7) then
        if (current_gen /= gen_per_batch) then
          call print_generation()
        else
          call print_batch_keff()
        end if
      end if

    elseif (run_mode == MODE_FIXEDSOURCE) then
      ! For fixed-source mode, we need to sample the external source
      if (path_source == '') then
        do i = 1, work
          call set_particle_seed((total_gen + overall_generation()) * &
               n_particles + work_index(rank) + i)
          call sample_external_source(source_bank(i))
        end do
      end if
    end if

  end subroutine finalize_generation

!===============================================================================
! FINALIZE_BATCH handles synchronization and accumulation of tallies,
! calculation of Shannon entropy, getting single-batch estimate of keff, and
! turning on tallies when appropriate
!===============================================================================

  subroutine finalize_batch()

#ifdef MPI
    integer :: mpi_err ! MPI error code
#endif

    ! Reduce tallies onto master process and accumulate
    call time_tallies % start()
    call accumulate_tallies()
    call time_tallies % stop()

    ! Reset global tally results
    if (current_batch <= n_inactive) then
      global_tallies(:,:) = ZERO
      n_realizations = 0
    end if

    if (run_mode == MODE_EIGENVALUE) then
      ! Perform CMFD calculation if on
      if (cmfd_on) call execute_cmfd()
    end if

    ! Check_triggers
    if (master) call check_triggers()
#ifdef MPI
    call MPI_BCAST(satisfy_triggers, 1, MPI_LOGICAL, 0, &
         mpi_intracomm, mpi_err)
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

  end subroutine finalize_batch

!===============================================================================
! REPLAY_BATCH_HISTORY displays keff and entropy for each generation within a
! batch using data read from a state point file
!===============================================================================

  subroutine replay_batch_history

    ! Write message at beginning
    if (current_batch == 1) then
      call write_message("Replaying history from state point...", 6)
    end if

    if (run_mode == MODE_EIGENVALUE) then
      do current_gen = 1, gen_per_batch
        call calculate_average_keff()

        ! print out batch keff
        if (verbosity >= 7) then
          if (current_gen < gen_per_batch) then
            if (master) call print_generation()
          else
            if (master) call print_batch_keff()
          end if
        end if
      end do
    end if

    ! Write message at end
    if (current_batch == restart_batch) then
      call write_message("Resuming simulation...", 6)
    end if

  end subroutine replay_batch_history

!===============================================================================
! INITIALIZE_SIMULATION
!===============================================================================

  subroutine initialize_simulation()

    ! Set up tally procedure pointers
    call init_tally_routines()

    ! Determine how much work each processor should do
    call calculate_work()

    ! Allocate source bank, and for eigenvalue simulations also allocate the
    ! fission bank
    call allocate_banks()

    ! Allocate tally results arrays if they're not allocated yet
    call configure_tallies()

!$omp parallel
    ! Allocate array for microscopic cross section cache
    allocate(micro_xs(n_nuclides))

    ! Allocate array for matching filter bins
    allocate(filter_matches(n_filters))
!$omp end parallel

    ! If this is a restart run, load the state point data and binary source
    ! file
    if (restart_run) then
      call load_state_point()
    else
      call initialize_source()
    end if

    ! Display header
    if (master) then
      if (run_mode == MODE_FIXEDSOURCE) then
        call header("FIXED SOURCE TRANSPORT SIMULATION", 3)
      elseif (run_mode == MODE_EIGENVALUE) then
        call header("K EIGENVALUE SIMULATION", 3)
        if (verbosity >= 7) call print_columns()
      end if
    end if

  end subroutine initialize_simulation

!===============================================================================
! FINALIZE_SIMULATION calculates tally statistics, writes tallies, and displays
! execution time and results
!===============================================================================

  subroutine finalize_simulation()

#ifdef MPI
    integer    :: i       ! loop index for tallies
    integer    :: n       ! size of arrays
    integer    :: mpi_err  ! MPI error code
    integer(8) :: temp
    real(8)    :: tempr(3) ! temporary array for communication
#endif

!$omp parallel
    deallocate(micro_xs, filter_matches)
!$omp end parallel

    ! Increment total number of generations
    total_gen = total_gen + n_batches*gen_per_batch

    ! Start finalization timer
    call time_finalize % start()

#ifdef MPI
    ! Broadcast tally results so that each process has access to results
    if (allocated(tallies)) then
      do i = 1, size(tallies)
        n = size(tallies(i) % obj % results)
        call MPI_BCAST(tallies(i) % obj % results, n, MPI_DOUBLE, 0, &
             mpi_intracomm, mpi_err)
      end do
    end if

    ! Also broadcast global tally results
    n = size(global_tallies)
    call MPI_BCAST(global_tallies, n, MPI_DOUBLE, 0, mpi_intracomm, mpi_err)

    ! These guys are needed so that non-master processes can calculate the
    ! combined estimate of k-effective
    tempr(1) = k_col_abs
    tempr(2) = k_col_tra
    tempr(3) = k_abs_tra
    call MPI_BCAST(tempr, 3, MPI_REAL8, 0, mpi_intracomm, mpi_err)
    k_col_abs = tempr(1)
    k_col_tra = tempr(2)
    k_abs_tra = tempr(3)

    if (check_overlaps) then
      call MPI_REDUCE(overlap_check_cnt, temp, n_cells, MPI_INTEGER8, &
           MPI_SUM, 0, mpi_intracomm, mpi_err)
      overlap_check_cnt = temp
    end if
#endif

    ! Write tally results to tallies.out
    if (output_tallies .and. master) call write_tallies()

    ! Stop timers and show timing statistics
    call time_finalize%stop()
    call time_total%stop()
    if (master) then
      if (verbosity >= 6) call print_runtime()
      if (verbosity >= 4) call print_results()
      if (check_overlaps) call print_overlap_check()
    end if

  end subroutine finalize_simulation

!===============================================================================
! CALCULATE_WORK determines how many particles each processor should simulate
!===============================================================================

  subroutine calculate_work()

    integer    :: i         ! loop index
    integer    :: remainder ! Number of processors with one extra particle
    integer(8) :: i_bank    ! Running count of number of particles
    integer(8) :: min_work  ! Minimum number of particles on each proc
    integer(8) :: work_i    ! Number of particles on rank i

    if (.not. allocated(work_index)) allocate(work_index(0:n_procs))

    ! Determine minimum amount of particles to simulate on each processor
    min_work = n_particles/n_procs

    ! Determine number of processors that have one extra particle
    remainder = int(mod(n_particles, int(n_procs,8)), 4)

    i_bank = 0
    work_index(0) = 0
    do i = 0, n_procs - 1
      ! Number of particles for rank i
      if (i < remainder) then
        work_i = min_work + 1
      else
        work_i = min_work
      end if

      ! Set number of particles
      if (rank == i) work = work_i

      ! Set index into source bank for rank i
      i_bank = i_bank + work_i
      work_index(i+1) = i_bank
    end do

  end subroutine calculate_work

!===============================================================================
! ALLOCATE_BANKS allocates memory for the fission and source banks
!===============================================================================

  subroutine allocate_banks()

    integer :: alloc_err  ! allocation error code

    ! Allocate source bank
    if (allocated(source_bank)) deallocate(source_bank)
    allocate(source_bank(work), STAT=alloc_err)

    ! Check for allocation errors
    if (alloc_err /= 0) then
      call fatal_error("Failed to allocate source bank.")
    end if

    if (run_mode == MODE_EIGENVALUE) then

#ifdef _OPENMP
      ! If OpenMP is being used, each thread needs its own private fission
      ! bank. Since the private fission banks need to be combined at the end of
      ! a generation, there is also a 'master_fission_bank' that is used to
      ! collect the sites from each thread.

      n_threads = omp_get_max_threads()

!$omp parallel
      thread_id = omp_get_thread_num()

      if (allocated(fission_bank)) deallocate(fission_bank)
      if (thread_id == 0) then
        allocate(fission_bank(3*work))
      else
        allocate(fission_bank(3*work/n_threads))
      end if
!$omp end parallel
      if (allocated(master_fission_bank)) deallocate(master_fission_bank)
      allocate(master_fission_bank(3*work), STAT=alloc_err)
#else
      if (allocated(fission_bank)) deallocate(fission_bank)
      allocate(fission_bank(3*work), STAT=alloc_err)
#endif

      ! Check for allocation errors
      if (alloc_err /= 0) then
        call fatal_error("Failed to allocate fission bank.")
      end if
    end if

  end subroutine allocate_banks

end module simulation
