module simulation_header

  use, intrinsic :: ISO_C_BINDING

  use bank_header
  use constants
  use settings, only: gen_per_batch
  use stl_vector, only: VectorReal

  implicit none

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Number of lost particles
  integer(C_INT), bind(C) :: n_lost_particles = 0

  real(8) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! SIMULATION VARIABLES

  integer(C_INT), bind(C) :: current_batch     ! current batch
  integer(C_INT), bind(C) :: current_gen       ! current generation within a batch
  integer(C_INT), bind(C) :: total_gen     = 0 ! total number of generations simulated
  logical(C_BOOL), bind(C) :: simulation_initialized = .false.
  logical :: need_depletion_rx ! need to calculate depletion reaction rx?

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  integer(C_INT64_T), bind(C) :: work         ! number of particles per processor
  integer(C_INT64_T), allocatable :: work_index(:) ! starting index in source bank for each process
  integer(C_INT64_T), bind(C) :: current_work ! index in source bank of current history simulated

  ! ============================================================================
  ! K-EIGENVALUE SIMULATION VARIABLES

  ! Temporary k-effective values
  type(VectorReal) :: k_generation ! single-generation estimates of k
  real(C_DOUBLE), bind(C, name='openmc_keff')     :: keff = ONE  ! average k over active batches
  real(C_DOUBLE), bind(C, name='openmc_keff_std') :: keff_std    ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

#ifdef _OPENMP
  integer(C_INT), bind(C) :: n_threads = NONE      ! number of OpenMP threads
  integer :: thread_id             ! ID of a given thread
#endif

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  integer :: restart_batch

  logical(C_BOOL), bind(C) :: trace

!$omp threadprivate(trace, thread_id, current_work)

  interface
    subroutine entropy_clear() bind(C)
    end subroutine
  end interface


contains

!===============================================================================
! OVERALL_GENERATION determines the overall generation number
!===============================================================================

  pure function overall_generation() result(gen) bind(C)
    integer(C_INT) :: gen
    gen = gen_per_batch*(current_batch - 1) + current_gen
  end function overall_generation

!===============================================================================
! FREE_MEMORY_SIMULATION deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_simulation()

    if (allocated(work_index)) deallocate(work_index)

    call k_generation % clear()
    call k_generation % shrink_to_fit()
    call entropy_clear()
  end subroutine free_memory_simulation

end module simulation_header
