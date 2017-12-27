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
  integer :: n_lost_particles = 0

  real(8) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! SIMULATION VARIABLES

  integer :: current_batch     ! current batch
  integer :: current_gen       ! current generation within a batch
  integer :: total_gen     = 0 ! total number of generations simulated
  logical(C_BOOL), bind(C) :: simulation_initialized = .false.
  logical :: need_depletion_rx ! need to calculate depletion reaction rx?

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  integer(8) :: work         ! number of particles per processor
  integer(8), allocatable :: work_index(:) ! starting index in source bank for each process
  integer(8) :: current_work ! index in source bank of current history simulated

  ! ============================================================================
  ! K-EIGENVALUE SIMULATION VARIABLES

  ! Temporary k-effective values
  type(VectorReal) :: k_generation ! single-generation estimates of k
  real(C_DOUBLE), bind(C) :: keff = ONE  ! average k over active batches
  real(C_DOUBLE), bind(C) :: keff_std    ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength

  ! Shannon entropy
  type(VectorReal)     :: entropy        ! shannon entropy at each generation
  real(8), allocatable :: entropy_p(:,:) ! % of source sites in each cell

  ! Uniform fission source weighting
  real(8), allocatable :: source_frac(:,:)

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

#ifdef _OPENMP
  integer :: n_threads = NONE      ! number of OpenMP threads
  integer :: thread_id             ! ID of a given thread
#endif

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  integer :: restart_batch

  ! Flag for enabling cell overlap checking during transport
  integer(8), allocatable  :: overlap_check_cnt(:)

  logical :: trace

  ! Number of distribcell maps
  integer :: n_maps

!$omp threadprivate(trace, thread_id, current_work)

contains

!===============================================================================
! OVERALL_GENERATION determines the overall generation number
!===============================================================================

  pure function overall_generation() result(gen)
    integer :: gen
    gen = gen_per_batch*(current_batch - 1) + current_gen
  end function overall_generation

!===============================================================================
! FREE_MEMORY_SIMULATION deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_simulation()
    if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)
    if (allocated(entropy_p)) deallocate(entropy_p)
    if (allocated(source_frac)) deallocate(source_frac)
    if (allocated(work_index)) deallocate(work_index)

    call k_generation % clear()
    call k_generation % shrink_to_fit()
    call entropy % clear()
    call entropy % shrink_to_fit()
  end subroutine free_memory_simulation

end module simulation_header
