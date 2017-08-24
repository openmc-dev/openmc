module simulation_header

  use bank_header
  use constants

  implicit none

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Number of lost particles
  integer :: n_lost_particles

  real(8) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! EIGENVALUE SIMULATION VARIABLES

  integer    :: current_batch     ! current batch
  integer    :: current_gen       ! current generation within a batch
  integer    :: total_gen     = 0 ! total number of generations simulated

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  integer(8) :: work         ! number of particles per processor
  integer(8), allocatable :: work_index(:) ! starting index in source bank for each process
  integer(8) :: current_work ! index in source bank of current history simulated

  ! Temporary k-effective values
  real(8), allocatable :: k_generation(:) ! single-generation estimates of k
  real(8) :: keff = ONE       ! average k over active batches
  real(8) :: keff_std         ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength

  ! Shannon entropy
  real(8), allocatable :: entropy(:)     ! shannon entropy at each generation
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

end module simulation_header
