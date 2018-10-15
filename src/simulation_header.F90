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
  integer(C_INT), bind(C) :: n_lost_particles

  real(C_DOUBLE), bind(C) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! SIMULATION VARIABLES

  integer(C_INT), bind(C) :: current_batch     ! current batch
  integer(C_INT), bind(C) :: current_gen       ! current generation within a batch
  integer(C_INT), bind(C) :: total_gen         ! total number of generations simulated
  logical(C_BOOL), bind(C) :: simulation_initialized
  logical(C_BOOL), bind(C) :: need_depletion_rx ! need to calculate depletion reaction rx?

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  logical(C_BOOL), bind(C) :: satisfy_triggers  ! whether triggers are satisfied

  integer(C_INT64_T), bind(C) :: work         ! number of particles per processor
  integer(C_INT64_T), bind(C) :: current_work ! index in source bank of current history simulated

  ! ============================================================================
  ! K-EIGENVALUE SIMULATION VARIABLES

  ! Temporary k-effective values
  real(C_DOUBLE), bind(C) :: keff      ! average k over active batches
  real(C_DOUBLE), bind(C) :: keff_std  ! standard deviation of average k
  real(C_DOUBLE), bind(C) :: k_col_abs ! sum over batches of k_collision * k_absorption
  real(C_DOUBLE), bind(C) :: k_col_tra ! sum over batches of k_collision * k_tracklength
  real(C_DOUBLE), bind(C) :: k_abs_tra ! sum over batches of k_absorption * k_tracklength

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

#ifdef _OPENMP
  integer(C_INT), bind(C) :: n_threads      ! number of OpenMP threads
  integer(C_INT), bind(C) :: thread_id      ! ID of a given thread
#endif

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  integer(C_INT), bind(C) :: restart_batch

  logical(C_BOOL), bind(C) :: trace

!$omp threadprivate(trace, thread_id, current_work)

  interface
    subroutine entropy_clear() bind(C)
    end subroutine

    pure function overall_generation() result(gen) bind(C)
      import C_INT
      integer(C_INT) :: gen
    end function overall_generation

    function k_generation(i) result(k) bind(C)
      import C_DOUBLE, C_INT
      integer(C_INT), value :: i
      real(C_DOUBLE) :: k
    end function

    function k_generation_size() result(sz) bind(C)
      import C_INT
      integer(C_INT) :: sz
    end function

    subroutine k_generation_clear() bind(C)
    end subroutine

    subroutine k_generation_reserve(i) bind(C)
      import C_INT
      integer(C_INT), value :: i
    end subroutine

    function work_index(rank) result(i) bind(C)
      import C_INT, C_INT64_T
      integer(C_INT), value :: rank
      integer(C_INT64_T) :: i
    end function
  end interface

contains

!===============================================================================
! FREE_MEMORY_SIMULATION deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_simulation()

    call k_generation_clear()
    call entropy_clear()
  end subroutine free_memory_simulation

end module simulation_header
