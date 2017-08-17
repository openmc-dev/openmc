module global

  use, intrinsic :: ISO_C_BINDING

#ifdef MPIF08
  use mpi_f08
#endif

  use bank_header,      only: Bank
  use cmfd_header
  use constants
  use stl_vector,       only: VectorInt

  ! Inherit module variables from other modules
  use geometry_header
  use material_header
  use mesh_header
  use mgxs_header
  use nuclide_header
  use plot_header
  use sab_header
  use surface_header
  use tally_filter_header
  use tally_header
  use timer_header
  use trigger_header
  use volume_header

  ! Inherit settings
  use settings

  implicit none

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Number of lost particles
  integer :: n_lost_particles

  real(8) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  ! Active tally lists
  type(VectorInt) :: active_analog_tallies
  type(VectorInt) :: active_tracklength_tallies
  type(VectorInt) :: active_current_tallies
  type(VectorInt) :: active_collision_tallies
  type(VectorInt) :: active_tallies
  type(VectorInt) :: active_surface_tallies

  ! Normalization for statistics
  integer :: n_realizations = 0 ! # of independent realizations
  real(8) :: total_weight       ! total starting particle weight in realization

  ! ============================================================================
  ! EIGENVALUE SIMULATION VARIABLES

  integer    :: current_batch     ! current batch
  integer    :: current_gen       ! current generation within a batch
  integer    :: total_gen     = 0 ! total number of generations simulated

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  ! Source and fission bank
  type(Bank), allocatable, target :: source_bank(:)
  type(Bank), allocatable, target :: fission_bank(:)
#ifdef _OPENMP
  type(Bank), allocatable, target :: master_fission_bank(:)
#endif
  integer(8) :: n_bank       ! # of sites in fission bank
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
  real(8), allocatable :: entropy(:)         ! shannon entropy at each generation
  real(8), allocatable :: entropy_p(:,:,:,:) ! % of source sites in each cell

  ! Uniform fission source weighting
  real(8), allocatable :: source_frac(:,:,:,:)

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

#ifdef _OPENMP
  integer :: n_threads = NONE      ! number of OpenMP threads
  integer :: thread_id             ! ID of a given thread
#endif

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  integer :: restart_batch

  ! Flag for enabling cell overlap checking during transport
  integer(8), allocatable  :: overlap_check_cnt(:)

  logical    :: trace

  ! Number of distribcell maps
  integer :: n_maps

!$omp threadprivate(fission_bank, n_bank, trace, thread_id, current_work)

contains

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory()

    integer :: i ! Loop Index

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(plots)) deallocate(plots)
    if (allocated(volume_calcs)) deallocate(volume_calcs)

    ! Deallocate geometry debugging information
    if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
    ! First call the clear routines
      do i = 1, size(nuclides)
        call nuclides(i) % clear()
      end do
      deallocate(nuclides)
    end if
    if (allocated(libraries)) deallocate(libraries)

    if (allocated(res_scat_nuclides)) deallocate(res_scat_nuclides)

    if (allocated(nuclides_MG)) deallocate(nuclides_MG)

    if (allocated(macro_xs)) deallocate(macro_xs)

    if (allocated(sab_tables)) deallocate(sab_tables)

    ! Deallocate external source
    if (allocated(external_source)) deallocate(external_source)

    ! Deallocate k and entropy
    if (allocated(k_generation)) deallocate(k_generation)
    if (allocated(entropy)) deallocate(entropy)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate tally-related arrays
    if (allocated(global_tallies)) deallocate(global_tallies)
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(filters)) deallocate(filters)
    if (allocated(tallies)) deallocate(tallies)

    ! Deallocate fission and source bank and entropy
!$omp parallel
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(tally_derivs)) deallocate(tally_derivs)
!$omp end parallel
#ifdef _OPENMP
    if (allocated(master_fission_bank)) deallocate(master_fission_bank)
#endif
    if (allocated(source_bank)) deallocate(source_bank)

    ! Deallocate array of work indices
    if (allocated(work_index)) deallocate(work_index)

    if (allocated(energy_bins)) deallocate(energy_bins)
    if (allocated(energy_bin_avg)) deallocate(energy_bin_avg)

    ! Deallocate CMFD
    call deallocate_cmfd(cmfd)

    ! Deallocate tally node lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_surface_tallies % clear()
    call active_tallies % clear()

    ! Deallocate track_identifiers
    if (allocated(track_identifiers)) deallocate(track_identifiers)

    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()
    call material_dict % clear()
    call mesh_dict % clear()
    call filter_dict % clear()
    call tally_dict % clear()
    call plot_dict % clear()
    call nuclide_dict % clear()
    call sab_dict % clear()
    call library_dict % clear()

    ! Clear statepoint and sourcepoint batch set
    call statepoint_batch % clear()
    call sourcepoint_batch % clear()

    ! Deallocate entropy mesh
    if (associated(entropy_mesh)) then
      if (allocated(entropy_mesh % lower_left)) &
           deallocate(entropy_mesh % lower_left)
      if (allocated(entropy_mesh % upper_right)) &
           deallocate(entropy_mesh % upper_right)
      if (allocated(entropy_mesh % width)) deallocate(entropy_mesh % width)
      deallocate(entropy_mesh)
    end if

    ! Deallocate ufs
    if (allocated(source_frac)) deallocate(source_frac)
    if (associated(ufs_mesh)) then
        if (allocated(ufs_mesh % lower_left)) deallocate(ufs_mesh % lower_left)
        if (allocated(ufs_mesh % upper_right)) &
             deallocate(ufs_mesh % upper_right)
        if (allocated(ufs_mesh % width)) deallocate(ufs_mesh % width)
        deallocate(ufs_mesh)
    end if

  end subroutine free_memory

!===============================================================================
! OVERALL_GENERATION determines the overall generation number
!===============================================================================

  pure function overall_generation() result(gen)
    integer :: gen
    gen = gen_per_batch*(current_batch - 1) + current_gen
  end function overall_generation

end module global
