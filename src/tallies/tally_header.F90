module tally_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use dict_header,         only: DictIntInt
  use hdf5_interface,      only: HID_T, HSIZE_T
  use message_passing,     only: n_procs, master
  use nuclide_header,      only: nuclide_dict
  use settings,            only: reduce_tallies, run_mode
  use stl_vector,          only: VectorInt
  use string,              only: to_lower, to_f_string, str_to_int, to_str
  use tally_filter_header, only: TallyFilterContainer, filters, n_filters
  use tally_filter

  implicit none

  interface
    function openmc_tally_set_filters(index, n, filter_indices) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT), value, intent(in) :: n
      integer(C_INT32_T), intent(in) :: filter_indices(n)
      integer(C_INT) :: err
    end function

    function openmc_tally_set_active(index, active) result(err) bind(C)
      import C_INT32_T, C_BOOL, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      logical(C_BOOL),    value, intent(in) :: active
      integer(C_INT) :: err
    end function

    function openmc_tally_get_active(index, active) result(err) bind(C)
      import C_INT32_T, C_BOOL, C_INT
      integer(C_INT32_T), value    :: index
      logical(C_BOOL), intent(out) :: active
      integer(C_INT) :: err
    end function

    function openmc_tally_get_nuclides(index, nuclides, n) result(err) bind(C)
      import C_INT32_T, C_PTR, C_INT
      integer(C_INT32_T), value :: index
      type(C_PTR), intent(out) :: nuclides
      integer(C_INT), intent(out) :: n
      integer(C_INT) :: err
    end function openmc_tally_get_nuclides

    function active_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function

    function active_tallies_data(i) result(tally) bind(C)
      import C_INT
      integer(C_INT), value :: i
      integer(C_INT) :: tally
    end function

    function active_analog_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function

    function active_tracklength_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function

    function active_collision_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function

    function active_meshsurf_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function

    function active_surface_tallies_size() result(size) bind(C)
      import C_INT
      integer(C_INT) :: size
    end function
  end interface

!===============================================================================
! TALLYOBJECT describes a user-specified tally. The region of phase space to
! tally in is given by the TallyFilters and the results are stored in a
! TallyResult array.
!===============================================================================

  type, public :: TallyObject
    type(C_PTR) :: ptr

    ! Basic data

    character(len=104) :: name = "" ! user-defined name
    real(8) :: volume               ! volume of region
    logical :: depletion_rx = .false. ! has depletion reactions, e.g. (n,2n)

    ! Values to score, e.g. flux, absorption, etc.
    integer              :: n_score_bins = 0
    integer, allocatable :: score_bins(:)

    ! Results for each bin -- the first dimension of the array is for scores
    ! (e.g. flux, total reaction rate, fission reaction rate, etc.) and the
    ! second dimension of the array is for the combination of filters
    ! (e.g. specific cell, specific energy group, etc.)

    integer :: total_score_bins
    real(C_DOUBLE), allocatable :: results(:,:,:)

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

  contains
    procedure :: accumulate => tally_accumulate
    procedure :: allocate_results => tally_allocate_results
    procedure :: read_results_hdf5 => tally_read_results_hdf5
    procedure :: write_results_hdf5 => tally_write_results_hdf5
    procedure :: id => tally_get_id
    procedure :: set_id => tally_set_id
    procedure :: type => tally_get_type
    procedure :: set_type => tally_set_type
    procedure :: estimator => tally_get_estimator
    procedure :: set_estimator => tally_set_estimator
    procedure :: n_filters => tally_get_n_filters
    procedure :: filter => tally_get_filter
    procedure :: stride => tally_get_stride
    procedure :: n_filter_bins => tally_get_n_filter_bins
    procedure :: n_nuclide_bins => tally_get_n_nuclide_bins
    procedure :: nuclide_bins => tally_get_nuclide_bins
    procedure :: set_nuclide_bins => tally_set_nuclide_bins
    procedure :: all_nuclides => tally_get_all_nuclides
    procedure :: energyout_filter => tally_get_energyout_filter
    procedure :: delayedgroup_filter => tally_get_delayedgroup_filter
    procedure :: deriv => tally_get_deriv
    procedure :: set_deriv => tally_set_deriv
  end type TallyObject

  type, public :: TallyContainer
    class(TallyObject), allocatable :: obj
  end type TallyContainer

  integer(C_INT32_T), public, bind(C) :: n_tallies = 0 ! # of tallies

  type(TallyContainer),  public, allocatable, target :: tallies(:)

  ! Dictionary that maps user IDs to indices in 'tallies'
  type(DictIntInt), public :: tally_dict

  ! The largest tally ID that has been specified in the system.  This is useful
  ! in case the code needs to find an ID for a new tally.
  integer :: largest_tally_id

  ! Global tallies
  !   1) collision estimate of k-eff
  !   2) absorption estimate of k-eff
  !   3) track-length estimate of k-eff
  !   4) leakage fraction

  real(C_DOUBLE), public, allocatable, target :: global_tallies(:,:)

  ! It is possible to protect accumulate operations on global tallies by using
  ! an atomic update. However, when multiple threads accumulate to the same
  ! global tally, it can cause a higher cache miss rate due to
  ! invalidation. Thus, we use threadprivate variables to accumulate global
  ! tallies and then reduce at the end of a generation.
  real(C_DOUBLE), public, bind(C) :: global_tally_collision
  real(C_DOUBLE), public, bind(C) :: global_tally_absorption
  real(C_DOUBLE), public, bind(C) :: global_tally_tracklength
  real(C_DOUBLE), public, bind(C) :: global_tally_leakage
!$omp threadprivate(global_tally_collision, global_tally_absorption, &
!$omp&              global_tally_tracklength, global_tally_leakage)

  ! Normalization for statistics
  integer(C_INT32_T), public, bind(C) :: n_realizations = 0 ! # of independent realizations
  real(C_DOUBLE), public, bind(C) :: total_weight       ! total starting particle weight in realization

contains

!===============================================================================
! ACCUMULATE_TALLY
!===============================================================================

  subroutine tally_accumulate(this)
    class(TallyObject), intent(inout) :: this

    integer :: i, j
    real(C_DOUBLE) :: val
    real(C_DOUBLE) :: total_source

    interface
      function total_source_strength() result(strength) bind(C)
        import C_DOUBLE
        real(C_DOUBLE) :: strength
      end function
    end interface

    ! Increment number of realizations
    if (reduce_tallies) then
      this % n_realizations = this % n_realizations + 1
    else
      this % n_realizations = this % n_realizations + n_procs
    end if

    if (master .or. (.not. reduce_tallies)) then
      ! Calculate total source strength for normalization
      if (run_mode == MODE_FIXEDSOURCE) then
        total_source = total_source_strength()
      else
        total_source = ONE
      end if

      ! Accumulate each result
      do j = 1, size(this % results, 3)
        do i = 1, size(this % results, 2)
          val = this % results(RESULT_VALUE, i, j)/total_weight * total_source
          this % results(RESULT_VALUE, i, j) = ZERO

          this % results(RESULT_SUM, i, j) = &
              this % results(RESULT_SUM, i, j) + val
          this % results(RESULT_SUM_SQ, i, j) = &
              this % results(RESULT_SUM_SQ, i, j) + val*val
        end do
      end do
    end if
  end subroutine tally_accumulate

  subroutine tally_write_results_hdf5(this, group_id)
    class(TallyObject), intent(in) :: this
    integer(HID_T),     intent(in) :: group_id

    integer(HSIZE_T) :: n_filter, n_score
    interface
      subroutine write_tally_results(group_id, n_filter, n_score, results) &
           bind(C)
        import HID_T, HSIZE_T, C_DOUBLE
        integer(HID_T), value :: group_id
        integer(HSIZE_T), value :: n_filter
        integer(HSIZE_T), value :: n_score
        real(C_DOUBLE), intent(in) :: results(*)
      end subroutine write_tally_results
    end interface

    n_filter = size(this % results, 3)
    n_score = size(this % results, 2)
    call write_tally_results(group_id, n_filter, n_score, this % results)
  end subroutine tally_write_results_hdf5

  subroutine tally_read_results_hdf5(this, group_id)
    class(TallyObject), intent(inout) :: this
    integer(HID_T),     intent(in) :: group_id

    integer(HSIZE_T) :: n_filter, n_score
    interface
      subroutine read_tally_results(group_id, n_filter, n_score, results) &
           bind(C)
        import HID_T, HSIZE_T, C_DOUBLE
        integer(HID_T), value :: group_id
        integer(HSIZE_T), value :: n_filter
        integer(HSIZE_T), value :: n_score
        real(C_DOUBLE), intent(out) :: results(*)
      end subroutine read_tally_results
    end interface

    n_filter = size(this % results, 3)
    n_score = size(this % results, 2)
    call read_tally_results(group_id, n_filter, n_score, this % results)
  end subroutine tally_read_results_hdf5

!===============================================================================
! ALLOCATE_RESULTS allocates and initializes the results component of the
! TallyObject derived type
!===============================================================================

  subroutine tally_allocate_results(this)
    class(TallyObject), intent(inout) :: this
    integer, parameter :: default_nuclide_bins(1) = (/-1/)

    ! If no nuclides were specified, add a single bin for total material
    if (this % n_nuclide_bins() == 0) then
      call this % set_nuclide_bins(default_nuclide_bins)
    end if

    ! Set total number of filter and scoring bins
    this % total_score_bins = this % n_score_bins * this % n_nuclide_bins()

    if (allocated(this % results)) then
      ! If results was already allocated but shape is wrong, then reallocate it
      ! to the correct shape
      if (this % total_score_bins /= size(this % results, 2) .or. &
           this % n_filter_bins() /= size(this % results, 3)) then
        deallocate(this % results)
        allocate(this % results(3, this % total_score_bins, this % n_filter_bins()))
      end if
    else
      allocate(this % results(3, this % total_score_bins, this % n_filter_bins()))
    end if

  end subroutine tally_allocate_results

  function tally_get_id(this) result(t)
    class(TallyObject) :: this
    integer(C_INT) :: t
    interface
      function tally_get_id_c(tally) result(t) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: t
      end function
    end interface
    t = tally_get_id_c(this % ptr)
  end function

  subroutine tally_set_id(this, t)
    class(TallyObject) :: this
    integer(C_INT) :: t
    interface
      subroutine tally_set_id_c(tally, t) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: t
      end subroutine
    end interface
    call tally_set_id_c(this % ptr, t)
  end subroutine

  function tally_get_type(this) result(t)
    class(TallyObject) :: this
    integer(C_INT) :: t
    interface
      function tally_get_type_c(tally) result(t) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: t
      end function
    end interface
    t = tally_get_type_c(this % ptr)
  end function

  subroutine tally_set_type(this, t)
    class(TallyObject) :: this
    integer(C_INT) :: t
    interface
      subroutine tally_set_type_c(tally, t) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: t
      end subroutine
    end interface
    call tally_set_type_c(this % ptr, t)
  end subroutine

  function tally_get_estimator(this) result(e)
    class(TallyObject) :: this
    integer(C_INT) :: e
    interface
      function tally_get_estimator_c(tally) result(e) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: e
      end function
    end interface
    e = tally_get_estimator_c(this % ptr)
  end function

  subroutine tally_set_estimator(this, e)
    class(TallyObject) :: this
    integer(C_INT) :: e
    interface
      subroutine tally_set_estimator_c(tally, e) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: e
      end subroutine
    end interface
    call tally_set_estimator_c(this % ptr, e)
  end subroutine

  function tally_get_n_filters(this) result(n)
    class(TallyObject) :: this
    integer(C_INT) :: n
    interface
      function tally_get_n_filters_c(tally) result(n) bind(C)
        import C_PTR, C_INT, C_INT32_T
        type(C_PTR), value :: tally
        integer(C_INT) :: n
      end function
    end interface
    n = tally_get_n_filters_c(this % ptr)
  end function

  function tally_get_filter(this, i) result(filt)
    class(TallyObject) :: this
    integer(C_INT) :: i
    integer(C_INT32_T) :: filt
    interface
      function tally_get_filter_c(tally, i) result(filt) bind(C)
        import C_PTR, C_INT, C_INT32_T
        type(C_PTR), value :: tally
        integer(C_INT), value :: i
        integer(C_INT32_T) :: filt
      end function
    end interface
    filt = tally_get_filter_c(this % ptr, i-1)
  end function

  function tally_get_stride(this, i) result(stride)
    class(TallyObject) :: this
    integer(C_INT) :: i
    integer(C_INT32_T) :: stride
    interface
      function tally_get_stride_c(tally, i) result(stride) bind(C)
        import C_PTR, C_INT, C_INT32_T
        type(C_PTR), value :: tally
        integer(C_INT), value :: i
        integer(C_INT32_T) :: stride
      end function
    end interface
    stride = tally_get_stride_c(this % ptr, i-1)
  end function

  function tally_get_n_filter_bins(this) result(n_filter_bins)
    class(TallyObject) :: this
    integer(C_INT32_T) :: n_filter_bins
    interface
      function tally_get_n_filter_bins_c(tally) result(n_filter_bins) bind(C)
        import C_PTR, C_INT, C_INT32_T
        type(C_PTR), value :: tally
        integer(C_INT32_T) :: n_filter_bins
      end function
    end interface
    n_filter_bins = tally_get_n_filter_bins_c(this % ptr)
  end function

  function tally_get_n_nuclide_bins(this) result(n)
    class(TallyObject) :: this
    integer(C_INT) :: n
    interface
      function tally_get_n_nuclide_bins_c(tally) result(n) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: n
      end function
    end interface
    n = tally_get_n_nuclide_bins_c(this % ptr)
  end function

  function tally_get_nuclide_bins(this, i) result(nuclide)
    class(TallyObject) :: this
    integer(C_INT) :: i, nuclide
    interface
      function tally_get_nuclide_bins_c(tally, i) result(nuclide) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: i
        integer(C_INT) :: nuclide
      end function
    end interface
    nuclide = tally_get_nuclide_bins_c(this % ptr, i)
  end function

  subroutine tally_set_nuclide_bins(this, bins, all_nuclides)
    class(TallyObject) :: this
    integer :: bins(:)
    logical, optional :: all_nuclides
    logical(C_BOOL) :: all_
    interface
      subroutine tally_set_nuclide_bins_c(this, n, bins, all_nuclides) bind(C)
        import C_PTR, C_INT, C_BOOL
        type(C_PTR), value :: this
        integer(C_INT), value :: n
        integer(C_INT) :: bins(n)
        logical(C_BOOL), value :: all_nuclides
      end subroutine
    end interface
    if (present(all_nuclides)) then
      all_ = logical(all_nuclides, kind=C_BOOL)
    else
      all_ = .false._C_BOOL
    end if
    call tally_set_nuclide_bins_c(this % ptr, size(bins), bins, all_)
  end subroutine

  function tally_get_all_nuclides(this) result(all_nuc)
    class(TallyObject) :: this
    logical(C_BOOL) :: all_nuc
    interface
      function tally_get_all_nuclides_c(tally) result(all_nuc) bind(C)
        import C_PTR, C_BOOL
        type(C_PTR), value :: tally
        logical(C_BOOL) :: all_nuc
      end function
    end interface
    all_nuc = tally_get_all_nuclides_c(this % ptr)
  end function

  function tally_get_energyout_filter(this) result(filt)
    class(TallyObject) :: this
    integer(C_INT) :: filt
    interface
      function tally_get_energyout_filter_c(tally) result(filt) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: filt
      end function
    end interface
    filt = tally_get_energyout_filter_c(this % ptr)
  end function

  function tally_get_delayedgroup_filter(this) result(filt)
    class(TallyObject) :: this
    integer(C_INT) :: filt
    interface
      function tally_get_delayedgroup_filter_c(tally) result(filt) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: filt
      end function
    end interface
    filt = tally_get_delayedgroup_filter_c(this % ptr)
  end function

  function tally_get_deriv(this) result(deriv)
    class(TallyObject) :: this
    integer(C_INT) :: deriv
    interface
      function tally_get_deriv_c(tally) result(deriv) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: deriv
      end function
    end interface
    deriv = tally_get_deriv_c(this % ptr)
  end function

  subroutine tally_set_deriv(this, deriv)
    class(TallyObject) :: this
    integer(C_INT) :: deriv
    interface
      subroutine tally_set_deriv_c(tally, deriv) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: deriv
      end subroutine
    end interface
    call tally_set_deriv_c(this % ptr, deriv)
  end subroutine

!===============================================================================
! CONFIGURE_TALLIES initializes several data structures related to tallies. This
! is called after the basic tally data has already been read from the
! tallies.xml file.
!===============================================================================

  subroutine allocate_tally_results() bind(C)

    integer :: i

    ! Allocate global tallies
    if (.not. allocated(global_tallies)) then
      allocate(global_tallies(3, N_GLOBAL_TALLIES))
    end if

    ! Allocate results arrays for tallies
    do i = 1, n_tallies
      call tallies(i) % obj % allocate_results()
    end do

  end subroutine allocate_tally_results

!===============================================================================
! FREE_MEMORY_TALLY deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_tally()
    interface
      subroutine free_memory_tally_c() bind(C)
      end subroutine free_memory_tally_c
    end interface

    call free_memory_tally_c()

    n_tallies = 0
    if (allocated(tallies)) deallocate(tallies)
    call tally_dict % clear()
    largest_tally_id = 0

    if (allocated(global_tallies)) deallocate(global_tallies)
  end subroutine free_memory_tally

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_tallies(n, index_start, index_end) result(err) bind(C)
    ! Extend the tallies array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    integer :: i
    type(TallyContainer), allocatable :: temp(:) ! temporary tallies array

    interface
      subroutine extend_tallies_c() bind(C)
      end subroutine
    end interface

    ! Extend the C++ tallies array first
    call extend_tallies_c()

    if (n_tallies == 0) then
      ! Allocate tallies array
      allocate(tallies(n))
    else
      ! Allocate tallies array with increased size
      allocate(temp(n_tallies + n))

      ! Move original tallies to temporary array
      do i = 1, n_tallies
        call move_alloc(tallies(i) % obj, temp(i) % obj)
      end do

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=tallies)
    end if

    ! Return indices in tallies array
    if (present(index_start)) index_start = n_tallies + 1
    if (present(index_end)) index_end = n_tallies + n
    n_tallies = n_tallies + n

    err = 0
  end function openmc_extend_tallies


  function openmc_get_tally_index(id, index) result(err) bind(C)
    ! Returns the index in the tallies array of a tally with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(tallies)) then
      if (tally_dict % has(id)) then
        index = tally_dict % get(id)
        err = 0
      else
        err = E_INVALID_ID
        call set_errmsg("No tally exists with ID=" // trim(to_str(id)) // ".")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory has not been allocated for tallies.")
    end if
  end function openmc_get_tally_index


  function openmc_global_tallies(ptr) result(err) bind(C)
    type(C_PTR), intent(out) :: ptr
    integer(C_INT) :: err

    if (.not. allocated(global_tallies)) then
      err = E_ALLOCATE
      call set_errmsg("Global tallies have not been allocated yet.")
    else
      err = 0
      ptr = C_LOC(global_tallies)
    end if
  end function openmc_global_tallies


  function openmc_tally_get_estimator(index, estimator) result(err) bind(C)
    ! Return the type of estimator of a tally
    integer(C_INT32_T), value    :: index
    integer(C_INT32_T), intent(out) :: estimator
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      estimator = tallies(index) % obj % estimator()
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_get_estimator


  function openmc_tally_get_id(index, id) result(err) bind(C)
    ! Return the ID of a tally
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      id = tallies(index) % obj % id()
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_get_id


  function openmc_tally_get_n_realizations(index, n) result(err) bind(C)
    ! Return the number of realizations for a tally
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      n = tallies(index) % obj % n_realizations
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_get_n_realizations


  function openmc_tally_get_scores(index, scores, n) result(err) bind(C)
    ! Return the list of nuclides assigned to a tally
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: scores
    integer(C_INT), intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index) % obj)
        if (allocated(t % score_bins)) then
          scores = C_LOC(t % score_bins(1))
          n = size(t % score_bins)
          err = 0
        else
          err = E_ALLOCATE
          call set_errmsg("Tally scores have not been allocated yet.")
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_get_scores


  function openmc_tally_reset(index) result(err) bind(C)
    ! Reset tally results and number of realizations
    integer(C_INT32_T), intent(in), value :: index
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index) % obj)
        t % n_realizations = 0
        if (allocated(t % results)) t % results(:, :, :) = ZERO
        err = 0
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_reset


  function openmc_tally_results(index, ptr, shape_) result(err) bind(C)
    ! Returns a pointer to a tally results array along with its shape. This
    ! allows a user to obtain in-memory tally results from Python directly.
    integer(C_INT32_T), intent(in), value :: index
    type(C_PTR),        intent(out) :: ptr
    integer(C_SIZE_T),  intent(out) :: shape_(3)
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index) % obj)
        if (allocated(t % results)) then
          ptr = C_LOC(t % results(1,1,1))

          ! Note that shape is reversed since it is assumed to be used from
          ! C/C++ code
          shape_(1) = size(t % results, 3)
          shape_(2) = size(t % results, 2)
          shape_(3) = size(t % results, 1)
          err = 0
        else
          err = E_ALLOCATE
          call set_errmsg("Tally results have not been allocated yet.")
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_results


  function openmc_tally_set_estimator(index, estimator) result(err) bind(C)
    ! Set the type of estimator a tally
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: estimator(*)
    integer(C_INT) :: err

    character(:), allocatable :: estimator_

    ! Convert C string to Fortran string
    estimator_ = to_f_string(estimator)

    err = 0
    if (index >= 1 .and. index <= size(tallies)) then
      select case (estimator_)
      case ('analog')
        call tallies(index) % obj % set_estimator(ESTIMATOR_ANALOG)
      case ('tracklength')
        call tallies(index) % obj % set_estimator(ESTIMATOR_TRACKLENGTH)
      case ('collision')
        call tallies(index) % obj % set_estimator(ESTIMATOR_COLLISION)
      case default
        err = E_INVALID_ARGUMENT
        call set_errmsg("Unknown tally estimator: " // trim(estimator_))
      end select
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in tally array is out of bounds.")
    end if
  end function openmc_tally_set_estimator


  function openmc_tally_set_id(index, id) result(err) bind(C)
    ! Set the ID of a tally
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_tallies) then
      if (allocated(tallies(index) % obj)) then
        if (tally_dict % has(id)) then
          call set_errmsg("Two or more tallies use the same unique ID: " &
               // to_str(id))
          err = E_INVALID_ID
        else
          call tallies(index) % obj % set_id(id)
          call tally_dict % set(id, index)
          if (id > largest_tally_id) largest_tally_id = id

          err = 0
        end if
      else
        err = E_ALLOCATE
        call set_errmsg("Tally type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_set_id


  function openmc_tally_set_nuclides(index, n, nuclides) result(err) bind(C)
    ! Sets the nuclides in the tally which results should be scored for
    integer(C_INT32_T), value  :: index
    integer(C_INT), value      :: n
    type(C_PTR),    intent(in) :: nuclides(n)
    integer(C_INT) :: err

    integer, allocatable :: bins(:)
    integer :: i
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: nuclide_

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index) % obj)
        allocate(bins(n))

        do i = 1, n
          ! Convert C string to Fortran string
          call c_f_pointer(nuclides(i), string, [10])
          nuclide_ = to_lower(to_f_string(string))

          select case (nuclide_)
          case ('total')
            bins(i) = -1
          case default
            if (nuclide_dict % has(nuclide_)) then
              bins(i) = nuclide_dict % get(nuclide_)
            else
              err = E_DATA
              call set_errmsg("Nuclide '" // trim(to_f_string(string)) // &
                   "' has not been loaded yet.")
              return
            end if
          end select
        end do

        call t % set_nuclide_bins(bins)

        err = 0
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_set_nuclides


  function openmc_tally_set_scores(index, n, scores) result(err) bind(C)
    ! Sets the scores in the tally
    integer(C_INT32_T), value  :: index
    integer(C_INT), value      :: n
    type(C_PTR),    intent(in) :: scores(n)
    integer(C_INT) :: err

    integer :: i
    integer :: MT
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: score_
    logical :: depletion_rx

    err = E_UNASSIGNED
    depletion_rx = .false.
    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index) % obj)
        if (allocated(t % score_bins)) deallocate(t % score_bins)
        allocate(t % score_bins(n))
        t % n_score_bins = n

        do i = 1, n
          ! Convert C string to Fortran string
          call c_f_pointer(scores(i), string, [20])
          score_ = to_lower(to_f_string(string))

          select case (score_)
          case ('flux')
            t % score_bins(i) = SCORE_FLUX
          case ('total', '(n,total)')
            t % score_bins(i) = SCORE_TOTAL
          case ('scatter')
            t % score_bins(i) = SCORE_SCATTER
          case ('nu-scatter')
            t % score_bins(i) = SCORE_NU_SCATTER
          case ('(n,2n)')
            t % score_bins(i) = N_2N
            depletion_rx = .true.
          case ('(n,3n)')
            t % score_bins(i) = N_3N
            depletion_rx = .true.
          case ('(n,4n)')
            t % score_bins(i) = N_4N
            depletion_rx = .true.
          case ('absorption')
            t % score_bins(i) = SCORE_ABSORPTION
          case ('fission', '18')
            t % score_bins(i) = SCORE_FISSION
          case ('nu-fission')
            t % score_bins(i) = SCORE_NU_FISSION
          case ('decay-rate')
            t % score_bins(i) = SCORE_DECAY_RATE
          case ('delayed-nu-fission')
            t % score_bins(i) = SCORE_DELAYED_NU_FISSION
          case ('prompt-nu-fission')
            t % score_bins(i) = SCORE_PROMPT_NU_FISSION
          case ('kappa-fission')
            t % score_bins(i) = SCORE_KAPPA_FISSION
          case ('inverse-velocity')
            t % score_bins(i) = SCORE_INVERSE_VELOCITY
          case ('fission-q-prompt')
            t % score_bins(i) = SCORE_FISS_Q_PROMPT
          case ('fission-q-recoverable')
            t % score_bins(i) = SCORE_FISS_Q_RECOV
          case ('current')
            t % score_bins(i) = SCORE_CURRENT
          case ('events')
            t % score_bins(i) = SCORE_EVENTS
          case ('elastic', '(n,elastic)')
            t % score_bins(i) = ELASTIC
          case ('(n,2nd)')
            t % score_bins(i) = N_2ND
          case ('(n,na)')
            t % score_bins(i) = N_2NA
          case ('(n,n3a)')
            t % score_bins(i) = N_N3A
          case ('(n,2na)')
            t % score_bins(i) = N_2NA
          case ('(n,3na)')
            t % score_bins(i) = N_3NA
          case ('(n,np)')
            t % score_bins(i) = N_NP
          case ('(n,n2a)')
            t % score_bins(i) = N_N2A
          case ('(n,2n2a)')
            t % score_bins(i) = N_2N2A
          case ('(n,nd)')
            t % score_bins(i) = N_ND
          case ('(n,nt)')
            t % score_bins(i) = N_NT
          case ('(n,nHe-3)')
            t % score_bins(i) = N_N3HE
          case ('(n,nd2a)')
            t % score_bins(i) = N_ND2A
          case ('(n,nt2a)')
            t % score_bins(i) = N_NT2A
          case ('(n,3nf)')
            t % score_bins(i) = N_3NF
          case ('(n,2np)')
            t % score_bins(i) = N_2NP
          case ('(n,3np)')
            t % score_bins(i) = N_3NP
          case ('(n,n2p)')
            t % score_bins(i) = N_N2P
          case ('(n,npa)')
            t % score_bins(i) = N_NPA
          case ('(n,n1)')
            t % score_bins(i) = N_N1
          case ('(n,nc)')
            t % score_bins(i) = N_NC
          case ('(n,gamma)')
            t % score_bins(i) = N_GAMMA
            depletion_rx = .true.
          case ('(n,p)')
            t % score_bins(i) = N_P
            depletion_rx = .true.
          case ('(n,d)')
            t % score_bins(i) = N_D
          case ('(n,t)')
            t % score_bins(i) = N_T
          case ('(n,3He)')
            t % score_bins(i) = N_3HE
          case ('(n,a)')
            t % score_bins(i) = N_A
            depletion_rx = .true.
          case ('(n,2a)')
            t % score_bins(i) = N_2A
          case ('(n,3a)')
            t % score_bins(i) = N_3A
          case ('(n,2p)')
            t % score_bins(i) = N_2P
          case ('(n,pa)')
            t % score_bins(i) = N_PA
          case ('(n,t2a)')
            t % score_bins(i) = N_T2A
          case ('(n,d2a)')
            t % score_bins(i) = N_D2A
          case ('(n,pd)')
            t % score_bins(i) = N_PD
          case ('(n,pt)')
            t % score_bins(i) = N_PT
          case ('(n,da)')
            t % score_bins(i) = N_DA
          case default
            ! Assume that user has specified an MT number
            MT = int(str_to_int(score_))

            if (MT /= ERROR_INT) then
              ! Specified score was an integer
              if (MT > 1) then
                t % score_bins(i) = MT
              else
                err = E_INVALID_ARGUMENT
                call set_errmsg("Negative MT number cannot be used as a score.")
              end if

            else
              err = E_INVALID_ARGUMENT
              call set_errmsg("Unknown score: " // trim(score_) // ".")
            end if

          end select
        end do

        err = 0
        t % depletion_rx = depletion_rx
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg('Index in tallies array is out of bounds.')
    end if
  end function openmc_tally_set_scores


  subroutine openmc_get_tally_next_id(id) bind(C)
    ! Returns an ID number that has not been used by any other tallies.
    integer(C_INT32_T), intent(out) :: id

    id = largest_tally_id + 1
  end subroutine openmc_get_tally_next_id

end module tally_header
