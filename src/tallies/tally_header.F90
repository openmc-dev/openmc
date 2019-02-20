module tally_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use dict_header,         only: DictIntInt
  use hdf5_interface,      only: HID_T, HSIZE_T
  use message_passing,     only: n_procs, master
  use settings,            only: reduce_tallies, run_mode
  use stl_vector,          only: VectorInt
  use string,              only: to_lower, to_f_string, str_to_int, to_str, to_c_string
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

    integer :: total_score_bins

  contains
    procedure :: id => tally_get_id
    procedure :: set_id => tally_set_id
    procedure :: type => tally_get_type
    procedure :: set_type => tally_set_type
    procedure :: estimator => tally_get_estimator
    procedure :: set_estimator => tally_set_estimator
    procedure :: depletion_rx => tally_get_depletion_rx
    procedure :: n_score_bins => tally_get_n_score_bins
    procedure :: score_bins => tally_get_score_bin
    procedure :: n_filters => tally_get_n_filters
    procedure :: filter => tally_get_filter
    procedure :: n_filter_bins => tally_get_n_filter_bins
    procedure :: n_nuclide_bins => tally_get_n_nuclide_bins
    procedure :: nuclide_bins => tally_get_nuclide_bins
    procedure :: energyout_filter => tally_get_energyout_filter
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

contains

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

  function tally_get_depletion_rx(this) result(drx)
    class(TallyObject) :: this
    logical(C_BOOL) :: drx
    interface
      function tally_get_depletion_rx_c(tally) result(drx) bind(C)
        import C_PTR, C_BOOL
        type(C_PTR), value :: tally
        logical(C_BOOl) :: drx
      end function
    end interface
    drx = tally_get_depletion_rx_c(this % ptr)
  end function

  function tally_get_n_score_bins(this) result(n)
    class(TallyObject) :: this
    integer(C_INT) :: n
    interface
      function tally_get_n_scores_c(tally) result(n) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT) :: n
      end function
    end interface
    n = tally_get_n_scores_c(this % ptr)
  end function

  function tally_get_score_bin(this, i) result(filt)
    class(TallyObject) :: this
    integer(C_INT) :: i
    integer(C_INT32_T) :: filt
    interface
      function tally_get_score_c(tally, i) result(filt) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally
        integer(C_INT), value :: i
        integer(C_INT) :: filt
      end function
    end interface
    filt = tally_get_score_c(this % ptr, i-1)
  end function

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


  subroutine openmc_get_tally_next_id(id) bind(C)
    ! Returns an ID number that has not been used by any other tallies.
    integer(C_INT32_T), intent(out) :: id

    id = largest_tally_id + 1
  end subroutine openmc_get_tally_next_id

end module tally_header
