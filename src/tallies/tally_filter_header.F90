module tally_filter_header

  use, intrinsic :: ISO_C_BINDING

  use constants,       only: MAX_LINE_LEN
  use dict_header,     only: DictIntInt
  use error
  use hdf5_interface,  only: HID_T
  use particle_header, only: Particle
  use stl_vector,      only: VectorInt, VectorReal
  use string,          only: to_str
  use xml_interface,   only: XMLNode

  implicit none
  private
  public :: free_memory_tally_filter
  public :: verify_filter
  public :: openmc_extend_filters
  public :: openmc_filter_get_id
  public :: openmc_filter_set_id
  public :: openmc_get_filter_index
  public :: openmc_get_filter_next_id
  public :: filter_match_pointer

  interface
    function filter_match_pointer(indx) bind(C) result(ptr)
      import C_PTR, C_INT
      integer(C_INT), intent(in), value :: indx
      type(C_PTR)                       :: ptr
    end function filter_match_pointer
  end interface

!===============================================================================
! TALLYFILTERMATCH stores every valid bin and weight for a filter
!===============================================================================

  type, public :: TallyFilterMatch
    type(C_PTR) :: ptr

    ! Index of the bin and weight being used in the current filter combination
    integer          :: i_bin

    ! Indicates whether all valid bins for this filter have been found
    logical          :: bins_present = .false.

  contains
    procedure :: bins_push_back
    procedure :: weights_push_back
    procedure :: bins_clear
    procedure :: weights_clear
    procedure :: bins_size
    procedure :: bins_data
    procedure :: weights_data
    procedure :: bins_set_data
  end type TallyFilterMatch

!===============================================================================
! TALLYFILTER describes a filter that limits what events score to a tally. For
! example, a cell filter indicates that only particles in a specified cell
! should score to the tally.
!===============================================================================

  type, public, abstract :: TallyFilter
    integer :: id
    integer :: n_bins = 0
  contains
    procedure(from_xml_),      deferred :: from_xml
    procedure(get_all_bins_),  deferred :: get_all_bins
    procedure(to_statepoint_), deferred :: to_statepoint
    procedure(text_label_),    deferred :: text_label
    procedure                           :: initialize => filter_initialize
  end type TallyFilter

  abstract interface

    subroutine from_xml_(this, node)
      import TallyFilter, XMLNode
      class(TallyFilter), intent(inout) :: this
      type(XMLNode), intent(in) :: node
    end subroutine from_xml_

!===============================================================================
! GET_NEXT_BIN gives the index for the next valid filter bin and a weight that
! will be applied to the flux.
!
! In principle, a filter can have multiple valid bins.  If current_bin =
! NO_BIN_FOUND, then this method should give the first valid bin.  Providing the
! first valid bin should then give the second valid bin, and so on.  When there
! are no valid bins left, the next_bin should be NO_VALID_BIN.

    subroutine get_all_bins_(this, p, estimator, match)
      import TallyFilter
      import Particle
      import TallyFilterMatch
      class(TallyFilter), intent(in)  :: this
      type(Particle),     intent(in)  :: p
      integer,            intent(in)  :: estimator
      type(TallyFilterMatch), intent(inout) :: match
    end subroutine get_all_bins_

!===============================================================================
! TO_STATEPOINT writes all the information needed to reconstruct the filter to
! the given filter_group.

    subroutine to_statepoint_(this, filter_group)
      import TallyFilter
      import HID_T
      class(TallyFilter), intent(in) :: this
      integer(HID_T),     intent(in) :: filter_group
    end subroutine to_statepoint_

!===============================================================================
! TEXT_LABEL returns a string describing the given filter bin.  For example, an
! energy filter might return the string "Incoming Energy [0.625E-6, 20.0)".
! This is used to write the tallies.out file.

    function text_label_(this, bin) result(label)
      import TallyFilter
      import MAX_LINE_LEN
      class(TallyFilter), intent(in) :: this
      integer,            intent(in) :: bin
      character(MAX_LINE_LEN)        :: label
    end function text_label_

  end interface

!===============================================================================
!===============================================================================

  type, public, abstract, extends(TallyFilter) :: CppTallyFilter
    type(C_PTR) :: ptr
  contains
    procedure :: from_xml_c
    procedure :: get_all_bins_c
    procedure :: to_statepoint_c
    procedure :: initialize_c
  end type CppTallyFilter

!===============================================================================
! TALLYFILTERCONTAINER contains an allocatable TallyFilter object for arrays of
! TallyFilters
!===============================================================================

  type, public :: TallyFilterContainer
    class(TallyFilter), allocatable :: obj
  end type TallyFilterContainer

  integer(C_INT32_T), public, bind(C) :: n_filters = 0 ! # of filters

  type(TallyFilterContainer), public, allocatable, target :: filters(:)
  type(TallyFilterMatch), public, allocatable :: filter_matches(:)
!$omp threadprivate(filter_matches)

  ! Dictionary that maps user IDs to indices in 'filters'
  type(DictIntInt), public :: filter_dict

  ! The largest filter ID that has been specified in the system.  This is useful
  ! in case the code needs to find an ID for a new filter.
  integer :: largest_filter_id

contains

!===============================================================================
!===============================================================================

  subroutine bins_push_back(this, val)
    class(TallyFilterMatch), intent(inout) :: this
    integer,                 intent(in)    :: val
    interface
      subroutine filter_match_bins_push_back(ptr, val) bind(C)
        import C_PTR, C_INT
        type(C_PTR),                value :: ptr
        integer(C_INT), intent(in), value :: val
      end subroutine
    end interface
    call filter_match_bins_push_back(this % ptr, val)
  end subroutine bins_push_back

  subroutine weights_push_back(this, val)
    class(TallyFilterMatch), intent(inout) :: this
    real(8),                 intent(in)    :: val
    interface
      subroutine filter_match_weights_push_back(ptr, val) bind(C)
        import C_PTR, C_DOUBLE
        type(C_PTR),                value :: ptr
        real(C_DOUBLE), intent(in), value :: val
      end subroutine
    end interface
    call filter_match_weights_push_back(this % ptr, val)
  end subroutine weights_push_back

  subroutine bins_clear(this)
    class(TallyFilterMatch), intent(inout) :: this
    interface
      subroutine filter_match_bins_clear(ptr) bind(C)
        import C_PTR
        type(C_PTR), value :: ptr
      end subroutine
    end interface
    call filter_match_bins_clear(this % ptr)
  end subroutine bins_clear

  subroutine weights_clear(this)
    class(TallyFilterMatch), intent(inout) :: this
    interface
      subroutine filter_match_weights_clear(ptr) bind(C)
        import C_PTR
        type(C_PTR), value :: ptr
      end subroutine
    end interface
    call filter_match_weights_clear(this % ptr)
  end subroutine weights_clear

  function bins_size(this) result(len)
    class(TallyFilterMatch), intent(inout) :: this
    integer                                :: len
    interface
      function filter_match_bins_size(ptr) bind(C) result(len)
        import C_PTR, C_INT
        type(C_PTR), value :: ptr
        integer(C_INT)     :: len
      end function
    end interface
    len = filter_match_bins_size(this % ptr)
  end function bins_size

  function bins_data(this, indx) result(val)
    class(TallyFilterMatch), intent(inout) :: this
    integer,                 intent(in)    :: indx
    integer                                :: val
    interface
      function filter_match_bins_data(ptr, indx) bind(C) result(val)
        import C_PTR, C_INT
        type(C_PTR),                value :: ptr
        integer(C_INT), intent(in), value :: indx
        integer(C_INT)                    :: val
      end function
    end interface
    val = filter_match_bins_data(this % ptr, indx)
  end function bins_data

  function weights_data(this, indx) result(val)
    class(TallyFilterMatch), intent(inout) :: this
    integer,                 intent(in)    :: indx
    real(8)                                :: val
    interface
      function filter_match_weights_data(ptr, indx) bind(C) result(val)
        import C_PTR, C_INT, C_DOUBLE
        type(C_PTR),                value :: ptr
        integer(C_INT), intent(in), value :: indx
        real(C_DOUBLE)                    :: val
      end function
    end interface
    val = filter_match_weights_data(this % ptr, indx)
  end function weights_data

  subroutine bins_set_data(this, indx, val)
    class(TallyFilterMatch), intent(inout) :: this
    integer,                 intent(in)    :: indx
    integer,                 intent(in)    :: val
    interface
      subroutine filter_match_bins_set_data(ptr, indx, val) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value             :: ptr
        integer(C_INT), value, intent(in) :: indx
        integer(C_INT), value, intent(in) :: val
      end subroutine
    end interface
    call filter_match_bins_set_data(this % ptr, indx, val)
  end subroutine bins_set_data

!===============================================================================
! INITIALIZE sets up any internal data, as necessary.  If this procedure is not
! overriden by the derived class, then it will do nothing by default.

  subroutine filter_initialize(this)
    class(TallyFilter), intent(inout) :: this
  end subroutine filter_initialize

!===============================================================================

  subroutine from_xml_c(this, node)
    class(CppTallyFilter), intent(inout) :: this
    class(XMLNode),        intent(in)    :: node
    interface
      subroutine filter_from_xml(filt, node) bind(C)
        import C_PTR
        type(C_PTR), value :: filt
        type(C_PTR) :: node
      end subroutine filter_from_xml
    end interface
    call filter_from_xml(this % ptr, node % ptr)
  end subroutine from_xml_c

  subroutine get_all_bins_c(this, p, estimator, match)
    class(CppTallyFilter), intent(in) :: this
    type(Particle),     intent(in)    :: p
    integer,            intent(in)    :: estimator
    type(TallyFilterMatch), intent(inout) :: match
    interface
      subroutine filter_get_all_bins(filt, p, estimator, match) bind(C)
        import C_PTR, Particle, C_INT
        type(C_PTR),                value :: filt
        type(Particle), intent(in)        :: p
        integer(C_INT), intent(in), value :: estimator
        type(C_PTR),                value :: match
      end subroutine filter_get_all_bins
    end interface
    call filter_get_all_bins(this % ptr, p, estimator, match % ptr)
  end subroutine get_all_bins_c

  subroutine to_statepoint_c(this, filter_group)
    class(CppTallyFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group
    interface
      subroutine filter_to_statepoint(filt, filter_group) bind(C)
        import C_PTR, HID_T
        type(C_PTR),                value :: filt
        integer(HID_T), intent(in), value :: filter_group
      end subroutine filter_to_statepoint
    end interface
    call filter_to_statepoint(this % ptr, filter_group)
  end subroutine to_statepoint_c

  subroutine initialize_c(this)
    class(CppTallyFilter), intent(inout) :: this
    interface
      subroutine filter_initialize(filt) bind(C)
        import C_PTR
        type(C_PTR), value :: filt
      end subroutine filter_initialize
    end interface
    call filter_initialize(this % ptr)
  end subroutine initialize_c

!===============================================================================
! FREE_MEMORY_TALLY_FILTER deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_tally_filter()
    n_filters = 0
    if (allocated(filters)) deallocate(filters)
    call filter_dict % clear()
    largest_filter_id = 0
  end subroutine free_memory_tally_filter

!===============================================================================
! VERIFY_FILTER makes sure that given a filter index, the size of the filters
! array is sufficient and a filter object has already been allocated.
!===============================================================================

  function verify_filter(index) result(err)
    integer(C_INT32_T), intent(in) :: index
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (.not. allocated(filters(index) % obj)) then
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function verify_filter

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_filters(n, index_start, index_end) result(err) bind(C)
    ! Creates or extends the filters array
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    integer :: i ! loop counter
    type(TallyFilterContainer), allocatable :: temp(:) ! temporary filters

    if (n_filters == 0) then
      ! Allocate filters array
      allocate(filters(n))
    else
      ! Move filters to temporary array
      allocate(temp(n_filters + n))
      do i = 1, n_filters
        call move_alloc(filters(i) % obj, temp(i) % obj)
      end do

      ! Move filters back from temporary array to filters array
      call move_alloc(temp, filters)
    end if

    ! Return indices in filters array
    if (present(index_start)) index_start = n_filters + 1
    if (present(index_end)) index_end = n_filters + n
    n_filters = n_filters + n

    err = 0
  end function openmc_extend_filters


  function openmc_filter_get_id(index, id) result(err) bind(C)
    ! Return the ID of a filter
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_filters) then
      id = filters(index) % obj % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_filter_get_id


  function openmc_filter_set_id(index, id) result(err) bind(C)
    ! Set the ID of a filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        filters(index) % obj % id = id
        call filter_dict % set(id, index)
        if (id > largest_filter_id) largest_filter_id = id

        err = 0
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_filter_set_id


  function openmc_get_filter_index(id, index) result(err) bind(C)
    ! Returns the index in the filters array of a filter with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(filters)) then
      if (filter_dict % has(id)) then
        index = filter_dict % get(id)
        err = 0
      else
        err = E_INVALID_ID
        call set_errmsg("No filter exists with ID=" // trim(to_str(id)) // ".")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory has not been allocated for filters.")
    end if
  end function openmc_get_filter_index


  subroutine openmc_get_filter_next_id(id) bind(C)
    ! Returns an ID number that has not been used by any other filters.
    integer(C_INT32_T), intent(out) :: id

    id = largest_filter_id + 1
  end subroutine openmc_get_filter_next_id

end module tally_filter_header
