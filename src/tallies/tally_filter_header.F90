module tally_filter_header

  use, intrinsic :: ISO_C_BINDING

  use constants,       only: MAX_LINE_LEN
  use dict_header,     only: DictIntInt
  use error
  use particle_header, only: Particle
  use stl_vector,      only: VectorInt, VectorReal
  use string,          only: to_str
  use xml_interface,   only: XMLNode

  use hdf5

  implicit none
  private
  public :: free_memory_tally_filter
  public :: openmc_extend_filters
  public :: openmc_filter_get_id
  public :: openmc_filter_set_id
  public :: openmc_get_filter_index
  public :: openmc_get_filter_next_id

!===============================================================================
! TALLYFILTERMATCH stores every valid bin and weight for a filter
!===============================================================================

  type, public :: TallyFilterMatch
    ! Index of the bin and weight being used in the current filter combination
    integer          :: i_bin
    type(VectorInt)  :: bins
    type(VectorReal) :: weights

    ! Indicates whether all valid bins for this filter have been found
    logical          :: bins_present = .false.
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
! INITIALIZE sets up any internal data, as necessary.  If this procedure is not
! overriden by the derived class, then it will do nothing by default.

  subroutine filter_initialize(this)
    class(TallyFilter), intent(inout) :: this
  end subroutine filter_initialize

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
