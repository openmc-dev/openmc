module tally_filter_cpp
  use, intrinsic :: ISO_C_BINDING

  use constants,       only: MAX_LINE_LEN
  use hdf5_interface,  only: HID_T
  use particle_header, only: Particle
  use tally_filter_header
  use xml_interface,   only: XMLNode

  implicit none

!===============================================================================
! CPPTALLYFILTER is a TallyFilter that is at least partially implemented in C++.
!
! Moving a tally to C++ is easier if done piece-by-piece, so this filter has
! *_cpp_inner procedures that allows the Fortran code to override e.g.
! get_all_bins and still call the C++ implementation for side-by-side comparison
!===============================================================================

  type, abstract, extends(TallyFilter) :: CppTallyFilter
    type(C_PTR) :: ptr
  contains
    procedure :: n_bins_cpp
    procedure :: from_xml_cpp_inner
    procedure :: initialize_cpp_inner
    procedure :: from_xml => from_xml_cpp_default
    procedure :: get_all_bins => get_all_bins_cpp
    procedure :: to_statepoint => to_statepoint_cpp
    procedure :: text_label => text_label_cpp
    procedure :: initialize => initialize_cpp_default
  end type CppTallyFilter

!===============================================================================
! Pure C++ filters
!===============================================================================

  type, extends(CppTallyFilter) :: AzimuthalFilter
  end type

  type, extends(CppTallyFilter) :: CellFilter
  end type

  type, extends(CppTallyFilter) :: CellbornFilter
  end type

  type, extends(CellFilter) :: CellFromFilter
  end type

  type, extends(CppTallyFilter) :: EnergyFunctionFilter
  end type

  type, extends(CppTallyFilter) :: MaterialFilter
  end type

  type, extends(CppTallyFilter) :: MuFilter
  end type

  type, extends(CppTallyFilter) :: PolarFilter
  end type

  type, extends(CppTallyFilter) :: SpatialLegendreFilter
  end type

  type, extends(CppTallyFilter) :: SurfaceFilter
    ! True if this filter is used for surface currents
    logical :: current = .false.
  end type

  type, extends(CppTallyFilter) :: UniverseFilter
  end type

  type, extends(CppTallyFilter) :: ZernikeFilter
  end type

  type, extends(ZernikeFilter) :: ZernikeRadialFilter
  end type

contains

!===============================================================================
! CppTallyFilter implementation
!===============================================================================

  function n_bins_cpp(this) result(n_bins)
    class(CppTallyFilter), intent(inout) :: this
    integer                              :: n_bins
    interface
      function filter_n_bins(filt) result(n_bins) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: n_bins
      end function filter_n_bins
    end interface
    n_bins = filter_n_bins(this % ptr)
  end function n_bins_cpp

  subroutine from_xml_cpp_inner(this, node)
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
  end subroutine from_xml_cpp_inner

  subroutine initialize_cpp_inner(this)
    class(CppTallyFilter), intent(inout) :: this
    interface
      subroutine filter_initialize(filt) bind(C)
        import C_PTR
        type(C_PTR), value :: filt
      end subroutine filter_initialize
    end interface
    call filter_initialize(this % ptr)
  end subroutine initialize_cpp_inner

  subroutine from_xml_cpp_default(this, node)
    class(CppTallyFilter), intent(inout) :: this
    type(XMLNode),         intent(in)    :: node
    call this % from_xml_cpp_inner(node)
    this % n_bins = this % n_bins_cpp()
  end subroutine from_xml_cpp_default

  subroutine get_all_bins_cpp(this, p, estimator, match)
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
  end subroutine get_all_bins_cpp

  subroutine to_statepoint_cpp(this, filter_group)
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
  end subroutine to_statepoint_cpp

  function text_label_cpp(this, bin) result(label)
    class(CppTallyFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label
    character(kind=C_CHAR)            :: label_(MAX_LINE_LEN+1)
    integer :: i
    interface
      subroutine filter_text_label(filt, bin, label) bind(C)
        import C_PTR, C_INT, C_CHAR
        type(C_PTR), value     :: filt
        integer(C_INT), value  :: bin
        character(kind=C_CHAR) :: label(*)
      end subroutine filter_text_label
    end interface
    call filter_text_label(this % ptr, bin, label_)
    label = " "
    do i = 1, MAX_LINE_LEN
      if (label_(i) == C_NULL_CHAR) exit
      label(i:i) = label_(i)
    end do
  end function text_label_cpp

  subroutine initialize_cpp_default(this)
    class(CppTallyFilter), intent(inout) :: this
    call this % initialize_cpp_inner()
  end subroutine initialize_cpp_default

!===============================================================================
! FILTER_FROM_F given a Fortran index, return a pointer to a C++ filter.
!===============================================================================

  function filter_from_f(index) result(filt) bind(C)
    integer(C_INT32_T), intent(in), value :: index
    type(C_PTR) :: filt

    filt = C_NULL_PTR
    select type(f => filters(index) % obj)
    class is (CppTallyFilter)
      filt = f % ptr
    end select
  end function

!===============================================================================
! FILTER_UPDATE_N_BINS given a Fortran index, updates filt % n_bins using C++.
!===============================================================================

  subroutine filter_update_n_bins(index) bind(C)
    integer(C_INT32_T), intent(in), value :: index
    select type(f => filters(index) % obj)
    class is (CppTallyFilter)
      f % n_bins = f % n_bins_cpp()
    end select
  end subroutine
end module tally_filter_cpp
