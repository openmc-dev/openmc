module sab_header

  use, intrinsic :: ISO_C_BINDING

  use dict_header, only: DictCharInt
  use hdf5_interface
  use stl_vector, only: VectorReal
  use string, only: to_c_string

  implicit none
  private

  public :: free_memory_sab

!===============================================================================
! SALPHABETA contains S(a,b) data for thermal neutron scattering, typically off
! of light isotopes such as water, graphite, Be, etc
!===============================================================================

  type, public :: SAlphaBeta
    type(C_PTR) :: ptr
  contains
    procedure :: from_hdf5
    procedure :: calculate_xs
    procedure :: free
    procedure :: has_nuclide
    procedure :: sample
    procedure :: threshold
  end type SAlphaBeta

  ! S(a,b) tables
  type(SAlphaBeta), public, allocatable, target :: sab_tables(:)
  integer(C_INT), public, bind(C) :: n_sab_tables
  type(DictCharInt), public :: sab_dict

  interface
    function sab_from_hdf5(group_id, temperature, n, method, &
         tolerance, minmax) result(ptr) bind(C)
      import HID_T, C_DOUBLE, C_INT, C_PTR
      integer(HID_T), value :: group_id
      real(C_DOUBLE), intent(in) :: temperature
      integer(C_INT), value :: n
      integer(C_INT), value :: method
      real(C_DOUBLE), value :: tolerance
      real(C_DOUBLE), intent(in) :: minmax(2)
      type(C_PTR) :: ptr
    end function

    subroutine sab_calculate_xs(ptr, E, sqrtkT, i_temp, elastic, &
         inelastic) bind(C)
      import C_PTR, C_DOUBLE, C_INT
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      real(C_DOUBLE), value :: sqrtkT
      integer(C_INT), intent(out) :: i_temp
      real(C_DOUBLE), intent(out) :: elastic
      real(C_DOUBLE), intent(out) :: inelastic
    end subroutine

    subroutine sab_free(ptr) bind(C)
      import C_PTR
      type(C_PTR), value :: ptr
    end subroutine

    function sab_has_nuclide(ptr, name) result(val) bind(C)
      import C_PTR, C_CHAR, C_BOOL
      type(C_PTR), value :: ptr
      character(kind=C_CHAR), intent(in) :: name(*)
      logical(C_BOOL) :: val
    end function

    subroutine sab_sample(ptr, micro_xs, E_in, E_out, mu) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: ptr
      type(C_PTR), value :: micro_xs
      real(C_DOUBLE), value :: E_in
      real(C_DOUBLE), intent(out) :: E_out
      real(C_DOUBLE), intent(out) :: mu
    end subroutine

    function sab_threshold(ptr) result(threshold) bind(C)
      import C_PTR, C_double
      type(C_PTR), value :: ptr
      real(C_DOUBLE) :: threshold
    end function
  end interface

contains

  subroutine from_hdf5(this, group_id, temperature, method, &
       tolerance, minmax)
    class(SAlphaBeta), intent(inout) :: this
    integer(HID_T),    intent(in)    :: group_id
    type(VectorReal),  intent(in)    :: temperature ! list of temperatures
    integer(C_INT),    intent(in)    :: method
    real(C_DOUBLE),    intent(in)    :: tolerance
    real(C_DOUBLE),    intent(in)    :: minmax(2)

    real(C_DOUBLE) :: dummy
    integer(C_INT) :: n

    n = temperature % size()
    if (n > 0) then
      this % ptr = sab_from_hdf5(group_id, temperature % data(1), n, method, tolerance, minmax)
    else
      ! In this case, temperatures % data(1) doesn't exist, so we just pass a
      ! dummy value
      this % ptr = sab_from_hdf5(group_id, dummy, n, method, tolerance, minmax)
    end if
  end subroutine from_hdf5

!===============================================================================
! SAB_CALCULATE_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range.
!===============================================================================

  subroutine calculate_xs(this, E, sqrtkT, i_temp, elastic, inelastic)
    class(SAlphaBeta), intent(in) :: this ! S(a,b) object
    real(C_DOUBLE), intent(in) :: E          ! energy
    real(C_DOUBLE), intent(in) :: sqrtkT     ! temperature
    integer, intent(out) :: i_temp    ! index in the S(a,b)'s temperature
    real(C_DOUBLE), intent(out) :: elastic   ! thermal elastic cross section
    real(C_DOUBLE), intent(out) :: inelastic ! thermal inelastic cross section

    call sab_calculate_xs(this % ptr, E, sqrtkT, i_temp, elastic, inelastic)
  end subroutine

  subroutine free(this)
    class(SAlphaBeta), intent(inout) :: this
    call sab_free(this % ptr)
  end subroutine

  function has_nuclide(this, name) result(val)
    class(SAlphaBeta), intent(in) :: this
    character(len=*), intent(in) :: name
    logical(C_BOOL) :: val

    val = sab_has_nuclide(this % ptr, to_c_string(name))
  end function

  subroutine sample(this, micro_xs, E_in, E_out, mu)
    class(SAlphaBeta), intent(in) :: this
    type(C_PTR), value :: micro_xs
    real(C_DOUBLE), value :: E_in
    real(C_DOUBLE), intent(out) :: E_out
    real(C_DOUBLE), intent(out) :: mu

    call sab_sample(this % ptr, micro_xs, E_in, E_out, mu)
  end subroutine

  function threshold(this)
    class(SAlphaBeta), intent(in) :: this
    real(C_DOUBLE) :: threshold
    threshold = sab_threshold(this % ptr)
  end function

!===============================================================================
! FREE_MEMORY_SAB deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_sab()
    integer :: i
    n_sab_tables = 0
    if (allocated(sab_tables)) then
      do i = 1, size(sab_tables)
        call sab_tables(i) % free()
      end do
      deallocate(sab_tables)
    end if
    call sab_dict % clear()
  end subroutine free_memory_sab

end module sab_header
