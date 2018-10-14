module tally_filter_sph_harm

  use, intrinsic :: ISO_C_BINDING

  use error
  use string,              only: to_str, to_lower, to_f_string
  use tally_filter_header

  implicit none
  private
  public :: openmc_sphharm_filter_get_order
  public :: openmc_sphharm_filter_get_cosine
  public :: openmc_sphharm_filter_set_order
  public :: openmc_sphharm_filter_set_cosine

  integer, public, parameter :: COSINE_SCATTER = 1
  integer, public, parameter :: COSINE_PARTICLE = 2

!===============================================================================
! SPHERICALHARMONICSFILTER gives spherical harmonics expansion moments of a
! tally score
!===============================================================================

  type, public, extends(CppTallyFilter) :: SphericalHarmonicsFilter
  contains
    procedure :: cosine
  end type SphericalHarmonicsFilter

contains

  function cosine(this) result(val)
    class(SphericalHarmonicsFilter) :: this
    integer                         :: val
    interface
      function sphharm_filter_get_cosine(filt) result(val) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT) :: val
      end function
    end interface
    val = sphharm_filter_get_cosine(this % ptr)
  end function cosine

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_sphharm_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    interface
      function sphharm_filter_get_order(filt) result(order) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: order
      end function sphharm_filter_get_order
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SphericalHarmonicsFilter)
        order = sphharm_filter_get_order(f % ptr)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spherical harmonics filter.")
      end select
    end if
  end function openmc_sphharm_filter_get_order


  function openmc_sphharm_filter_get_cosine(index, cosine) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value :: index
    character(kind=C_CHAR), intent(out) :: cosine(*)
    integer(C_INT) :: err

    integer :: i
    character(10) :: cosine_

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SphericalHarmonicsFilter)
        select case (f % cosine())
        case (COSINE_SCATTER)
          cosine_ = 'scatter'
        case (COSINE_PARTICLE)
          cosine_ = 'particle'
        end select

        ! Convert to C string
        do i = 1, len_trim(cosine_)
          cosine(i) = cosine_(i:i)
        end do
        cosine(len_trim(cosine_) + 1) = C_NULL_CHAR

      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spherical harmonics filter.")
      end select
    end if
  end function openmc_sphharm_filter_get_cosine


  function openmc_sphharm_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    interface
      subroutine sphharm_filter_set_order(filt, order) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: order
      end subroutine
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SphericalHarmonicsFilter)
        call sphharm_filter_set_order(f % ptr, order)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spherical harmonics filter.")
      end select
    end if
  end function openmc_sphharm_filter_set_order


  function openmc_sphharm_filter_set_cosine(index, cosine) result(err) bind(C)
    ! Set the cosine parameter
    integer(C_INT32_T), value :: index
    character(kind=C_CHAR), intent(in) :: cosine(*)
    integer(C_INT) :: err

    character(:), allocatable :: cosine_

    interface
      subroutine sphharm_filter_set_cosine(filt, cosine) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: cosine
      end subroutine
    end interface

    ! Convert C string to Fortran string
    cosine_ = to_f_string(cosine)

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SphericalHarmonicsFilter)
        select case (cosine_)
        case ('scatter')
          call sphharm_filter_set_cosine(f % ptr, COSINE_SCATTER)
        case ('particle')
          call sphharm_filter_set_cosine(f % ptr, COSINE_PARTICLE)
        end select

      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spherical harmonics filter.")
      end select
    end if
  end function openmc_sphharm_filter_set_cosine

end module tally_filter_sph_harm
