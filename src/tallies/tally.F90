module tally

  use, intrinsic :: ISO_C_BINDING

  use bank_header
  use constants
  use material_header
  use message_passing
  use settings
  use simulation_header
  use tally_derivative_header
  use tally_filter
  use tally_header

  implicit none

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_tally_allocate(index, type) result(err) bind(C)
    ! Set the type of the tally
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    integer(C_INT32_T) :: empty(0)
    character(:), allocatable :: type_

    interface
      function tally_pointer(indx) bind(C) result(ptr)
        import C_INT, C_PTR
        integer(C_INT), value :: indx
        type(C_PTR)           :: ptr
      end function
    end interface

    ! Convert C string to Fortran string
    type_ = to_f_string(type)

    err = 0
    if (index >= 1 .and. index <= n_tallies) then
      if (allocated(tallies(index) % obj)) then
        err = E_ALLOCATE
        call set_errmsg("Tally type has already been set.")
      else
        select case (type_)
        case ('generic')
          allocate(TallyObject :: tallies(index) % obj)
        case default
          err = E_UNASSIGNED
          call set_errmsg("Unknown tally type: " // trim(type_))
        end select

        ! Assign the pointer to the C++ tally
        tallies(index) % obj % ptr = tally_pointer(index - 1)

        ! When a tally is allocated, set it to have 0 filters
        err = openmc_tally_set_filters(index, 0, empty)
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in tallies array is out of bounds.")
    end if
  end function openmc_tally_allocate

end module tally
