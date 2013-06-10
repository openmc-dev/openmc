module ndpp

  use ace_header,       only: Nuclide, GrpTransfer
  use constants
  use global
  use output,           only: write_message

  implicit none

contains

!===============================================================================
! READ_NDPP reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ndpp()

    ! Loop through, ensure that 
    
  end subroutine read_ndpp

!===============================================================================
! GET_INT returns an array of integers read from the current position in the XSS
! array
!===============================================================================

  function get_int(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    integer             :: array(n_values) ! array of values

    array = int(XSS(XSS_index:XSS_index + n_values - 1))
    XSS_index = XSS_index + n_values

  end function get_int

!===============================================================================
! GET_REAL returns an array of real(8)s read from the current position in the
! XSS array
!===============================================================================

  function get_real(n_values) result(array)

    integer, intent(in) :: n_values        ! number of values to read
    real(8)             :: array(n_values) ! array of values

    array = XSS(XSS_index:XSS_index + n_values - 1)
    XSS_index = XSS_index + n_values

  end function get_real

end module ace
