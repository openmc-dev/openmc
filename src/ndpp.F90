module ndpp

  use ace_header,       only: Nuclide, GrpTransfer
  use constants
  use error,            only: fatal_error, warning
  use global
  use output,           only: write_message

  implicit none

contains

!===============================================================================
! READ_NDPP_LIB reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_ndpp_lib()

    ! Loop through each nuclide. Make sure each has correct energy group
    ! structure and was obtained at the same temperature as the nuclide in 
    ! question. Then, read the data per the library type set in tallies.xml 
    ! (currently stored in global % integrated_scatt_lib.
    
  end subroutine read_ndpp_lib


end module ndpp
