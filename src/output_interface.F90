module output_interface

  use constants

#ifdef HDF5
  use hdf5_interface
#elif MPI
  use mpi_interface
#endif

  implicit none

#if MPI
  integer :: fh_state_point ! mpi file for state point
  integer :: fh_source      ! mpi file for source
#endif

  interface write_data
    module procedure write_double
    module procedure write_double_1Darray
    module procedure write_integer
    module procedure write_integer_1Darray
    module procedure write_long
    module procedure write_string
  end interface write_data


contains

!===============================================================================
! FILE_CREATE
!===============================================================================

  subroutine file_create(filename, fh_str)

    character(MAX_FILE_LEN) :: filename
    character(MAX_WORD_LEN) :: fh_str

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'state_point') then
      if (master) call hdf5_file_create(trim(trim(filename) // '.h5'), &
                       hdf5_fh)
   else
     call hdf5_parallel_file_create(trim(filename) // '.h5', &
          hdf5_fh)
    endif
# else
   if (trim(fh_str) == 'state_point') then
     call hdf5_file_create(trim(filename) // '.h5', &
          hdf5_fh)
   else
     call hdf5_file_create(trim(filename) // '.h5', &
          hdf5_fh)
   endif
# endif
#elif MPI
   if (trim(fh_str) == 'state_point') then
     call mpi_file_create(trim(filename) // '.binary/' , fh_state_point) 
   else
     call mpi_file_create(trim(filename) // '.binary', fh_source)
   endif
#else
   if (trim(fh_str) == 'state_point') then
     open(UNIT=UNIT_STATE, FILE=trim(filename) // '.binary', STATUS='replace', &
          ACCESS='stream')
   else
     open(UNIT=UNIT_SOURCE, FILE=trim(filename) // '.binary', STATUS='replace', &
          ACCESS='stream')
   endif
#endif

  end subroutine file_create

!===============================================================================
! FILE_CLOSE
!===============================================================================

  subroutine file_close(fh_str)

    character(MAX_WORD_LEN) :: fh_str

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'state_point') then
      if (master) call hdf5_file_close(hdf5_fh)
   else
     call hdf5_parallel_file_close(hdf5_fh)
    endif
# else
     call hdf5_file_close(hdf5_fh)
   endif
# endif
#elif MPI
   if (trim(fh_str) == 'state_point') then
     call mpi_close_file(fh_state_point) 
   else
     call mpi_close_file(fh_source)
   endif
#else
   if (trim(fh_str) == 'state_point') then
     close(UNIT=UNIT_STATE)
   else
     close(UNIT=UNIT_SOURCE)
   endif
#endif

  end subroutine file_close

!===============================================================================
! WRITE_INTEGER
!===============================================================================

  subroutine write_double(buffer, name, group)

    real(8),        intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_double

!===============================================================================
! WRITE_INTEGER_1Darray
!===============================================================================

  subroutine write_double_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    real(8),        intent(in)           :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_double_1Darray

!===============================================================================
! WRITE_INTEGER
!===============================================================================

  subroutine write_integer(buffer, name, group)

    integer,        intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_integer

!===============================================================================
! WRITE_INTEGER_1Darray
!===============================================================================

  subroutine write_integer_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    integer,        intent(in)           :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_integer_1Darray

!===============================================================================
! WRITE_LONG
!===============================================================================

  subroutine write_long(buffer, name, group)

    integer(8),     intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_long(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_long


!===============================================================================
! WRITE_STRING
!===============================================================================

  subroutine write_string(buffer, name, group)

    character(*), intent(in)           :: buffer
    character(*), intent(in)           :: name
    character(*), intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_string(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_string

end module output_interface
