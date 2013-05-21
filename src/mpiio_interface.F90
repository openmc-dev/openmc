module mpiio_interface

#ifdef MPI
  use mpi

  implicit none

  integer :: mpi_fh  ! MPI file handle
  integer :: mpi_err ! MPI error code

contains

!===============================================================================
! MPI_CREATE_FILE creates a file using MPI file I/O
!===============================================================================

  subroutine mpi_create_file(filename, fh)

    character(*), intent(in)    :: filename ! name of file to create
    integer,      intent(inout) :: fh       ! file handle

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err) 

  end subroutine mpi_create_file

!===============================================================================
! MPI_OPEN_FILE opens a file using MPI file I/O
!===============================================================================

  subroutine mpi_open_file(filename, fh, mode)

    character(*), intent(in)    :: filename ! name of file to open
    character(*), intent(in)    :: mode     ! open 'r' read, 'w' write
    integer,      intent(inout) :: fh       ! file handle

    integer :: open_mode

    ! Determine access mode
    open_mode = MPI_MODE_RDONLY
    if (mode == 'w') then
      open_mode = MPI_MODE_WRONLY
    end if

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         open_mode, MPI_INFO_NULL, fh, mpi_err) 

  end subroutine mpi_open_file

!===============================================================================
! MPI_CLOSE_FILE closes a file using MPI file I/O
!===============================================================================

  subroutine mpi_close_file(fh)

    integer, intent(inout) :: fh ! file handle

    call MPI_FILE_CLOSE(fh, mpi_err)

  end subroutine mpi_close_file

!===============================================================================
! MPI_WRITE_INTEGER writes integer scalar data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer(fh, buffer)

    integer, intent(in) :: fh     ! file handle
    integer, intent(in) :: buffer ! data to write

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err) 

  end subroutine mpi_write_integer

!===============================================================================
! MPI_WRITE_INTEGER_1DARRAY writes integer 1-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer_1Darray(fh, buffer, length)

    integer, intent(in) :: fh        ! file handle
    integer, intent(in) :: length    ! length of array
    integer, intent(in) :: buffer(:) ! data to write

    call MPI_FILE_WRITE(fh, buffer, length, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_integer_1Darray

!===============================================================================
! MPI_WRITE_LONG writes long integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_write_long(fh, buffer)

    integer,    intent(in) :: fh     ! file handle
    integer(8), intent(in) :: buffer ! data to write

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_long

!===============================================================================
! MPI_WRITE_DOUBLE writes double precision scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_write_double(fh, buffer)

    integer, intent(in) :: fh     ! file handle
    real(8), intent(in) :: buffer ! data to write

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_double

!===============================================================================
! MPI_WRITE_DOUBLE_1DARRAY writes double precision 1-D array using MPI file I/O
!===============================================================================

  subroutine mpi_write_double_1Darray(fh, buffer, length)

    integer, intent(in) :: fh        ! file handle
    integer, intent(in) :: length    ! length of array
    real(8), intent(in) :: buffer(:) ! data to write

    call MPI_FILE_WRITE(fh, buffer, length, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_double_1Darray

!===============================================================================
! MPI_WRITE_STRING writes string data using MPI file I/O
!===============================================================================

  subroutine mpi_write_string(fh, buffer, length)

    character(*), intent(in) :: buffer ! data to write
    integer,      intent(in) :: fh     ! file handle
    integer,      intent(in) :: length ! length of data

    call MPI_FILE_WRITE(fh, buffer, length, MPI_CHARACTER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_string

!===============================================================================
! MPI_READ_INTEGER reads integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer(fh, buffer)

    integer, intent(in)    :: fh     ! file handle
    integer, intent(inout) :: buffer ! read data to here

    call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err) 

  end subroutine mpi_read_integer

!===============================================================================
! MPI_READ_INTEGER_1DARRAY reads integer 1-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer_1Darray(fh, buffer, length)

    integer, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length    ! length of array
    integer, intent(inout) :: buffer(:) ! read data to here

    call MPI_FILE_READ(fh, buffer, length, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_integer_1Darray

!===============================================================================
! MPI_READ_LONG reads long integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_long(fh, buffer)

    integer,    intent(in)    :: fh     ! file handle
    integer(8), intent(inout) :: buffer ! read data to here

    call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_long

!===============================================================================
! MPI_READ_DOUBLE reads double precision scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_double(fh, buffer)

    integer, intent(in)    :: fh     ! file handle
    real(8), intent(inout) :: buffer ! read data to here

    call MPI_FILE_READ(fh, buffer, 1, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_double

!===============================================================================
! MPI_READ_DOUBLE_1DARRAY reads double precision 1-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_double_1Darray(fh, buffer, length)

    integer, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length    ! length of array
    real(8), intent(inout) :: buffer(:) ! read data to here

    call MPI_FILE_READ(fh, buffer, length, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_double_1Darray

!===============================================================================
! MPI_READ_STRING reads string data using MPI file I/O
!===============================================================================

  subroutine mpi_read_string(fh, buffer, length)

    character(*), intent(inout) :: buffer ! read data to here
    integer,      intent(in)    :: fh     ! file handle
    integer,      intent(in)    :: length ! length of string

    call MPI_FILE_READ(fh, buffer, length, MPI_CHARACTER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_string

#endif
end module mpiio_interface
