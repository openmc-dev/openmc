module mpi_interface

  use constants
  use global

#ifdef MPI
  use mpi


  implicit none

! integer :: mpi_fh

contains

!===============================================================================
! MPI_CREATE_FILE
!===============================================================================

  subroutine mpi_create_file(filename, fh)

    character(*) :: filename
    integer, intent(in)     :: fh

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err) 

  end subroutine mpi_create_file

!===============================================================================
! MPI_OPEN_FILE
!===============================================================================

  subroutine mpi_open_file(filename, fh, mode)

    character(*), intent(in) :: filename
    character(*), intent(in) :: mode
    integer,      intent(in) :: fh

    integer :: open_mode

    ! Determine access mode
    open_mode = MPI_MODE_RDONLY
    if (mode == 'rw') then
      open_mode = MPI_MODE_RDWR
    end if

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         open_mode, MPI_INFO_NULL, fh, mpi_err) 

  end subroutine mpi_open_file

!===============================================================================
! MPI_CLOSE_FILE
!===============================================================================

  subroutine mpi_close_file(fh)

    integer, intent(in) :: fh

    call MPI_FILE_CLOSE(fh, mpi_err)

  end subroutine mpi_close_file

!===============================================================================
! MPI_WRITE_INTEGER
!===============================================================================

  subroutine mpi_write_integer(fh, buffer)

    integer, intent(in) :: fh
    integer, intent(in) :: buffer

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err) 

  end subroutine mpi_write_integer

!===============================================================================
! MPI_WRITE_INTEGER_1DARRAY
!===============================================================================

  subroutine mpi_write_integer_1Darray(fh, buffer, length)

    integer, intent(in) :: fh
    integer, intent(in) :: length
    integer, intent(in) :: buffer(:)

    call MPI_FILE_WRITE(fh, buffer, length, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_integer_1Darray

!===============================================================================
! MPI_WRITE_LONG
!===============================================================================

  subroutine mpi_write_long(fh, buffer)

    integer,    intent(in) :: fh
    integer(8), intent(in) :: buffer

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_long

!===============================================================================
! MPI_WRITE_DOUBLE
!===============================================================================

  subroutine mpi_write_double(fh, buffer)

    integer, intent(in) :: fh
    real(8), intent(in) :: buffer

    call MPI_FILE_WRITE(fh, buffer, 1, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_double

!===============================================================================
! MPI_WRITE_DOUBLE_1DARRAY
!===============================================================================

  subroutine mpi_write_double_1Darray(fh, buffer, length)

    integer, intent(in) :: fh
    integer, intent(in) :: length
    real(8), intent(in) :: buffer(:)

    call MPI_FILE_WRITE(fh, buffer, length, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_double_1Darray

!===============================================================================
! MPI_WRITE_STRING
!===============================================================================

  subroutine mpi_write_string(fh, buffer, length)

    character(*), intent(in) :: buffer
    integer,      intent(in) :: fh
    integer,      intent(in) :: length

    call MPI_FILE_WRITE(fh, buffer, length, MPI_CHARACTER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_write_string

!===============================================================================
! MPI_READ_INTEGER
!===============================================================================

  subroutine mpi_read_integer(fh, buffer)

    integer, intent(in) :: fh
    integer, intent(inout) :: buffer

    call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err) 

  end subroutine mpi_read_integer

!===============================================================================
! MPI_READ_INTEGER_1DARRAY
!===============================================================================

  subroutine mpi_read_integer_1Darray(fh, buffer, length)

    integer, intent(in) :: fh
    integer, intent(in) :: length
    integer, intent(inout) :: buffer(:)

    call MPI_FILE_READ(fh, buffer, length, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_integer_1Darray

!===============================================================================
! MPI_READ_LONG
!===============================================================================

  subroutine mpi_read_long(fh, buffer)

    integer,    intent(in) :: fh
    integer(8), intent(inout) :: buffer

    call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_long

!===============================================================================
! MPI_READ_DOUBLE
!===============================================================================

  subroutine mpi_read_double(fh, buffer)

    integer, intent(in) :: fh
    real(8), intent(inout) :: buffer

    call MPI_FILE_READ(fh, buffer, 1, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_double

!===============================================================================
! MPI_READ_DOUBLE_1DARRAY
!===============================================================================

  subroutine mpi_read_double_1Darray(fh, buffer, length)

    integer, intent(in) :: fh
    integer, intent(in) :: length
    real(8), intent(inout) :: buffer(:)

    call MPI_FILE_READ(fh, buffer, length, MPI_REAL8, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_double_1Darray

!===============================================================================
! MPI_READ_STRING
!===============================================================================

  subroutine mpi_read_string(fh, buffer, length)

    character(*), intent(inout) :: buffer
    integer,      intent(in) :: fh
    integer,      intent(in) :: length

    call MPI_FILE_READ(fh, buffer, length, MPI_CHARACTER, &
         MPI_STATUS_IGNORE, mpi_err)

  end subroutine mpi_read_string

#endif
end module mpi_interface
