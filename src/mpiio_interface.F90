module mpiio_interface

#ifdef MPI
#ifndef HDF5
  use message_passing

  implicit none

#ifdef MPIF08
#define FH_TYPE type(MPI_File)
#else
#define FH_TYPE integer
#endif

  integer :: mpiio_err ! MPI error code

  ! Generic HDF5 write procedure interface
  interface mpi_write_data
    module procedure mpi_write_double
    module procedure mpi_write_double_1Darray
    module procedure mpi_write_double_2Darray
    module procedure mpi_write_double_3Darray
    module procedure mpi_write_double_4Darray
    module procedure mpi_write_integer
    module procedure mpi_write_integer_1Darray
    module procedure mpi_write_integer_2Darray
    module procedure mpi_write_integer_3Darray
    module procedure mpi_write_integer_4Darray
    module procedure mpi_write_long
    module procedure mpi_write_string
  end interface mpi_write_data

  ! Generic HDF5 read procedure interface
  interface mpi_read_data
    module procedure mpi_read_double
    module procedure mpi_read_double_1Darray
    module procedure mpi_read_double_2Darray
    module procedure mpi_read_double_3Darray
    module procedure mpi_read_double_4Darray
    module procedure mpi_read_integer
    module procedure mpi_read_integer_1Darray
    module procedure mpi_read_integer_2Darray
    module procedure mpi_read_integer_3Darray
    module procedure mpi_read_integer_4Darray
    module procedure mpi_read_long
    module procedure mpi_read_string
  end interface mpi_read_data

contains

!===============================================================================
! MPI_CREATE_FILE creates a file using MPI file I/O
!===============================================================================

  subroutine mpi_create_file(filename, fh)

    character(*), intent(in)    :: filename ! name of file to create
    FH_TYPE,      intent(inout) :: fh       ! file handle

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpiio_err)

  end subroutine mpi_create_file

!===============================================================================
! MPI_OPEN_FILE opens a file using MPI file I/O
!===============================================================================

  subroutine mpi_open_file(filename, fh, mode)

    character(*), intent(in)    :: filename ! name of file to open
    character(*), intent(in)    :: mode     ! open 'r' read, 'w' write
    FH_TYPE,      intent(inout) :: fh       ! file handle

    integer :: open_mode

    ! Determine access mode
    open_mode = MPI_MODE_RDONLY
    if (mode == 'w') then
      open_mode = ior(MPI_MODE_APPEND, MPI_MODE_WRONLY)
    end if

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         open_mode, MPI_INFO_NULL, fh, mpiio_err)

  end subroutine mpi_open_file

!===============================================================================
! MPI_CLOSE_FILE closes a file using MPI file I/O
!===============================================================================

  subroutine mpi_close_file(fh)

    FH_TYPE, intent(inout) :: fh ! file handle

    call MPI_FILE_CLOSE(fh, mpiio_err)

  end subroutine mpi_close_file

!===============================================================================
! MPI_WRITE_INTEGER writes integer scalar data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer(fh, buffer, collect)

    FH_TYPE, intent(in) :: fh      ! file handle
    integer, intent(in) :: buffer  ! data to write
    logical, intent(in) :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_integer

!===============================================================================
! MPI_READ_INTEGER reads integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer(fh, buffer, collect)

    FH_TYPE, intent(in)    :: fh      ! file handle
    integer, intent(inout) :: buffer  ! read data to here
    logical, intent(in)    :: collect ! collective I/O

   if (collect) then
     call MPI_FILE_READ_ALL(fh, buffer, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, mpiio_err)
   else
     call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER, &
          MPI_STATUS_IGNORE, mpiio_err)
   end if

  end subroutine mpi_read_integer

!===============================================================================
! MPI_WRITE_INTEGER_1DARRAY writes integer 1-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer_1Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length    ! length of array
    integer, intent(in) :: buffer(:) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, length, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, length, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_integer_1Darray

!===============================================================================
! MPI_READ_INTEGER_1DARRAY reads integer 1-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer_1Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length    ! length of array
    integer, intent(inout) :: buffer(:) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, length, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, length, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_integer_1Darray

!===============================================================================
! MPI_WRITE_INTEGER_2DARRAY writes integer 2-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer_2Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(2)  ! length of array
    integer, intent(in) :: buffer(length(1),length(2)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_integer_2Darray

!===============================================================================
! MPI_READ_INTEGER_2DARRAY reads integer 2-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer_2Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(2)    ! length of array
    integer, intent(inout) :: buffer(length(1),length(2)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_integer_2Darray

!===============================================================================
! MPI_WRITE_INTEGER_3DARRAY writes integer 3-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer_3Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(3)  ! length of array
    integer, intent(in) :: buffer(length(1),length(2),&
                           length(3)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_integer_3Darray

!===============================================================================
! MPI_READ_INTEGER_3DARRAY reads integer 3-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer_3Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(3)    ! length of array
    integer, intent(inout) :: buffer(length(1),length(2), &
                              length(3)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_integer_3Darray

!===============================================================================
! MPI_WRITE_INTEGER_4DARRAY writes integer 4-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_integer_4Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(4)  ! length of array
    integer, intent(in) :: buffer(length(1),length(2),&
                           length(3),length(4)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_integer_4Darray

!===============================================================================
! MPI_READ_INTEGER_4DARRAY reads integer 4-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_integer_4Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(4)    ! length of array
    integer, intent(inout) :: buffer(length(1),length(2), &
                              length(3),length(4)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_integer_4Darray

!===============================================================================
! MPI_WRITE_DOUBLE writes integer scalar data using MPI File I/O
!===============================================================================

  subroutine mpi_write_double(fh, buffer, collect)

    FH_TYPE, intent(in) :: fh      ! file handle
    real(8), intent(in) :: buffer  ! data to write
    logical, intent(in) :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, 1, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, 1, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_double

!===============================================================================
! MPI_READ_DOUBLE reads integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_double(fh, buffer, collect)

    FH_TYPE, intent(in)    :: fh      ! file handle
    real(8), intent(inout) :: buffer  ! read data to here
    logical, intent(in)    :: collect ! collective I/O

   if (collect) then
     call MPI_FILE_READ_ALL(fh, buffer, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, mpiio_err)
   else
     call MPI_FILE_READ(fh, buffer, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, mpiio_err)
   end if

  end subroutine mpi_read_double

!===============================================================================
! MPI_WRITE_DOUBLE_1DARRAY writes integer 1-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_double_1Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length    ! length of array
    real(8), intent(in) :: buffer(:) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, length, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, length, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_double_1Darray

!===============================================================================
! MPI_READ_DOUBLE_1DARRAY reads integer 1-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_double_1Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length    ! length of array
    real(8), intent(inout) :: buffer(:) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, length, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, length, MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_double_1Darray

!===============================================================================
! MPI_WRITE_DOUBLE_2DARRAY writes integer 2-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_double_2Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(2)  ! length of array
    real(8), intent(in) :: buffer(length(1),length(2)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_double_2Darray

!===============================================================================
! MPI_READ_DOUBLE_2DARRAY reads integer 2-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_double_2Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(2)    ! length of array
    real(8), intent(inout) :: buffer(length(1),length(2)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_double_2Darray

!===============================================================================
! MPI_WRITE_DOUBLE_3DARRAY writes integer 3-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_double_3Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(3)  ! length of array
    real(8), intent(in) :: buffer(length(1),length(2),&
                           length(3)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_double_3Darray

!===============================================================================
! MPI_READ_DOUBLE_3DARRAY reads integer 3-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_double_3Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(3)    ! length of array
    real(8), intent(inout) :: buffer(length(1),length(2), &
                              length(3)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_double_3Darray

!===============================================================================
! MPI_WRITE_DOUBLE_4DARRAY writes integer 4-D array data using MPI File I/O
!===============================================================================

  subroutine mpi_write_double_4Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in) :: fh        ! file handle
    integer, intent(in) :: length(4)  ! length of array
    real(8), intent(in) :: buffer(length(1),length(2),&
                           length(3),length(4)) ! data to write
    logical, intent(in) :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_double_4Darray

!===============================================================================
! MPI_READ_DOUBLE_4DARRAY reads integer 4-D array using MPI file I/O
!===============================================================================

  subroutine mpi_read_double_4Darray(fh, buffer, length, collect)

    FH_TYPE, intent(in)    :: fh        ! file handle
    integer, intent(in)    :: length(4)    ! length of array
    real(8), intent(inout) :: buffer(length(1),length(2), &
                              length(3),length(4)) ! read data to here
    logical, intent(in)    :: collect   ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, product(length), MPI_REAL8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_double_4Darray

!===============================================================================
! MPI_WRITE_LONG writes long integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_write_long(fh, buffer, collect)

    FH_TYPE,    intent(in) :: fh      ! file handle
    integer(8), intent(in) :: buffer  ! data to write
    logical,    intent(in) :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_long

!===============================================================================
! MPI_READ_LONG reads long integer scalar data using MPI file I/O
!===============================================================================

  subroutine mpi_read_long(fh, buffer, collect)

    FH_TYPE,    intent(in)    :: fh      ! file handle
    integer(8), intent(inout) :: buffer  ! read data to here
    logical,    intent(in)    :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_long

!===============================================================================
! MPI_WRITE_STRING writes string data using MPI file I/O
!===============================================================================

  subroutine mpi_write_string(fh, buffer, length, collect)

    character(*), intent(in) :: buffer  ! data to write
    FH_TYPE,      intent(in) :: fh      ! file handle
    integer,      intent(in) :: length  ! length of data
    logical,      intent(in) :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_WRITE_ALL(fh, buffer, length, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_WRITE(fh, buffer, length, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_write_string

!===============================================================================
! MPI_READ_STRING reads string data using MPI file I/O
!===============================================================================

  subroutine mpi_read_string(fh, buffer, length, collect)

    character(*), intent(inout) :: buffer  ! read data to here
    FH_TYPE,      intent(in)    :: fh      ! file handle
    integer,      intent(in)    :: length  ! length of string
    logical,      intent(in)    :: collect ! collective I/O

    if (collect) then
      call MPI_FILE_READ_ALL(fh, buffer, length, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpiio_err)
    else
      call MPI_FILE_READ(fh, buffer, length, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpiio_err)
    end if

  end subroutine mpi_read_string

#endif
#endif
end module mpiio_interface
