module mpi_interface

  use constants

#ifdef MPI
  use mpi
#endif

  implicit none

  ! define array writing interface
! interface mpi_write_data
!   module procedure mpi_write_integer
!   module procedure mpi_write_double
!   module procedure mpi_write_double_1Darray
!   module procedure mpi_write_integer_1Darray
! end interface mpi_write_data

  integer :: mpi_err

contains

!===============================================================================
! MPI_FILE_CREATE
!===============================================================================

  subroutine mpi_file_create(filename, file_id)

    character(MAX_FILE_LEN) :: filename
    integer                 :: file_id

    ! Create the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, file_id, mpi_err) 

  end subroutine mpi_file_create

!===============================================================================
! MPI_CLOSE_FILE
!===============================================================================

  subroutine mpi_close_file(file_id)

    integer :: file_id

    call MPI_FILE_CLOSE(file_id, mpi_err)

  end subroutine mpi_close_file

end module mpi_interface
