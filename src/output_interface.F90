module output_interface

  use constants
  use error,         only: warning, fatal_error
  use global
  use tally_header,  only: TallyResult

#ifdef HDF5
  use hdf5_interface
#endif
#ifdef MPI
  use mpiio_interface
#endif

  implicit none
  private

  type, public :: BinaryOutput
    private
    ! Compilation specific data
#ifdef HDF5
    integer(HID_T) :: hdf5_fh
    integer(HID_T) :: hdf5_grp
#else
    integer :: unit_fh
# endif
    logical :: serial ! Serial I/O when using MPI/PHDF5
   contains
    generic, public :: write_data =>  write_double, &
                                      write_double_1Darray, &
                                      write_double_2Darray, &
                                      write_double_3Darray, &
                                      write_double_4Darray, &
                                      write_integer, &
                                      write_integer_1Darray, &
                                      write_integer_2Darray, &
                                      write_integer_3Darray, &
                                      write_integer_4Darray, &
                                      write_long, &
                                      write_string
    generic, public :: read_data => read_double, &
                                    read_double_1Darray, &
                                    read_double_2Darray, &
                                    read_double_3Darray, &
                                    read_double_4Darray, &
                                    read_integer, &
                                    read_integer_1Darray, &
                                    read_integer_2Darray, &
                                    read_integer_3Darray, &
                                    read_integer_4Darray, &
                                    read_long, &
                                    read_string
    procedure :: write_double => write_double
    procedure :: write_double_1Darray => write_double_1Darray
    procedure :: write_double_2Darray => write_double_2Darray
    procedure :: write_double_3Darray => write_double_3Darray
    procedure :: write_double_4Darray => write_double_4Darray
    procedure :: write_integer => write_integer
    procedure :: write_integer_1Darray => write_integer_1Darray
    procedure :: write_integer_2Darray => write_integer_2Darray
    procedure :: write_integer_3Darray => write_integer_3Darray
    procedure :: write_integer_4Darray => write_integer_4Darray
    procedure :: write_long => write_long
    procedure :: write_string => write_string
    procedure :: read_double => read_double
    procedure :: read_double_1Darray => read_double_1Darray
    procedure :: read_double_2Darray => read_double_2Darray
    procedure :: read_double_3Darray => read_double_3Darray
    procedure :: read_double_4Darray => read_double_4Darray
    procedure :: read_integer => read_integer
    procedure :: read_integer_1Darray => read_integer_1Darray
    procedure :: read_integer_2Darray => read_integer_2Darray
    procedure :: read_integer_3Darray => read_integer_3Darray
    procedure :: read_integer_4Darray => read_integer_4Darray
    procedure :: read_long => read_long
    procedure :: read_string => read_string
    procedure, public :: file_create => file_create
    procedure, public :: file_open => file_open
    procedure, public :: file_close => file_close
    procedure, public :: write_tally_result => write_tally_result
    procedure, public :: read_tally_result => read_tally_result
    procedure, public :: write_source_bank => write_source_bank
    procedure, public :: read_source_bank => read_source_bank
#ifdef HDF5
    procedure, public :: write_attribute_string => write_attribute_string
    procedure, public :: open_group => open_group
    procedure, public :: close_group => close_group
#endif
  end type BinaryOutput

contains

!===============================================================================
! FILE_CREATE creates a new file to write data to
!===============================================================================

  subroutine file_create(self, filename, serial)

    character(*),      intent(in) :: filename  ! name of file to be created
    logical, optional, intent(in) :: serial    ! processor rank to write from
    class(BinaryOutput) :: self

    ! Check for serial option
    if (present(serial)) then
      self % serial = serial
    else
      self % serial = .true.
     end if

#ifdef HDF5
# ifdef MPI
    if (self % serial) then
      call hdf5_file_create(filename, self % hdf5_fh)
    else
      call hdf5_file_create_parallel(filename, self % hdf5_fh)
    endif
# else
    call hdf5_file_create(filename, self % hdf5_fh)
# endif
#elif MPI
    if (self % serial) then
      open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
           STATUS='replace', ACCESS='stream')
    else
      call mpi_create_file(filename, self % unit_fh)
    end if
#else
    open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
         STATUS='replace', ACCESS='stream')
#endif

  end subroutine file_create

!===============================================================================
! FILE_OPEN opens an existing file for reading or read/writing
!===============================================================================

  subroutine file_open(self, filename, mode, serial)

    character(*),      intent(in) :: filename ! name of file to be opened
    character(*),      intent(in) :: mode     ! file access mode 
    logical, optional, intent(in) :: serial    ! processor rank to write from
    class(BinaryOutput) :: self

    ! Check for serial option
    if (present(serial)) then
      self % serial = serial
    else
      self % serial = .true.
     end if

#ifdef HDF5
# ifdef MPI
    if (self % serial) then
      call hdf5_file_open(filename, self % hdf5_fh, mode)
    else
      call hdf5_file_open_parallel(filename, self % hdf5_fh, mode)
    endif
# else
    call hdf5_file_open(filename, self % hdf5_fh, mode)
# endif
#elif MPI
    if (self % serial) then
      ! Check for read/write mode to open, default is read only
      if (mode == 'w') then
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
             STATUS='old', ACCESS='stream', POSITION='append')
      else
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
             STATUS='old', ACCESS='stream')
      end if
    else
      call mpi_open_file(filename, self % unit_fh, mode)
    end if
#else
    ! Check for read/write mode to open, default is read only
    if (mode == 'w') then
      open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
           STATUS='old', ACCESS='stream', POSITION='append')
    else
      open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
           STATUS='old', ACCESS='stream')
    end if
#endif

  end subroutine file_open

!===============================================================================
! FILE_CLOSE closes a file
!===============================================================================

  subroutine file_close(self)

    class(BinaryOutput) :: self

#ifdef HDF5
# ifdef MPI
     call hdf5_file_close(self % hdf5_fh)
# else
     call hdf5_file_close(self % hdf5_fh)
# endif
#elif MPI
     if (self % serial) then
       close(UNIT=self % unit_fh)
     else
       call mpi_close_file(self % unit_fh)
     end if
#else
     close(UNIT=self % unit_fh)
#endif

  end subroutine file_close

!===============================================================================
! OPEN_GROUP call hdf5 routine to open a group within binary output context
!===============================================================================

#ifdef HDF5
  subroutine open_group(self, group)

    character(*), intent(in) :: group ! HDF5 group name
    class(BinaryOutput) :: self

    call hdf5_open_group(self % hdf5_fh, group, self % hdf5_grp)

  end subroutine open_group
#endif

!===============================================================================
! CLOSE_GROUP call hdf5 routine to close a group within binary output context
!===============================================================================

#ifdef HDF5
  subroutine close_group(self)

    class(BinaryOutput) :: self

    call hdf5_close_group(self % hdf5_grp)

  end subroutine close_group
#endif

!===============================================================================
! WRITE_DOUBLE writes double precision scalar data
!===============================================================================

  subroutine write_double(self, buffer, name, group, collect)

    real(8),      intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_double(self % hdf5_grp, name_, buffer)
    else
      call hdf5_write_double_parallel(self % hdf5_grp, name_, buffer, collect_)
    end if
# else
    call hdf5_write_double(self % hdf5_grp, name_, buffer)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer
    else
      call mpi_write_double(self % unit_fh, buffer, collect_)
    end if
#else
    write(self % unit_fh) buffer
#endif

  end subroutine write_double

!===============================================================================
! READ_DOUBLE reads double precision scalar data
!===============================================================================

  subroutine read_double(self, buffer, name, group, collect)

    real(8),      intent(inout)        :: buffer  ! read data to here 
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_double(self % hdf5_grp, name_, buffer)
    else 
      call hdf5_read_double_parallel(self % hdf5_grp, name_, buffer, collect_)
    end if
# else
    call hdf5_read_double(self % hdf5_grp, name_, buffer)
# endif
    ! Check if HDf5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer
    else
      call mpi_read_double(self % unit_fh, buffer, collect_)
    end if
#else
    read(self % unit_fh) buffer
#endif

  end subroutine read_double

!===============================================================================
! WRITE_DOUBLE_1DARRAY writes double precision 1-D array data
!===============================================================================

  subroutine write_double_1Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length    ! length of array to write
    real(8),      intent(in)           :: buffer(:) ! data to write
    character(*), intent(in)           :: name      ! name of data
    character(*), intent(in), optional :: group     ! HDF5 group name
    logical,      intent(in), optional :: collect   ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_double_1Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length)
    else
      call mpi_write_double_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length)
#endif

  end subroutine write_double_1Darray

!===============================================================================
! READ_DOUBLE_1DARRAY reads double precision 1-D array data
!===============================================================================

  subroutine read_double_1Darray(self, buffer, name, group, length, collect)

    integer,        intent(in)           :: length    ! length of array to read
    real(8),        intent(inout)        :: buffer(:) ! read data to here
    character(*),   intent(in)           :: name      ! name of data
    character(*),   intent(in), optional :: group     ! HDF5 group name
    logical,        intent(in), optional :: collect   ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_double_1Darray_parallel(self % hdf5_grp, name_, buffer, &
         length, collect_)
    end if
# else
    call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length)
    else
      call mpi_read_double_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length)
#endif

  end subroutine read_double_1Darray

!===============================================================================
! WRITE_DOUBLE_2DARRAY writes double precision 2-D array data
!===============================================================================

  subroutine write_double_2Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(2) ! dimension of array
    real(8),      intent(in)           :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_double_2Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_double_2Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_double_2Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length(1),1:length(2))
    else
      call mpi_write_double_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2))
#endif

  end subroutine write_double_2Darray

!===============================================================================
! READ_DOUBLE_2DARRAY reads double precision 2-D array data
!===============================================================================

  subroutine read_double_2Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(2) ! dimension of array
    real(8),      intent(inout)        :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_double_2Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_double_2Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_double_2Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2))
    else
      call mpi_read_double_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2))
#endif

  end subroutine read_double_2Darray

!===============================================================================
! WRITE_DOUBLE_3DARRAY writes double precision 3-D array data
!===============================================================================

  subroutine write_double_3Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(3) ! length of each dimension
    real(8),      intent(in)           :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_double_3Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_double_3Darray_parallel(self % hdf5_grp, name_, buffer, length, &
         collect_)
    end if
# else
    call hdf5_write_double_3Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    else
      call mpi_write_double_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
#endif

  end subroutine write_double_3Darray

!===============================================================================
! READ_DOUBLE_3DARRAY reads double precision 3-D array data
!===============================================================================

  subroutine read_double_3Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(3) ! length of each dimension
    real(8),      intent(inout)        :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_double_3Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_double_3Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_double_3Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    else
      call mpi_read_double_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
#endif

  end subroutine read_double_3Darray

!===============================================================================
! WRITE_DOUBLE_4DARRAY writes double precision 4-D array data
!===============================================================================

  subroutine write_double_4Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(4) ! length of each dimension
    real(8),      intent(in)           :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_double_4Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_double_4Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    ! Write the data in serial
    call hdf5_write_double_4Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                         1:length(4))
    else
      call mpi_write_double_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                       1:length(4))
#endif

  end subroutine write_double_4Darray

!===============================================================================
! READ_DOUBLE_4DARRAY reads double precision 4-D array data
!===============================================================================

  subroutine read_double_4Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(4) ! length of each dimension
    real(8),      intent(inout)        :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_double_4Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_double_4Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_double_4Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                        1:length(4))
    else
      call mpi_read_double_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                      1:length(4))
#endif

  end subroutine read_double_4Darray

!===============================================================================
! WRITE_INTEGER writes integer precision scalar data
!===============================================================================

  subroutine write_integer(self, buffer, name, group, collect)

    integer,      intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_integer(self % hdf5_grp, name_, buffer)
    else
      call hdf5_write_integer_parallel(self % hdf5_grp, name_, buffer, collect_)
    end if
# else
    call hdf5_write_integer(self % hdf5_grp, name_, buffer)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer
    else
      call mpi_write_integer(self % unit_fh, buffer, collect_)
    end if
#else
    write(self % unit_fh) buffer
#endif

  end subroutine write_integer

!===============================================================================
! READ_INTEGER reads integer precision scalar data
!===============================================================================

  subroutine read_integer(self, buffer, name, group, collect)

    integer,      intent(inout)        :: buffer  ! read data to here 
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_integer(self % hdf5_grp, name_, buffer)
    else
      call hdf5_read_integer_parallel(self % hdf5_grp, name_, buffer, collect_)
    end if
# else
    call hdf5_read_integer(self % hdf5_grp, name_, buffer)
# endif
    ! Check if HDf5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer
    else
      call mpi_read_integer(self % unit_fh, buffer, collect_)
    end if
#else
    read(self % unit_fh) buffer
#endif

  end subroutine read_integer

!===============================================================================
! WRITE_INTEGER_1DARRAY writes integer precision 1-D array data
!===============================================================================

  subroutine write_integer_1Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length    ! length of array to write
    integer,      intent(in)           :: buffer(:) ! data to write
    character(*), intent(in)           :: name      ! name of data
    character(*), intent(in), optional :: group     ! HDF5 group name
    logical,      intent(in), optional :: collect   ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_integer_1Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_integer_1Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_integer_1Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length)
    else
      call mpi_write_integer_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length)
#endif

  end subroutine write_integer_1Darray

!===============================================================================
! READ_INTEGER_1DARRAY reads integer precision 1-D array data
!===============================================================================

  subroutine read_integer_1Darray(self, buffer, name, group, length, collect)

    integer,        intent(in)           :: length    ! length of array to read
    integer,        intent(inout)        :: buffer(:) ! read data to here
    character(*),   intent(in)           :: name      ! name of data
    character(*),   intent(in), optional :: group     ! HDF5 group name
    logical,        intent(in), optional :: collect   ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_integer_1Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_integer_1Darray_parallel(self % hdf5_grp, name_, buffer, &
           length, collect_)
    end if
# else
    ! Read the data in serial
    call hdf5_read_integer_1Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length)
    else
      call mpi_read_integer_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length)
#endif

  end subroutine read_integer_1Darray

!===============================================================================
! WRITE_INTEGER_2DARRAY writes integer precision 2-D array data
!===============================================================================

  subroutine write_integer_2Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(2) ! dimension of array
    integer,      intent(in)           :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_integer_2Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_integer_2Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_integer_2Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then 
      write(self % unit_fh) buffer(1:length(1),1:length(2))
    else
      call mpi_write_integer_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2))
#endif

  end subroutine write_integer_2Darray

!===============================================================================
! READ_INTEGER_2DARRAY reads integer precision 2-D array data
!===============================================================================

  subroutine read_integer_2Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(2) ! dimension of array
    integer,      intent(inout)        :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_integer_2Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_integer_2Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_integer_2Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2))
    else
      call mpi_read_integer_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2))
#endif

  end subroutine read_integer_2Darray

!===============================================================================
! WRITE_INTEGER_3DARRAY writes integer precision 3-D array data
!===============================================================================

  subroutine write_integer_3Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(3) ! length of each dimension
    integer,      intent(in)           :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_integer_3Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_integer_3Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_integer_3Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    else
      call mpi_write_integer_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
#endif

  end subroutine write_integer_3Darray

!===============================================================================
! READ_INTEGER_3DARRAY reads integer precision 3-D array data
!===============================================================================

  subroutine read_integer_3Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(3) ! length of each dimension
    integer,      intent(inout)        :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_integer_3Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_integer_3Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_integer_3Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    else
      call mpi_read_integer_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
#endif

  end subroutine read_integer_3Darray

!===============================================================================
! WRITE_INTEGER_4DARRAY writes integer precision 4-D array data
!===============================================================================

  subroutine write_integer_4Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(4) ! length of each dimension
    integer,      intent(in)           :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_integer_4Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_write_integer_4Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_write_integer_4Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                         1:length(4))
    else
      call mpi_write_integer_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                       1:length(4))
#endif

  end subroutine write_integer_4Darray

!===============================================================================
! READ_INTEGER_4DARRAY reads integer precision 4-D array data
!===============================================================================

  subroutine read_integer_4Darray(self, buffer, name, group, length, collect)

    integer,      intent(in)           :: length(4) ! length of each dimension
    integer,      intent(inout)        :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_integer_4Darray(self % hdf5_grp, name_, buffer, length)
    else
      call hdf5_read_integer_4Darray_parallel(self % hdf5_grp, name_, buffer, length, &
           collect_)
    end if
# else
    call hdf5_read_integer_4Darray(self % hdf5_grp, name_, buffer, length)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                        1:length(4))
    else
      call mpi_read_integer_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
                      1:length(4))
#endif

  end subroutine read_integer_4Darray

!===============================================================================
! WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine write_long(self, buffer, name, group, collect)

    integer(8),   intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name of data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_long(self % hdf5_grp, name_, buffer, hdf5_integer8_t)
    else
      call hdf5_write_long_parallel(self % hdf5_grp, name_, buffer, &
           hdf5_integer8_t, collect_)
    end if
# else
    call hdf5_write_long(self % hdf5_grp, name_, buffer, hdf5_integer8_t)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer
    else
      call mpi_write_long(self % unit_fh, buffer, collect_)
    end if
#else
    write(self % unit_fh) buffer
#endif

  end subroutine write_long

!===============================================================================
! READ_LONG reads long integer scalar data
!===============================================================================

  subroutine read_long(self, buffer, name, group, collect)

    integer(8),   intent(inout)        :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name of data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    logical :: collect_

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_long(self % hdf5_grp, name_, buffer, hdf5_integer8_t)
    else
      call hdf5_read_long_parallel(self % hdf5_grp, name_, buffer, &
           hdf5_integer8_t, collect_)
    end if
# else
    call hdf5_read_long(self % hdf5_grp, name_, buffer, hdf5_integer8_t)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer
    else
      call mpi_read_long(self % unit_fh, buffer, collect_)
    end if
#else
    read(self % unit_fh) buffer
#endif

  end subroutine read_long

!===============================================================================
! WRITE_STRING writes string data
!===============================================================================

  subroutine write_string(self, buffer, name, group, collect)

    character(*), intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name of data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    integer :: n
    logical :: collect_

    ! Get string length
    n = len_trim(buffer)

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if

#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_write_string(self % hdf5_grp, name_, buffer, n)
    else
      call hdf5_write_string_parallel(self % hdf5_grp, name_, buffer, n, collect_)
    end if
# else
    ! Write the data
    call hdf5_write_string(self % hdf5_grp, name_, buffer, n)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      write(self % unit_fh) buffer
    else
      call mpi_write_string(self % unit_fh, buffer, n, collect_)
    end if
#else
    write(self % unit_fh) buffer
#endif

  end subroutine write_string

!===============================================================================
! READ_STRING reads string data
!===============================================================================

  subroutine read_string(self, buffer, name, group, collect)

    character(*), intent(inout)        :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name of data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name
    integer :: n
    logical :: collect_

    ! Get string length
    n = len(buffer)

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

    ! Set up collective vs. independent I/O
    if (present(collect)) then
      collect_ = collect
    else
      collect_ = .true.
    end if


#ifdef HDF5
    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif
# ifdef MPI
    if (self % serial) then
      call hdf5_read_string(self % hdf5_grp, name_, buffer, n)
    else
      call hdf5_read_string_parallel(self % hdf5_grp, name_, buffer, n, collect_)
    end if
# else
    call hdf5_read_string(self % hdf5_grp, name_, buffer, n)
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      read(self % unit_fh) buffer
    else
      call mpi_read_string(self % unit_fh, buffer, n, collect_)
    end if
#else
    read(self % unit_fh) buffer
#endif

  end subroutine read_string

!===============================================================================
! WRITE_ATTRIBUTE_STRING
!===============================================================================

#ifdef HDF5
  subroutine write_attribute_string(self, var, attr_type, attr_str, group)

    character(*), intent(in)           :: var       ! variable name for attr
    character(*), intent(in)           :: attr_type ! attr identifier type
    character(*), intent(in)           :: attr_str  ! string for attr id type
    character(*), intent(in), optional :: group     ! HDF5 group name
    class(BinaryOutput) :: self

    ! Check if HDF5 group should be created/opened
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    endif

    ! Write the attribute string
    call hdf5_write_attribute_string(self % hdf5_grp, var, attr_type, attr_str)

    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)

  end subroutine write_attribute_string
#endif

!===============================================================================
! WRITE_TALLY_RESULT writes an OpenMC TallyResult type
!===============================================================================

  subroutine write_tally_result(self, buffer, name, group, n1, n2)

    character(*),      intent(in), optional :: group   ! HDF5 group name
    character(*),      intent(in)           :: name    ! name of data
    integer,           intent(in)           :: n1, n2  ! TallyResult dims
    type(TallyResult), intent(in), target   :: buffer(n1, n2) ! data to write
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name

#ifndef HDF5
    integer :: j,k ! iteration counters
#endif

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

#ifdef HDF5

    ! Open up sub-group if present
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    end if

    ! Set overall size of vector to write
    dims1(1) = n1*n2 

    ! Create up a dataspace for size
    call h5screate_simple_f(1, dims1, dspace, hdf5_err)

    ! Create the dataset
    call h5dcreate_f(self % hdf5_grp, name_, hdf5_tallyresult_t, dspace, dset, &
         hdf5_err)

    ! Set pointer to first value and write
    f_ptr = c_loc(buffer(1,1))
    call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

    ! Close ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    if (present(group)) then
      call hdf5_close_group(self % hdf5_grp)
    end if

#else

    ! Write out tally buffer
    do k = 1, n2
      do j = 1, n1
        write(self % unit_fh) buffer(j,k) % sum
        write(self % unit_fh) buffer(j,k) % sum_sq
      end do
    end do

#endif 
   
  end subroutine write_tally_result

!===============================================================================
! READ_TALLY_RESULT reads OpenMC TallyResult data
!===============================================================================

  subroutine read_tally_result(self, buffer, name, group, n1, n2)

    character(*),      intent(in), optional  :: group  ! HDF5 group name
    character(*),      intent(in)            :: name   ! name of data
    integer,           intent(in)            :: n1, n2 ! TallyResult dims
    type(TallyResult), intent(inout), target :: buffer(n1, n2) ! read data here
    class(BinaryOutput) :: self

    character(len=MAX_WORD_LEN) :: name_  ! HDF5 dataset name
    character(len=MAX_WORD_LEN) :: group_ ! HDF5 group name

#ifndef HDF5
# ifndef MPI
    integer :: j,k ! iteration counters
# endif
#endif

    ! Set name
    name_ = trim(name)

    ! Set group
    if (present(group)) then
      group_ = trim(group)
    end if

#ifdef HDF5

    ! Open up sub-group if present
    if (present(group)) then
      call hdf5_open_group(self % hdf5_fh, group_, self % hdf5_grp)
    else
      self % hdf5_grp = self % hdf5_fh
    end if

    ! Open the dataset
    call h5dopen_f(self % hdf5_grp, name, dset, hdf5_err)

    ! Set pointer to first value and write
    f_ptr = c_loc(buffer(1,1))
    call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

    ! Close ids
    call h5dclose_f(dset, hdf5_err)
    if (present(group)) call hdf5_close_group(self % hdf5_grp)

# elif MPI

    ! Write out tally buffer
    call MPI_FILE_READ(self % unit_fh, buffer, n1*n2, MPI_TALLYRESULT, &
         MPI_STATUS_IGNORE, mpiio_err)

#else

    ! Read tally result
    do k = 1, n2
      do j = 1, n1
        read(self % unit_fh) buffer(j,k) % sum
        read(self % unit_fh) buffer(j,k) % sum_sq
      end do
    end do

#endif 
   
  end subroutine read_tally_result

!===============================================================================
! WRITE_SOURCE_BANK writes OpenMC source_bank data
!===============================================================================

  subroutine write_source_bank(self)

    class(BinaryOutput) :: self

#ifdef MPI
# ifndef HDF5
    integer(MPI_OFFSET_KIND) :: offset           ! offset of data
    integer                  :: size_offset_kind ! the data offset kind
    integer                  :: size_bank        ! size of bank to write
# endif
# ifdef HDF5
    integer(8)               :: offset(1)        ! source data offset
# endif
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims1(1) = n_particles
    hdf5_rank = 1

    ! Create that dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create the dataset for that dataspace
    call h5dcreate_f(self % hdf5_fh, "source_bank", hdf5_bank_t, dspace, dset, hdf5_err)

    ! Close the dataspace
    call h5sclose_f(dspace, hdf5_err)

    ! Create another data space but for each proc individually
    dims1(1) = work
    call h5screate_simple_f(hdf5_rank, dims1, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims1, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Write data to file in parallel
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id = dspace, mem_space_id = memspace, &
         xfer_prp = plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

# else

    ! Set size
    dims1(1) = work
    hdf5_rank = 1

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(self % hdf5_fh, "source_bank", hdf5_bank_t, &
         dspace, dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Write dataset to file
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

# endif 

#elif MPI

    ! Get current offset for master 
    if (master) call MPI_FILE_GET_POSITION(self % unit_fh, offset, mpiio_err)

    ! Determine offset on master process and broadcast to all processors
    call MPI_SIZEOF(offset, size_offset_kind, mpi_err)
    select case (size_offset_kind)
    case (4)
      call MPI_BCAST(offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    case (8)
      call MPI_BCAST(offset, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_err)
    end select

    ! Set the proper offset for source data on this processor
    call MPI_TYPE_SIZE(MPI_BANK, size_bank, mpi_err)
    offset = offset + size_bank*work_index(rank)

    ! Write all source sites
    call MPI_FILE_WRITE_AT(self % unit_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpiio_err)

#else

    ! Write out source sites
    write(self % unit_fh) source_bank

#endif

  end subroutine write_source_bank

!===============================================================================
! READ_SOURCE_BANK reads OpenMC source_bank data
!===============================================================================

  subroutine read_source_bank(self)

    class(BinaryOutput) :: self

#ifdef MPI
# ifndef HDF5
    integer(MPI_OFFSET_KIND) :: offset           ! offset of data
    integer                  :: size_bank        ! size of bank to read
# endif
# ifdef HDF5
    integer(8)               :: offset(1)        ! offset of data
# endif
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims1(1) = n_particles
    hdf5_rank = 1

    ! Open the dataset
    call h5dopen_f(self % hdf5_fh, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims1(1) = work
    call h5screate_simple_f(hdf5_rank, dims1, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims1, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Read data from file in parallel
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id = dspace, mem_space_id = memspace, &
         xfer_prp = plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

# else

    ! Open dataset
    call h5dopen_f(self % hdf5_fh, "source_bank", dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Read dataset from file
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)

# endif 

#elif MPI

    ! Go to the end of the file to set file pointer
    offset = 0
    call MPI_FILE_SEEK(self % unit_fh, offset, MPI_SEEK_END, &
         mpiio_err)

    ! Get current offset (will be at EOF) 
    call MPI_FILE_GET_POSITION(self % unit_fh, offset, mpiio_err)

    ! Get the size of the source bank on all procs
    call MPI_TYPE_SIZE(MPI_BANK, size_bank, mpi_err)

    ! Calculate offset where the source bank will begin
    offset = offset - n_particles*size_bank

    ! Set the proper offset for source data on this processor
    offset = offset + size_bank*work_index(rank)

    ! Write all source sites
    call MPI_FILE_READ_AT(self % unit_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpiio_err)

#else

    ! Write out source sites
    read(self % unit_fh) source_bank

#endif

  end subroutine read_source_bank

end module output_interface
