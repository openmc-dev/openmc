module output_interface

  use constants

#ifdef HDF5
  use h5lt
  use hdf5_interface
#endif
#ifdef MPI
  use mpiio_interface
#endif

  implicit none
  private

  type, public :: BinaryOutput
    ! Compilation specific data
#ifdef HDF5
    integer(HID_T) :: hdf5_fh
    integer(HID_T) :: hdf5_grp
#else
    integer :: unit_fh
# endif
    logical :: serial ! Serial I/O when using MPI/PHDF5
    logical :: direct_access ! If the file is opened with access='direct'
    integer :: record_len ! RECL for direct unformatted access
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

  subroutine file_create(self, filename, serial, record_len)

    character(*),      intent(in) :: filename  ! name of file to be created
    logical, optional, intent(in) :: serial    ! processor rank to write from
    integer, optional, intent(in) :: record_len ! RECL
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
      if (present(record_len)) then
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
             STATUS='replace', ACCESS='direct', RECL=record_len)
      else
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
             STATUS='replace', ACCESS='stream')
      end if
    else
      call mpi_create_file(filename, self % unit_fh)
    end if
#else
    if (present(record_len)) then
      open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
           STATUS='replace', ACCESS='direct', RECL=record_len)
    else
      open(NEWUNIT=self % unit_fh, FILE=filename, ACTION="write", &
           STATUS='replace', ACCESS='stream')
    end if
#endif

  end subroutine file_create

!===============================================================================
! FILE_OPEN opens an existing file for reading or read/writing
!===============================================================================

  subroutine file_open(self, filename, mode, serial, direct_access, record_len)

    character(*),      intent(in) :: filename ! name of file to be opened
    character(*),      intent(in) :: mode     ! file access mode 
    logical, optional, intent(in) :: serial   ! processor rank to write from
    logical, optional, intent(in) :: direct_access
    integer, optional, intent(in) :: record_len
    class(BinaryOutput) :: self

    ! Check for serial option
    if (present(serial)) then
      self % serial = serial
    else
      self % serial = .true.
    end if

    ! Check for direct_access option
    if (present(direct_access)) then
      self % direct_access = direct_access
      self % record_len = 4
    else
      self % direct_access = .false.
    end if

    ! Check for record_len option
    if (present(record_len)) then
      self % record_len = record_len
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
        if (self % direct_access) then
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
               STATUS='old', ACCESS='direct', RECL=self % record_len)
        else
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
               STATUS='old', ACCESS='stream', POSITION='append')
        end if
      else
        if (self % direct_access) then
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
               STATUS='old', ACCESS='direct', RECL=self % record_len)
        else
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
               STATUS='old', ACCESS='stream')
        endif
      end if
    else
      call mpi_open_file(filename, self % unit_fh, mode)
    end if
#else
    ! Check for read/write mode to open, default is read only
    if (mode == 'w') then
        if (self % direct_access) then
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
               STATUS='old', ACCESS='direct', RECL=self % record_len)
        else
          open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='write', &
               STATUS='old', ACCESS='stream', POSITION='append')
        end if
    else
      if (self % direct_access) then
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
             STATUS='old', ACCESS='direct', RECL=self % record_len)
      else
        open(NEWUNIT=self % unit_fh, FILE=filename, ACTION='read', &
             STATUS='old', ACCESS='stream')
      end if
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

  subroutine write_double(self, buffer, name, group, collect, record)

    real(8),      intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record     ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer
      else
        write(self % unit_fh) buffer
      endif
    else
      call mpi_write_double(self % unit_fh, buffer, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer
    else
      write(self % unit_fh) buffer
    endif
#endif

  end subroutine write_double

!===============================================================================
! READ_DOUBLE reads double precision scalar data
!===============================================================================

  subroutine read_double(self, buffer, name, group, collect, record)

    real(8),      intent(inout)        :: buffer  ! read data to here 
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer
      else
        read(self % unit_fh) buffer
      end if
    else
      call mpi_read_double(self % unit_fh, buffer, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer
    else
      read(self % unit_fh) buffer
    end if
#endif

  end subroutine read_double

!===============================================================================
! WRITE_DOUBLE_1DARRAY writes double precision 1-D array data
!===============================================================================

  subroutine write_double_1Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length    ! length of array to write
    real(8),      intent(in)           :: buffer(:) ! data to write
    character(*), intent(in)           :: name      ! name of data
    character(*), intent(in), optional :: group     ! HDF5 group name
    logical,      intent(in), optional :: collect   ! collective I/O
    integer,      intent(in), optional :: record    ! REC
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
      if (present(record)) then
        call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length, record)
      else
        call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length)
      end if
    else
      if (present(record)) then
        call hdf5_write_double_1Darray_parallel(self % hdf5_grp, name_, buffer, length, &
             collect_, record)
      else
        call hdf5_write_double_1Darray_parallel(self % hdf5_grp, name_, buffer, length, &
             collect_)
      end if
    end if
# else
    if (present(record)) then
      call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length, record)
    else
      call hdf5_write_double_1Darray(self % hdf5_grp, name_, buffer, length)
    end if
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length)
      else
        write(self % unit_fh) buffer(1:length)
      endif
    else
      call mpi_write_double_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length)
    else
      write(self % unit_fh) buffer(1:length)
    endif
#endif

  end subroutine write_double_1Darray

!===============================================================================
! READ_DOUBLE_1DARRAY reads double precision 1-D array data
!===============================================================================

  subroutine read_double_1Darray(self, buffer, name, group, length, collect, &
                                 record)

    integer,        intent(in)           :: length    ! length of array to read
    real(8),        intent(inout)        :: buffer(:) ! read data to here
    character(*),   intent(in)           :: name      ! name of data
    character(*),   intent(in), optional :: group     ! HDF5 group name
    logical,        intent(in), optional :: collect   ! collective I/O
    integer,        intent(in), optional :: record    ! REC
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
      if (present(record)) then
        call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length, record)
      else
        call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length)
      end if
    else
      if (present(record)) then
        call hdf5_read_double_1Darray_parallel(self % hdf5_grp, name_, buffer, &
           length, collect_, record)
      else
        call hdf5_read_double_1Darray_parallel(self % hdf5_grp, name_, buffer, &
           length, collect_)
      end if
    end if
# else
    if (present(record)) then
      call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length, record)
    else
      call hdf5_read_double_1Darray(self % hdf5_grp, name_, buffer, length)
    end if
# endif
    ! Check if HDF5 group should be closed
    if (present(group)) call hdf5_close_group(self % hdf5_grp)
#elif MPI
    if (self % serial) then
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length)
      else
        read(self % unit_fh) buffer(1:length)
      end if
    else
      call mpi_read_double_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length)
    else
      read(self % unit_fh) buffer(1:length)
    end if
#endif

  end subroutine read_double_1Darray

!===============================================================================
! WRITE_DOUBLE_2DARRAY writes double precision 2-D array data
!===============================================================================

  subroutine write_double_2Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(2) ! dimension of array
    real(8),      intent(in)           :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record     ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2))
      end if
    else
      call mpi_write_double_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2))
    end if
#endif

  end subroutine write_double_2Darray

!===============================================================================
! READ_DOUBLE_2DARRAY reads double precision 2-D array data
!===============================================================================

  subroutine read_double_2Darray(self, buffer, name, group, length, collect, &
                                 record)

    integer,      intent(in)           :: length(2) ! dimension of array
    real(8),      intent(inout)        :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2))
      end if
    else
      call mpi_read_double_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2))
    end if
#endif

  end subroutine read_double_2Darray

!===============================================================================
! WRITE_DOUBLE_3DARRAY writes double precision 3-D array data
!===============================================================================

  subroutine write_double_3Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(3) ! length of each dimension
    real(8),      intent(in)           :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
      end if
    else
      call mpi_write_double_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    end if
#endif

  end subroutine write_double_3Darray

!===============================================================================
! READ_DOUBLE_3DARRAY reads double precision 3-D array data
!===============================================================================

  subroutine read_double_3Darray(self, buffer, name, group, length, collect, &
                                 record)

    integer,      intent(in)           :: length(3) ! length of each dimension
    real(8),      intent(inout)        :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
      end if
    else
      call mpi_read_double_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    end if
#endif

  end subroutine read_double_3Darray

!===============================================================================
! WRITE_DOUBLE_4DARRAY writes double precision 4-D array data
!===============================================================================

  subroutine write_double_4Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(4) ! length of each dimension
    real(8),      intent(in)           :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
            1:length(3), 1:length(4))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
            1:length(4))
      end if
    else
      call mpi_write_double_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
          1:length(3), 1:length(4))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
          1:length(4))
    end if
#endif

  end subroutine write_double_4Darray

!===============================================================================
! READ_DOUBLE_4DARRAY reads double precision 4-D array data
!===============================================================================

  subroutine read_double_4Darray(self, buffer, name, group, length, collect, &
                                 record)

    integer,      intent(in)           :: length(4) ! length of each dimension
    real(8),      intent(inout)        :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
            1:length(3), 1:length(4))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
            1:length(4))
      end if
    else
      call mpi_read_double_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
          1:length(3), 1:length(4))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
          1:length(4))
    end if
#endif

  end subroutine read_double_4Darray

!===============================================================================
! WRITE_INTEGER writes integer precision scalar data
!===============================================================================

  subroutine write_integer(self, buffer, name, group, collect, record)

    integer,      intent(in)           :: buffer  ! data to write
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer
      else
        write(self % unit_fh) buffer
      end if
    else
      call mpi_write_integer(self % unit_fh, buffer, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer
    else
      write(self % unit_fh) buffer
    end if
#endif

  end subroutine write_integer

!===============================================================================
! READ_INTEGER reads integer precision scalar data
!===============================================================================

  subroutine read_integer(self, buffer, name, group, collect, record)

    integer,      intent(inout)        :: buffer  ! read data to here 
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), optional :: group   ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer
      else
        read(self % unit_fh) buffer
      end if
    else
      call mpi_read_integer(self % unit_fh, buffer, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer
    else
      read(self % unit_fh) buffer
    end if
#endif

  end subroutine read_integer

!===============================================================================
! WRITE_INTEGER_1DARRAY writes integer precision 1-D array data
!===============================================================================

  subroutine write_integer_1Darray(self, buffer, name, group, length, collect, &
                                   record)

    integer,      intent(in)           :: length    ! length of array to write
    integer,      intent(in)           :: buffer(:) ! data to write
    character(*), intent(in)           :: name      ! name of data
    character(*), intent(in), optional :: group     ! HDF5 group name
    logical,      intent(in), optional :: collect   ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length)
      else
        write(self % unit_fh) buffer(1:length)
      end if
    else
      call mpi_write_integer_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length)
    else
      write(self % unit_fh) buffer(1:length)
    end if
#endif

  end subroutine write_integer_1Darray

!===============================================================================
! READ_INTEGER_1DARRAY reads integer precision 1-D array data
!===============================================================================

  subroutine read_integer_1Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,        intent(in)           :: length    ! length of array to read
    integer,        intent(inout)        :: buffer(:) ! read data to here
    character(*),   intent(in)           :: name      ! name of data
    character(*),   intent(in), optional :: group     ! HDF5 group name
    logical,        intent(in), optional :: collect   ! collective I/O
    integer,        intent(in), optional :: record    ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length)
      else
        read(self % unit_fh) buffer(1:length)
      end if
    else
      call mpi_read_integer_1Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length)
    else
      read(self % unit_fh) buffer(1:length)
    end if
#endif

  end subroutine read_integer_1Darray

!===============================================================================
! WRITE_INTEGER_2DARRAY writes integer precision 2-D array data
!===============================================================================

  subroutine write_integer_2Darray(self, buffer, name, group, length, collect, &
                                   record)

    integer,      intent(in)           :: length(2) ! dimension of array
    integer,      intent(in)           :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2))
      end if
    else
      call mpi_write_integer_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2))
    end if
#endif

  end subroutine write_integer_2Darray

!===============================================================================
! READ_INTEGER_2DARRAY reads integer precision 2-D array data
!===============================================================================

  subroutine read_integer_2Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(2) ! dimension of array
    integer,      intent(inout)        :: buffer(length(1),length(2)) ! the data
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2))
      end if
    else
      call mpi_read_integer_2Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2))
    end if
#endif

  end subroutine read_integer_2Darray

!===============================================================================
! WRITE_INTEGER_3DARRAY writes integer precision 3-D array data
!===============================================================================

  subroutine write_integer_3Darray(self, buffer, name, group, length, collect, &
                                   record)

    integer,      intent(in)           :: length(3) ! length of each dimension
    integer,      intent(in)           :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
      end if
    else
      call mpi_write_integer_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    end if
#endif

  end subroutine write_integer_3Darray

!===============================================================================
! READ_INTEGER_3DARRAY reads integer precision 3-D array data
!===============================================================================

  subroutine read_integer_3Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(3) ! length of each dimension
    integer,      intent(inout)        :: buffer(length(1),length(2),length(3))        
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
      end if
    else
      call mpi_read_integer_3Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2),1:length(3))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3))
    end if
#endif

  end subroutine read_integer_3Darray

!===============================================================================
! WRITE_INTEGER_4DARRAY writes integer precision 4-D array data
!===============================================================================

  subroutine write_integer_4Darray(self, buffer, name, group, length, collect, &
                                   record)

    integer,      intent(in)           :: length(4) ! length of each dimension
    integer,      intent(in)           :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
            1:length(3), 1:length(4))
      else
        write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
            1:length(4))
      end if
    else
      call mpi_write_integer_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      write(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
          1:length(3), 1:length(4))
    else
      write(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
          1:length(4))
    end if
#endif

  end subroutine write_integer_4Darray

!===============================================================================
! READ_INTEGER_4DARRAY reads integer precision 4-D array data
!===============================================================================

  subroutine read_integer_4Darray(self, buffer, name, group, length, collect, &
                                  record)

    integer,      intent(in)           :: length(4) ! length of each dimension
    integer,      intent(inout)        :: buffer(length(1),length(2),&
                                                 length(3),length(4))
    character(*), intent(in)           :: name ! name of data
    character(*), intent(in), optional :: group ! HDF5 group name
    logical,      intent(in), optional :: collect ! collective I/O
    integer,      intent(in), optional :: record  ! REC
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
      if (present(record)) then
        read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
            1:length(3), 1:length(4))
      else
        read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
            1:length(4))
      end if
    else
      call mpi_read_integer_4Darray(self % unit_fh, buffer, length, collect_)
    end if
#else
    if (present(record)) then
      read(self % unit_fh, REC=record) buffer(1:length(1),1:length(2), &
          1:length(3), 1:length(4))
    else
      read(self % unit_fh) buffer(1:length(1),1:length(2),1:length(3), &
          1:length(4))
    end if
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

end module output_interface
