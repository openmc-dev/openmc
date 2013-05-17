module output_interface

  use constants
  use global
  use tally_header,  only: TallyResult

#ifdef HDF5
  use hdf5_interface
#elif MPI
  use mpi_interface
#endif

  implicit none

  integer :: mpi_fh

  interface write_data
    module procedure write_double
    module procedure write_double_1Darray
    module procedure write_double_2Darray
    module procedure write_double_3Darray
    module procedure write_integer
    module procedure write_integer_1Darray
    module procedure write_integer_2Darray
    module procedure write_integer_3Darray
    module procedure write_long
    module procedure write_string
  end interface write_data

  interface read_data
    module procedure read_double
    module procedure read_double_1Darray
    module procedure read_integer
    module procedure read_integer_1Darray
    module procedure read_long
    module procedure read_string
  end interface read_data

contains

!===============================================================================
! FILE_CREATE
!===============================================================================

  subroutine file_create(filename, fh_str)

    character(*) :: filename
    character(*) :: fh_str

#ifdef HDF5
  filename = trim(filename) // '.h5'
#else
  filename = trim(filename) // '.binary'
#endif

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'serial') then
      if(master) call hdf5_file_create(filename, hdf5_fh)
    else
      call hdf5_parallel_file_create(filename, hdf5_fh)
    endif
# else
    call hdf5_file_create(filename, hdf5_fh)
# endif
#elif MPI
    call mpi_create_file(filename, mpi_fh)
#else
    open(UNIT=UNIT_OUTPUT, FILE=filename, ACTION="write", &
         STATUS='replace', ACCESS='stream')
#endif

  end subroutine file_create

!===============================================================================
! FILE_OPEN
!===============================================================================

  subroutine file_open(filename, fh_str, mode)

    character(*) :: filename
    character(*) :: fh_str
    character(*) :: mode

#ifdef HDF5
  filename = trim(filename) // '.h5'
#else
  filename = trim(filename) // '.binary'
#endif

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'serial') then
      if (master) call hdf5_file_open(filename, hdf5_fh, mode)
    else
      call hdf5_parallel_file_open(filename, hdf5_fh, mode)
    endif
# else
    call hdf5_file_open(filename, hdf5_fh, mode)
# endif
#elif MPI
    call mpi_open_file(filename, mpi_fh, mode)
#else
    if (mode == 'rw') then
      open(UNIT=UNIT_OUTPUT, FILE=filename, ACTION='write', &
           STATUS='old', ACCESS='stream')
    else
      open(UNIT=UNIT_OUTPUT, FILE=filename, ACTION='read', &
           STATUS='old', ACCESS='stream')
    end if
#endif

  end subroutine file_open

!===============================================================================
! FILE_CLOSE
!===============================================================================

  subroutine file_close(fh_str)

    character(*) :: fh_str

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'serial') then
     if(master) call hdf5_file_close(hdf5_fh)
   else
     call hdf5_file_close(hdf5_fh)
    endif
# else
     call hdf5_file_close(hdf5_fh)
# endif
#elif MPI
     call mpi_close_file(mpi_fh)
#else
     close(UNIT=UNIT_OUTPUT)
#endif

  end subroutine file_close

!===============================================================================
! WRITE_DOUBLE
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
    call mpi_write_double(mpi_fh, buffer)
#else
    write(UNIT_OUTPUT) buffer
#endif

  end subroutine write_double

!===============================================================================
! READ_INTEGER
!===============================================================================

  subroutine read_double(buffer, name, group)

    real(8),        intent(inout)        :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_double(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_read_double(mpi_fh, buffer)
#else
    read(UNIT_OUTPUT) buffer
#endif

  end subroutine read_double

!===============================================================================
! WRITE_DOUBLE_1DARRAY
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
    call mpi_write_double_1Darray(mpi_fh, buffer, length)
#else
    write(UNIT_OUTPUT) buffer(1:length)
#endif

  end subroutine write_double_1Darray

!===============================================================================
! READ_DOUBLE_1DARRAY
!===============================================================================

  subroutine read_double_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    real(8),        intent(inout)        :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_double_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_read_double_1Darray(mpi_fh, buffer, length)
#else
    read(UNIT_OUTPUT) buffer(1:length)
#endif

  end subroutine read_double_1Darray

!===============================================================================
! WRITE_DOUBLE_2DARRAY
!===============================================================================

  subroutine write_double_2Darray(buffer, name, group, length)

    integer,        intent(in)           :: length(2)
    real(8),        intent(in)           :: buffer(length(1),length(2))
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double_2Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#endif

  end subroutine write_double_2Darray

!===============================================================================
! WRITE_DOUBLE_3DARRAY
!===============================================================================

  subroutine write_double_3Darray(buffer, name, group, length)

    integer,        intent(in)           :: length(3)
    real(8),        intent(in)           :: buffer(length(1),length(2), &
                                                   length(3))
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double_3Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#endif

  end subroutine write_double_3Darray

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
    call mpi_write_integer(mpi_fh, buffer)
#else
    write(UNIT_OUTPUT) buffer
#endif

  end subroutine write_integer

!===============================================================================
! READ_INTEGER
!===============================================================================

  subroutine read_integer(buffer, name, group)

    integer,        intent(inout)        :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_integer(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_read_integer(mpi_fh, buffer)
#else
    read(UNIT_OUTPUT) buffer
#endif

  end subroutine read_integer

!===============================================================================
! WRITE_INTEGER_1DARRAY
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
    call mpi_write_integer_1Darray(mpi_fh, buffer, length)
#else
    write(UNIT_OUTPUT) buffer(1:length)
#endif

  end subroutine write_integer_1Darray

!===============================================================================
! READ_INTEGER_1DARRAY
!===============================================================================

  subroutine read_integer_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    integer,        intent(inout)        :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_integer_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_read_integer_1Darray(mpi_fh, buffer, length)
#else
    read(UNIT_OUTPUT) buffer(1:length)
#endif

  end subroutine read_integer_1Darray

!===============================================================================
! WRITE_INTEGER_2DARRAY
!===============================================================================

  subroutine write_integer_2Darray(buffer, name, group, length)

    integer,        intent(in)           :: length(2)
    integer,        intent(in)           :: buffer(length(1),length(2))
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer_2Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#endif

  end subroutine write_integer_2Darray

!===============================================================================
! WRITE_DOUBLE_3DARRAY
!===============================================================================

  subroutine write_integer_3Darray(buffer, name, group, length)

    integer,        intent(in)           :: length(3)
    integer,        intent(in)           :: buffer(length(1),length(2), &
                                                   length(3))
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer_3Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#endif

  end subroutine write_integer_3Darray

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
    call hdf5_write_long(temp_group, name, buffer, hdf5_integer8_t)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_write_long(mpi_fh, buffer)
#else
    write(UNIT_OUTPUT) buffer
#endif

  end subroutine write_long

!===============================================================================
! READ_LONG
!===============================================================================

  subroutine read_long(buffer, name, group)

    integer(8),     intent(inout)        :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_long(temp_group, name, buffer, hdf5_integer8_t)
    if (present(group)) call hdf5_close_group()
#elif MPI
    call mpi_read_long(mpi_fh, buffer)
#else
    read(UNIT_OUTPUT) buffer
#endif

  end subroutine read_long

!===============================================================================
! WRITE_STRING
!===============================================================================

  subroutine write_string(buffer, name, group)

    character(*), intent(in)           :: buffer
    character(*), intent(in)           :: name
    character(*), intent(in), optional :: group

#ifndef HDF5
# ifdef MPI
    integer :: n
# endif
#endif

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_string(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI
    n = len(buffer)
    call mpi_write_string(mpi_fh, buffer, n)
#else
    write(UNIT_OUTPUT) buffer
#endif

  end subroutine write_string

!===============================================================================
! WRITE_STRING
!===============================================================================

  subroutine read_string(buffer, name, group)

    character(*), intent(inout)        :: buffer
    character(*), intent(in)           :: name
    character(*), intent(in), optional :: group

#ifndef HDF5
# ifdef MPI
    integer :: n
# endif
#endif

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_read_string(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI
    n = len(buffer)
    call mpi_read_string(mpi_fh, buffer, n)
#else
    read(UNIT_OUTPUT) buffer
#endif

  end subroutine read_string

!===============================================================================
! WRITE_ATTRIBUTE_STRING
!===============================================================================

  subroutine write_attribute_string(var, attr_type, attr_str, group)

    character(*), intent(in)           :: var
    character(*), intent(in)           :: attr_type
    character(*), intent(in)           :: attr_str
    character(*), intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_attribute_string(temp_group, var, attr_type, attr_str)
#endif

  end subroutine write_attribute_string

!===============================================================================
! WRITE_TALLY_RESULT
!===============================================================================

  subroutine write_tally_result(buffer, name, group, n1, n2)

    character(*),      intent(in), optional :: group
    character(*),      intent(in)           :: name
    integer,           intent(in)           :: n1, n2
    type(TallyResult), intent(in), target   :: buffer(n1, n2)

#ifdef HDF5
    integer          :: hdf5_err
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr
#elif MPI
#else
    integer :: j,k
#endif

#ifdef HDF5

    ! Open up sub-group if present
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    end if

    ! Set overall size of vector to write
    dims(1) = n1*n2 

    ! Create up a dataspace for size
    call h5screate_simple_f(1, dims, dspace, hdf5_err)

    ! Create the dataset
    call h5dcreate_f(temp_group, name, hdf5_tallyresult_t, dspace, dset, &
         hdf5_err)

    ! Set pointer to first value and write
    f_ptr = c_loc(buffer(1,1))
    call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

    ! Close ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    if (present(group)) then
      call hdf5_close_group()
    end if

#elif MPI

    ! Write out tally buffer
    call MPI_FILE_WRITE(mpi_fh, buffer, n1*n2, MPI_TALLYRESULT, &
         MPI_STATUS_IGNORE, mpi_err)

#else

    ! Write out tally buffer
    do k = 1, n2
      do j = 1, n1
        write(UNIT_OUTPUT) buffer(j,k) % sum
        write(UNIT_OUTPUT) buffer(j,k) % sum_sq
      end do
    end do

#endif 
   
  end subroutine write_tally_result

!===============================================================================
! READ_TALLY_RESULT 
!===============================================================================

  subroutine read_tally_result(buffer, name, group, n1, n2)

    character(*),      intent(in), optional  :: group
    character(*),      intent(in)            :: name
    integer,           intent(in)            :: n1, n2
    type(TallyResult), intent(inout), target :: buffer(n1, n2)

#ifdef HDF5
    integer          :: hdf5_err
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr
#elif MPI
#else
    integer :: j,k
#endif

#ifdef HDF5

    ! Open up sub-group if present
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    end if

    ! Set overall size of vector to read
    dims(1) = n1*n2 

    ! Open the dataset
    call h5dopen_f(temp_group, name, dset, hdf5_err)

    ! Set pointer to first value and write
    f_ptr = c_loc(buffer(1,1))
    call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

    ! Close ids
    call h5dclose_f(dset, hdf5_err)
    if (present(group)) then
      call hdf5_close_group()
    end if

#elif MPI

    ! Write out tally buffer
    call MPI_FILE_READ(mpi_fh, buffer, n1*n2, MPI_TALLYRESULT, &
         MPI_STATUS_IGNORE, mpi_err)

#else

    ! Write out tally buffer
    do k = 1, n2
      do j = 1, n1
        read(UNIT_OUTPUT) buffer(j,k) % sum
        read(UNIT_OUTPUT) buffer(j,k) % sum_sq
      end do
    end do

#endif 
   
  end subroutine read_tally_result


!===============================================================================
! WRITE_SOURCE_BANK
!===============================================================================

  subroutine write_source_bank()

#ifdef HDF5
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dset
    integer(HID_T)   :: dspace
    integer(HID_T)   :: memspace
    integer(HID_T)   :: plist
    integer          :: rank
    integer(8)       :: offset(1)
    type(c_ptr)      :: f_ptr
#elif MPI
    integer(MPI_OFFSET_KIND) :: offset
    integer                  :: size_offset_kind
    integer                  :: size_bank
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims(1) = n_particles
    rank = 1

    ! Create that dataspace
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)

    ! Create the dataset for that dataspace
    call h5dcreate_f(hdf5_fh, "source_bank", hdf5_bank_t, dspace, dset, hdf5_err)

    ! Close the dataspace
    call h5sclose_f(dspace, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(rank, dims, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = bank_first - 1_8
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

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
    dims(1) = work

    ! Create dataspace
    call h5screate_simple_f(1, dims, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(hdf5_fh, "source_bank", hdf5_bank_t, &
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
    if (master) call MPI_FILE_GET_POSITION(mpi_fh, offset, mpi_err)

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
    offset = offset + size_bank*maxwork*rank

    ! Write all source sites
    call MPI_FILE_WRITE_AT(mpi_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpi_err)

#else

    ! Write out source sites
    write(UNIT_OUTPUT) source_bank

#endif

  end subroutine write_source_bank

!===============================================================================
! READ_SOURCE_BANK
!===============================================================================

  subroutine read_source_bank()

#ifdef HDF5
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dset
    integer(HID_T)   :: dspace
    integer(HID_T)   :: memspace
    integer(HID_T)   :: plist
    integer          :: rank
    integer(8)       :: offset(1)
    type(c_ptr)      :: f_ptr
#elif MPI
    integer(MPI_OFFSET_KIND) :: offset
    integer                  :: size_offset_kind
    integer                  :: size_bank
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims(1) = n_particles
    rank = 1

    ! Open the dataset
    call h5dopen_f(hdf5_fh, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(rank, dims, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = bank_first - 1_8
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

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

    ! Set size
    dims(1) = work

    ! Open dataset
    call h5dopen_f(hdf5_fh, "source_bank", dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Read dataset from file
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)

# endif 

#elif MPI

    ! Get current offset for master 
    if (master) call MPI_FILE_GET_POSITION(mpi_fh, offset, mpi_err)

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
    offset = offset + size_bank*maxwork*rank

    ! Write all source sites
    call MPI_FILE_READ_AT(mpi_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpi_err)

#else

    ! Write out source sites
    read(UNIT_OUTPUT) source_bank

#endif

  end subroutine read_source_bank

end module output_interface
