module hdf5_interface

#ifdef HDF5

  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING

#ifdef MPI
  use mpi
#endif

  implicit none

  integer(HID_T)  :: hdf5_fh
  integer(HID_T)  :: temp_group 
  integer         :: hdf5_err

  interface hdf5_write_data
    module procedure hdf5_write_double
    module procedure hdf5_write_double_1Darray
    module procedure hdf5_write_integer
    module procedure hdf5_write_integer_1Darray
    module procedure hdf5_write_long
    module procedure hdf5_write_string
  end interface hdf5_write_data

  interface hdf5_read_data
    module procedure hdf5_read_double
    module procedure hdf5_read_double_1Darray
    module procedure hdf5_read_integer
    module procedure hdf5_read_integer_1Darray
    module procedure hdf5_read_long
    module procedure hdf5_read_string
  end interface hdf5_read_data

contains

!===============================================================================
! HDF5_FILE_CREATE
!===============================================================================

  subroutine hdf5_file_create(filename, file_id)

    character(*) :: filename
    integer(HID_T)          :: file_id

    ! Create the file
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)

  end subroutine hdf5_file_create

!===============================================================================
! HDF5_FILE_OPEN
!===============================================================================

  subroutine hdf5_file_open(filename, file_id, mode)

    character(*) :: filename
    character(*)            :: mode
    integer(HID_T)          :: file_id
    integer                 :: open_mode

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'rw') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Open file
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)

  end subroutine hdf5_file_open

!===============================================================================
! HDF5_FILE_CLOSE
!===============================================================================

  subroutine hdf5_file_close(file_id)

    integer(HID_T) :: file_id

    ! Close the file
    call h5fclose_f(file_id, hdf5_err)

  end subroutine hdf5_file_close

#ifdef MPI

!===============================================================================
! HDF5_PARALLEL_FILE_CREATE
!===============================================================================

  subroutine hdf5_parallel_file_create(filename, file_id)

    character(*) :: filename
    integer(HID_T)          :: file_id
    integer(HID_T)          :: plist_id

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_err)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Create the file collectively
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                     access_prp = plist_id)

    ! Close the property list
    call h5pclose_f(plist_id, hdf5_err)

  end subroutine hdf5_parallel_file_create

!===============================================================================
! HDF5_PARALLEL_FILE_OPEN
!===============================================================================

  subroutine hdf5_parallel_file_open(filename, file_id, mode)

    character(*) :: filename
    character(*)            :: mode
    integer(HID_T)          :: file_id
    integer(HID_T)          :: plist_id
    integer                 :: open_mode

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_err)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'rw') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Create the file collectively
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err, &
                   access_prp = plist_id)

    ! Close the property list
    call h5pclose_f(plist_id, hdf5_err)

  end subroutine hdf5_parallel_file_open

#endif

!===============================================================================
! HDF5_OPEN_GROUP
!===============================================================================

  subroutine hdf5_open_group(group)

    character(*) :: group

    logical :: status

    ! Check if group exists
    call h5ltpath_valid_f(hdf5_fh, trim(group), .true., status, hdf5_err) 

    ! Either create or open group
    if (status) then
      call h5gopen_f(hdf5_fh, trim(group), temp_group, hdf5_err)
    else
      call h5gcreate_f(hdf5_fh, trim(group), temp_group, hdf5_err)
    end if

  end subroutine hdf5_open_group

!===============================================================================
! HDF5_CLOSE_GROUP
!===============================================================================

  subroutine hdf5_close_group()

    call h5gclose_f(temp_group, hdf5_err)

  end subroutine hdf5_close_group

!===============================================================================
! HDF5_WRITE_INTEGER
!===============================================================================

  subroutine hdf5_write_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_integer

!===============================================================================
! HDF5_WRITE_INTEGER_1DARRAY
!===============================================================================

  subroutine hdf5_write_integer_1Darray(group, name, buffer, len)

    integer,        intent(in) :: len
    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    integer,        intent(in) :: buffer(:)

    integer          :: rank
    integer(HSIZE_T) :: dims(1)

    rank = 1
    dims(1) = len

    call h5ltmake_dataset_int_f(group, name, rank, dims, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_1Darray

!===============================================================================
! HDF5_WRITE_LONG
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in) :: group
    character(*),       intent(in) :: name
    integer(8), target, intent(in) :: buffer
    integer(HID_T),     intent(in) :: long_type 

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr

    ! Create dataspace and dataset
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)
    call h5dcreate_f(group, name, long_type, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, long_type, f_ptr, hdf5_err)

    ! Close dataspace and dataset for long integer
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_long

!===============================================================================
! HDF5_WRITE_DOUBLE
!===============================================================================

  subroutine hdf5_write_double(group, name, buffer)

    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_double

!===============================================================================
! HDF5_WRITE_DOUBLE_1DARRAY
!===============================================================================

  subroutine hdf5_write_double_1Darray(group, name, buffer, len)

    integer,        intent(in) :: len
    integer(HID_T), intent(in) :: group
    character(*),   intent(in) :: name
    real(8),        intent(in) :: buffer(:)

    integer          :: rank
    integer(HSIZE_T) :: dims(1)

    rank = 1
    dims(1) = len

    call h5ltmake_dataset_double_f(group, name, rank, dims, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_1Darray

!===============================================================================
! HDF5_WRITE_STRING
!===============================================================================

  subroutine hdf5_write_string(group, name, buffer)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    character(*),   intent(in)    :: buffer

    call h5ltmake_dataset_string_f(group, name, buffer, hdf5_err)

  end subroutine hdf5_write_string
      
!===============================================================================
! HDF5_READ_INTEGER
!===============================================================================

  subroutine hdf5_read_integer(group, name, buffer)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    integer,        intent(inout) :: buffer

    integer          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_int_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_integer

!===============================================================================
! HDF5_READ_INTEGER_1DARRAY
!===============================================================================

  subroutine hdf5_read_integer_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    integer,        intent(inout) :: buffer(:)
    integer,        intent(in)    :: length

    integer(HSIZE_T) :: dims(1)

    dims(1) = length
    call h5ltread_dataset_int_f(group, name, buffer, dims, hdf5_err)

  end subroutine hdf5_read_integer_1Darray


!===============================================================================
! HDF5_READ_LONG
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in)  :: group
    character(*),       intent(in)  :: name
    integer(8), target, intent(out) :: buffer
    integer(HID_T),     intent(in)  :: long_type

    integer(HID_T) :: dset
    type(c_ptr)    :: f_ptr

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get pointer to buffer
    f_ptr = c_loc(buffer)

    ! Read data from dataset
    call h5dread_f(dset, long_type, f_ptr, hdf5_err)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)

  end subroutine hdf5_read_long

!===============================================================================
! HDF5_READ_DOUBLE
!===============================================================================

  subroutine hdf5_read_double(group, name, buffer)

    integer(HID_T), intent(in)  :: group
    character(*),   intent(in)  :: name
    real(8),        intent(out) :: buffer

    real(8)          :: buffer_copy(1)
    integer(HSIZE_T) :: dims(1) = (/1/)

    call h5ltread_dataset_double_f(group, name, buffer_copy, dims, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_double

!===============================================================================
! HDF5_READ_DOUBLE_1DARRAY
!===============================================================================

  subroutine hdf5_read_double_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)  :: group
    integer,        intent(in)  :: length
    character(*),   intent(in)  :: name
    real(8),        intent(out) :: buffer(:)

    integer(HSIZE_T) :: dims(1)

    dims(1) = length
    call h5ltread_dataset_double_f(group, name, buffer, dims, hdf5_err)

  end subroutine hdf5_read_double_1Darray

!===============================================================================
! HDF5_READ_STRING
!===============================================================================

  subroutine hdf5_read_string(group, name, buffer)

    integer(HID_T), intent(in)    :: group
    character(*),   intent(in)    :: name
    character(*),   intent(inout) :: buffer

    call h5ltread_dataset_string_f(group, name, buffer, hdf5_err)

  end subroutine hdf5_read_string

#endif

end module hdf5_interface
