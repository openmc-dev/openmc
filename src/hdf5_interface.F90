module hdf5_interface

#ifdef HDF5

  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING
  
#ifdef MPI
  use mpi
#endif

  implicit none

  integer          :: hdf5_err   ! HDF5 error code
  integer          :: hdf5_rank  ! rank of data
  integer(HID_T)   :: dset       ! data set handle
  integer(HID_T)   :: dspace     ! data or file space handle
  integer(HID_T)   :: hdf5_fh    ! HDF5 file handle
  integer(HID_T)   :: temp_group ! temporary HDF5 group handle
  integer(HID_T)   :: memspace   ! data space handle for individual procs
  integer(HID_T)   :: plist      ! property list handle
  integer(HSIZE_T) :: dims1(1)   ! dims type for 1-D array
  integer(HSIZE_T) :: dims2(2)   ! dims type for 2-D array
  integer(HSIZE_T) :: dims3(3)   ! dims type for 3-D array
  type(c_ptr)      :: f_ptr      ! pointer to data

  ! Generic HDF5 write procedure interface
  interface hdf5_write_data
    module procedure hdf5_write_double
    module procedure hdf5_write_double_1Darray
    module procedure hdf5_write_integer
    module procedure hdf5_write_integer_1Darray
    module procedure hdf5_write_long
    module procedure hdf5_write_string
  end interface hdf5_write_data

  ! Generic HDF5 read procedure interface
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
! HDF5_FILE_CREATE creates HDF5 file
!===============================================================================

  subroutine hdf5_file_create(filename, file_id)

    character(*),   intent(in)    :: filename ! name of file
    integer(HID_T), intent(inout) :: file_id  ! file handle

    ! Create the file
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)

  end subroutine hdf5_file_create

!===============================================================================
! HDF5_FILE_OPEN opens HDF5 file
!===============================================================================

  subroutine hdf5_file_open(filename, file_id, mode)

    character(*),  intent(in)      :: filename ! name of file
    character(*),  intent(in)      :: mode     ! access mode to file
    integer(HID_T), intent(inout)  :: file_id  ! file handle

    integer :: open_mode ! HDF5 open mode

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'w') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Open file
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)

  end subroutine hdf5_file_open

!===============================================================================
! HDF5_FILE_CLOSE closes HDF5 file
!===============================================================================

  subroutine hdf5_file_close(file_id)

    integer(HID_T), intent(inout) :: file_id ! file handle

    ! Close the file
    call h5fclose_f(file_id, hdf5_err)

  end subroutine hdf5_file_close

#ifdef MPI

!===============================================================================
! HDF5_PARALLEL_FILE_CREATE creates HDF5 file with parallel I/O
!===============================================================================

  subroutine hdf5_parallel_file_create(filename, file_id)

    character(*),   intent(in)    :: filename ! name of file
    integer(HID_T), intent(inout) :: file_id  ! file handle

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
    call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Create the file collectively
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                     access_prp = plist)

    ! Close the property list
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_file_create

!===============================================================================
! HDF5_PARALLEL_FILE_OPEN opens HDF5 file with parallel I/O
!===============================================================================

  subroutine hdf5_parallel_file_open(filename, file_id, mode)

    character(*),  intent(in)     :: filename ! name of file
    character(*),  intent(in)     :: mode     ! access mode
    integer(HID_T), intent(inout) :: file_id  ! file handle

    integer        :: open_mode ! HDF5 access mode

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
    call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'w') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Create the file collectively
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err, &
                   access_prp = plist)

    ! Close the property list
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_file_open

#endif

!===============================================================================
! HDF5_OPEN_GROUP creates/opens HDF5 group to temp_group
!===============================================================================

  subroutine hdf5_open_group(group)

    character(*), intent(in) :: group ! name of group

    logical :: status ! does the group exist

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
! HDF5_CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine hdf5_close_group()

    ! Close the group
    call h5gclose_f(temp_group, hdf5_err)

  end subroutine hdf5_close_group

!===============================================================================
! HDF5_WRITE_INTEGER writes integer scalar data
!===============================================================================

  subroutine hdf5_write_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group  ! name of group
    character(*),   intent(in) :: name   ! name of data
    integer,        intent(in) :: buffer ! data to write

    ! Set rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_integer

!===============================================================================
! HDF5_WRITE_INTEGER_1DARRAY writes integer 1-D array
!===============================================================================

  subroutine hdf5_write_integer_1Darray(group, name, buffer, len)

    integer,        intent(in) :: len       ! length of array to write
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(:) ! data to write

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = len

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_1Darray

!===============================================================================
! HDF5_WRITE_INTEGER_2DARRAY write integer 2-D array
!===============================================================================

  subroutine hdf5_write_integer_2Darray(group, name, buffer, length)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(length(1),length(2)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 2
    dims2 = length

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims2, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_2Darray

!===============================================================================
! HDF5_WRITE_INTEGER_3DARRAY writes integer 3-D array
!===============================================================================

  subroutine hdf5_write_integer_3Darray(group, name, buffer, length)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(length(1),length(2), length(3)) ! data

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims3, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_3Darray

!===============================================================================
! HDF5_WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in) :: group     ! name of group
    character(*),       intent(in) :: name      ! name of data
    integer(8), target, intent(in) :: buffer    ! data to write
    integer(HID_T),     intent(in) :: long_type ! HDF5 long type

    ! Set up rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    ! Create dataspace and dataset
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)
    call h5dcreate_f(group, name, long_type, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, long_type, f_ptr, hdf5_err)

    ! Close dataspace and dataset for long integer
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_long

!===============================================================================
! HDF5_WRITE_DOUBLE writes double precision scalar data
!===============================================================================

  subroutine hdf5_write_double(group, name, buffer)

    integer(HID_T), intent(in) :: group  ! name of group
    character(*),   intent(in) :: name   ! name of data
    real(8),        intent(in) :: buffer ! data to write

    ! Set up rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1 

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_double

!===============================================================================
! HDF5_WRITE_DOUBLE_1DARRAY writes double precision 1-D array
!===============================================================================

  subroutine hdf5_write_double_1Darray(group, name, buffer, len)

    integer,        intent(in) :: len       ! length of array
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(:) ! data to write

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = len

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_1Darray

!===============================================================================
! HDF5_WRITE_DOUBLE_2DARRAY writes double precision 2-D array
!===============================================================================

  subroutine hdf5_write_double_2Darray(group, name, buffer, length)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(length(1),length(2)) ! data to write

    ! Set rank and dimensions of data
    hdf5_rank = 2
    dims2 = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims2, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_2Darray

!===============================================================================
! HDF5_WRITE_DOUBLE_3DARRAY writes double precision 3-D aray
!===============================================================================

  subroutine hdf5_write_double_3Darray(group, name, buffer, length)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(length(1),length(2), length(3)) ! data

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims3, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_3Darray

!===============================================================================
! HDF5_WRITE_STRING writes string data
!===============================================================================

  subroutine hdf5_write_string(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    character(*),   intent(in)    :: buffer ! data to write
    integer,        intent(in)    :: length

    character(len=length), dimension(1) :: str_tmp

!   Fortran 2003 implementation not compatible with IBM compiler Feb 2013
!   type(c_ptr), dimension(1), target :: wdata
!   character(len=length, kind=c_char), dimension(1), target :: c_str
!   dims1(1) = 1
!   call h5screate_simple_f(1, dims1, dspace, hdf5_err)
!   call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)
!   c_str(1) = buffer
!   wdata(1) = c_loc(c_str(1))
!   f_ptr = c_loc(wdata(1))

    ! Number of strings to write
    dims1(1) = 1

    ! Insert null character at end of string when writing
    call h5tset_strpad_f(H5T_STRING, H5T_STR_NULLPAD_F, hdf5_err)

    ! Create the dataspace and dataset
    call h5screate_simple_f(1, dims1, dspace, hdf5_err)
    call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)

    ! Set up dimesnions of string to write
    dims2 = (/length, 1/) ! full array of strings to write 
    dims1(1) = length     ! length of string

    ! Copy over string buffer to a rank 1 array
    str_tmp(1) = buffer

    ! Write the variable dataset
    call h5dwrite_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace)

    ! Close all
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_string

!===============================================================================
! HDF5_WRITE_ATTRIBUTE_STRING writes a string attribute to a variables
!===============================================================================

  subroutine hdf5_write_attribute_string(group, var, attr_type, attr_str)

    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: var       ! name of varaible to set attr
    character(*),   intent(in) :: attr_type ! the attr type id
    character(*),   intent(in) :: attr_str  ! attribute sting

    call h5ltset_attribute_string_f(group, var, attr_type, attr_str, hdf5_err)

  end subroutine hdf5_write_attribute_string

!===============================================================================
! HDF5_READ_INTEGER reads integer scalar data
!===============================================================================

  subroutine hdf5_read_integer(group, name, buffer)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    integer,        intent(inout) :: buffer ! read data to here 

    integer :: buffer_copy(1) ! need an array for read

    ! Set up dimensions
    dims1(1) = 1

    ! Read data
    call h5ltread_dataset_int_f(group, name, buffer_copy, dims1, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_integer

!===============================================================================
! HDF5_READ_INTEGER_1DARRAY reads integer 1-D array
!===============================================================================

  subroutine hdf5_read_integer_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(inout) :: buffer(:) ! read data to here
    integer,        intent(in)    :: length    ! length of array

    ! Set dimensions
    dims1(1) = length

    ! Read data
    call h5ltread_dataset_int_f(group, name, buffer, dims1, hdf5_err)

  end subroutine hdf5_read_integer_1Darray


!===============================================================================
! HDF5_READ_LONG read long integer scalar data
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in)  :: group     ! name of group
    character(*),       intent(in)  :: name      ! name of data
    integer(8), target, intent(out) :: buffer    ! read data to here
    integer(HID_T),     intent(in)  :: long_type ! long integer type

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
! HDF5_READ_DOUBLE reads double precision scalar data
!===============================================================================

  subroutine hdf5_read_double(group, name, buffer)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    real(8),        intent(inout) :: buffer ! read data to here

    real(8) :: buffer_copy(1) ! need an array for read

    ! Set up dimensions
    dims1(1) = 1

    ! Read data
    call h5ltread_dataset_double_f(group, name, buffer_copy, dims1, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_double

!===============================================================================
! HDF5_READ_DOUBLE_1DARRAY reads double precision 1-D array
!===============================================================================

  subroutine hdf5_read_double_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group     ! name of group
    integer,        intent(in)    :: length    ! length of array
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(inout) :: buffer(:) ! read data to here

    ! Set dimensions of data
    dims1(1) = length

    ! Read data
    call h5ltread_dataset_double_f(group, name, buffer, dims1, hdf5_err)

  end subroutine hdf5_read_double_1Darray

!===============================================================================
! HDF5_READ_STRING reads string data
!===============================================================================

  subroutine hdf5_read_string(group, name, buffer)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    character(*),   intent(inout) :: buffer ! read data to here

    call h5ltread_dataset_string_f(group, name, buffer, hdf5_err)

  end subroutine hdf5_read_string

# ifdef MPI
!===============================================================================
! HDF5_PARALLEL_READ_INTEGER reads interger scalar data in parallel
!===============================================================================

  subroutine hdf5_parallel_read_integer(group, name, buffer, p_type)

    integer(HID_T),  intent(in)    :: group  ! name of group
    character(*),    intent(in)    :: name   ! name of data
    integer, target, intent(inout) :: buffer ! read data to here
    integer,         intent(in)    :: p_type ! independent or collective I/O

    ! Set up dimension
    dims1(1) = 1

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_integer

!===============================================================================
! HDF5_PARALLEL_READ_INTEGER_1DARRAY reads integer 1-D array in parallel
!===============================================================================

  subroutine hdf5_parallel_read_integer_1Darray(group, name, buffer, length, &
             p_type)

    integer(HID_T),  intent(in)    :: group     ! name of group
    character(*),    intent(in)    :: name      ! name of data
    integer,         intent(in)    :: length    ! length of array
    integer, target, intent(inout) :: buffer(length) ! read data to here
    integer,         intent(in)    :: p_type    ! independent or collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer(1))
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_integer_1Darray


!===============================================================================
! HDF5_PARALLEL_READ_LONG read long integer scalar data in parallel
!===============================================================================

  subroutine hdf5_parallel_read_long(group, name, buffer, long_type, p_type)

    integer(HID_T),     intent(in)  :: group     ! name of group
    character(*),       intent(in)  :: name      ! name of data
    integer(8), target, intent(out) :: buffer    ! read data to here
    integer(HID_T),     intent(in)  :: long_type ! long integer type
    integer,            intent(in)  :: p_type    ! independent or collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, long_type, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_long

!===============================================================================
! HDF5_PARALLEL_READ_DOUBLE reads double precision scalar data in parallel
!===============================================================================

  subroutine hdf5_parallel_read_double(group, name, buffer, p_type)

    integer(HID_T),  intent(in)    :: group  ! name of group
    character(*),    intent(in)    :: name   ! name of data
    real(8), target, intent(inout) :: buffer ! read data to here
    integer,         intent(in)    :: p_type ! independent or collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_double

!===============================================================================
! HDF5_PARLLEL_READ_DOUBLE_1DARRAY reads double precision 1-D array in parallel
!===============================================================================

  subroutine hdf5_parallel_read_double_1Darray(group, name, buffer, length, &
             p_type)

    integer(HID_T),  intent(in)    :: group     ! name of group
    integer,         intent(in)    :: length    ! length of array
    integer,         intent(in)    :: p_type    ! indepedent or collective I/O
    character(*),    intent(in)    :: name      ! name of data
    real(8), target, intent(inout) :: buffer(length) ! read data to here

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer(1))
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_double_1Darray

!===============================================================================
! HDF5_PARALLEL_READ_STRING reads string data
!===============================================================================

  subroutine hdf5_parallel_read_string(group, name, buffer, length, p_type)

    integer(HID_T),       intent(in)    :: group  ! name of group
    character(*),         intent(in)    :: name   ! name of data
    character(*), target, intent(inout) :: buffer ! read data to here
    integer,              intent(in)    :: p_type ! independent or collective IO
    integer,              intent(in)    :: length ! length of string

    character(len=length), dimension(1) :: str_tmp
    ! Fortran 2003 implementation not compatible with IBM Feb 2013 compiler
!    type(c_ptr), dimension(1), target :: buf_ptr
!    character(len=length, kind=c_char), pointer :: chr_ptr
!    f_ptr = c_loc(buf_ptr(1))
!    call h5dread_f(dset, H5T_STRING, f_ptr, hdf5_err, xfer_prp=plist)
!    call c_f_pointer(buf_ptr(1), chr_ptr) 
!    buffer = chr_ptr
!    nullify(chr_ptr)

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    call h5pset_dxpl_mpio_f(plist, p_type, hdf5_err)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get dataspace to read
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Set dimensions
    dims2 = (/length, 1/)
    dims1(1) = length

    ! Read in the data
    call h5dread_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace, xfer_prp = plist)

    ! Copy over buffer
    buffer = str_tmp(1)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_parallel_read_string

# endif

#endif

end module hdf5_interface
