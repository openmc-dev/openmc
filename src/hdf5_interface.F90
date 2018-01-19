module hdf5_interface

!==============================================================================
! HDF5_INTERFACE -- This module provides the high-level procedures which greatly
! simplify writing/reading different types of data to HDF5 files. In order to
! get it to work with gfotran 4.6, all the write_<type>_ND subroutines had to be
! split into two procedures, one accepting an assumed-shape array and another
! one with an explicit-shape array since in gfortran 4.6 C_LOC does not work
! with an assumed-shape array. When we move to gfortran 4.9+, these procedures
! can be combined into one simply accepting an assumed-shape array.
!==============================================================================

  use, intrinsic :: ISO_C_BINDING

  use hdf5
  use h5lt

  use error, only: fatal_error
#ifdef PHDF5
  use message_passing, only: mpi_intracomm, MPI_INFO_NULL
#endif

  implicit none
  private

  integer(HID_T), public :: hdf5_bank_t        ! Compound type for Bank
  integer(HID_T), public :: hdf5_integer8_t    ! type for integer(8)

  interface write_dataset
    module procedure write_double
    module procedure write_double_1D
    module procedure write_double_2D
    module procedure write_double_3D
    module procedure write_double_4D
    module procedure write_integer
    module procedure write_integer_1D
    module procedure write_integer_2D
    module procedure write_integer_3D
    module procedure write_integer_4D
    module procedure write_long
    module procedure write_string
    module procedure write_string_1D
  end interface write_dataset

  interface read_dataset
    module procedure read_double
    module procedure read_double_1D
    module procedure read_double_2D
    module procedure read_double_3D
    module procedure read_double_4D
    module procedure read_integer
    module procedure read_integer_1D
    module procedure read_integer_2D
    module procedure read_integer_3D
    module procedure read_integer_4D
    module procedure read_long
    module procedure read_string
    module procedure read_string_1D
    module procedure read_complex_2D
  end interface read_dataset

  interface read_attribute
    module procedure read_attribute_double
    module procedure read_attribute_double_1D
    module procedure read_attribute_double_2D
    module procedure read_attribute_integer
    module procedure read_attribute_integer_1D
    module procedure read_attribute_integer_2D
    module procedure read_attribute_string
    module procedure read_attribute_string_1D
    module procedure read_attribute_logical
  end interface read_attribute

  interface write_attribute
    module procedure write_attribute_double
    module procedure write_attribute_double_1D
    module procedure write_attribute_integer
    module procedure write_attribute_integer_1D
    module procedure write_attribute_string
  end interface write_attribute

  public :: write_dataset
  public :: read_dataset
  public :: attribute_exists
  public :: write_attribute
  public :: read_attribute
  public :: file_create
  public :: file_open
  public :: file_close
  public :: create_group
  public :: object_exists
  public :: open_group
  public :: close_group
  public :: open_dataset
  public :: close_dataset
  public :: get_shape
  public :: get_ndims
  public :: get_groups
  public :: get_datasets
  public :: get_name

contains

!===============================================================================
! FILE_CREATE creates HDF5 file
!===============================================================================

  function file_create(filename, parallel) result(file_id)
    character(*),      intent(in)    :: filename ! name of file
    logical, optional, intent(in)    :: parallel ! whether to write in serial
    integer(HID_T) :: file_id

    integer(HID_T) :: plist      ! property list handle
    integer        :: hdf5_err   ! HDF5 error code
    logical :: parallel_

    ! Check for serial option
    parallel_ = .false.
#ifdef PHDF5
    if (present(parallel)) parallel_ = parallel
#endif

    if (parallel_) then
      ! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
#ifdef PHDF5
#ifdef OPENMC_MPIF08
      call h5pset_fapl_mpio_f(plist, mpi_intracomm%MPI_VAL, &
           MPI_INFO_NULL%MPI_VAL, hdf5_err)
#else
      call h5pset_fapl_mpio_f(plist, mpi_intracomm, MPI_INFO_NULL, hdf5_err)
#endif
#endif

      ! Create the file collectively
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                       access_prp = plist)

      ! Close the property list
      call h5pclose_f(plist, hdf5_err)
    else
      ! Create the file
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)
    end if

  end function file_create

!===============================================================================
! FILE_OPEN opens HDF5 file
!===============================================================================

  function file_open(filename, mode, parallel) result(file_id)
    character(*),      intent(in)    :: filename ! name of file
    character(*),      intent(in)    :: mode     ! access mode to file
    logical, optional, intent(in)    :: parallel ! whether to write in serial
    integer(HID_T) :: file_id

    logical :: parallel_
    integer(HID_T) :: plist     ! property list handle
    integer        :: hdf5_err  ! HDF5 error code
    integer        :: open_mode ! HDF5 open mode

    ! Check for serial option
    parallel_ = .false.
#ifdef PHDF5
    if (present(parallel)) parallel_ = parallel
#endif

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (mode == 'w') open_mode = H5F_ACC_RDWR_F

    if (parallel_) then
      ! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
#ifdef PHDF5
#ifdef OPENMC_MPIF08
      call h5pset_fapl_mpio_f(plist, mpi_intracomm%MPI_VAL, &
           MPI_INFO_NULL%MPI_VAL, hdf5_err)
#else
      call h5pset_fapl_mpio_f(plist, mpi_intracomm, MPI_INFO_NULL, hdf5_err)
#endif
#endif

      ! Open the file collectively
      call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err, &
                     access_prp = plist)

      ! Close the property list
      call h5pclose_f(plist, hdf5_err)
    else
      ! Open file
      call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)
    end if

  end function file_open

!===============================================================================
! FILE_CLOSE closes HDF5 file
!===============================================================================

  subroutine file_close(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: hdf5_err

    call h5fclose_f(file_id, hdf5_err)
  end subroutine file_close

!===============================================================================
! GET_GROUPS Gets a list of all the groups in a given location.
!===============================================================================

  subroutine get_groups(object_id, names)
    integer(HID_T), intent(in)  :: object_id
    character(len=150), allocatable, intent(out) :: names(:)

    integer :: n_members, i, group_count, type
    integer :: hdf5_err
    character(len=150) :: name

    ! Get number of members in this location
    call h5gn_members_f(object_id, './', n_members, hdf5_err)

    ! Get the number of groups
    group_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_GROUP_F) then
        group_count = group_count + 1
      end if
    end do

    ! Now we can allocate the storage for the ids
    allocate(names(group_count))
    group_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_GROUP_F) then
        group_count = group_count + 1
        names(group_count) = trim(name)
      end if
    end do

  end subroutine get_groups

!===============================================================================
! CHECK_ATTRIBUTE Checks to see if an attribute exists in the object
!===============================================================================

  function attribute_exists(object_id, name) result(exists)
    integer(HID_T), intent(in) :: object_id
    character(*),   intent(in) :: name ! name of group
    logical :: exists

    integer :: hdf5_err ! HDF5 error code

    ! Check if attribute exists
    call h5aexists_by_name_f(object_id, '.', trim(name), exists, hdf5_err)

  end function attribute_exists

!===============================================================================
! CHECK_GROUP Checks to see if a group exists in the object
!===============================================================================

  function object_exists(object_id, name) result(exists)
    integer(HID_T), intent(in) :: object_id
    character(*),   intent(in) :: name ! name of group
    logical :: exists

    integer :: hdf5_err ! HDF5 error code

    ! Check if group exists
    call h5ltpath_valid_f(object_id, trim(name), .true., exists, hdf5_err)

  end function object_exists

!===============================================================================
! GET_DATASETS Gets a list of all the datasets in a given location.
!===============================================================================

  subroutine get_datasets(object_id, names)
    integer(HID_T), intent(in)  :: object_id
    character(len=150), allocatable, intent(out) :: names(:)

    integer :: n_members, i, dset_count, type
    integer :: hdf5_err
    character(len=150) :: name


    ! Get number of members in this location
    call h5gn_members_f(object_id, './', n_members, hdf5_err)

    ! Get the number of datasets
    dset_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_DATASET_F ) then
        dset_count = dset_count + 1
      end if
    end do

    ! Now we can allocate the storage for the ids
    allocate(names(dset_count))
    dset_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_DATASET_F ) then
        dset_count = dset_count + 1
        names(dset_count) = trim(name)
      end if
    end do

  end subroutine get_datasets

!===============================================================================
! GET_NAME Obtains the name of the current group in group_id
!===============================================================================

  function get_name(group_id, name_len_) result(name)
    integer(HID_T),            intent(in) :: group_id
    integer(SIZE_T), optional, intent(in) :: name_len_

    character(len=150) :: name ! name of group
    integer(SIZE_T) :: name_len, name_file_len
    integer :: hdf5_err ! HDF5 error code

    if (present(name_len_)) then
      name_len = name_len_
    else
      name_len = 150
    end if

    call h5iget_name_f(group_id, name, name_len, name_file_len, hdf5_err)
  end function get_name

!===============================================================================
! OPEN_GROUP opens an existing HDF5 group
!===============================================================================

  function open_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    logical :: exists   ! does the group exist
    integer :: hdf5_err ! HDF5 error code

    ! Check if group exists
    exists = object_exists(group_id, name)

    ! open group if it exists
    if (exists) then
      call h5gopen_f(group_id, trim(name), newgroup_id, hdf5_err)
    else
      call fatal_error("The group '" // trim(name) // "' does not exist.")
    end if
  end function open_group

!===============================================================================
! CREATE_GROUP creates a new HDF5 group
!===============================================================================

  function create_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    integer :: hdf5_err ! HDF5 error code
    logical :: exists   ! does the group exist

    ! Check if group exists
    exists = object_exists(group_id, name)

    ! create group
    if (exists) then
      call fatal_error("The group '" // trim(name) // "' already exists.")
    else
      call h5gcreate_f(group_id, trim(name), newgroup_id, hdf5_err)
    end if
  end function create_group

!===============================================================================
! CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine close_group(group_id)
    integer(HID_T), intent(inout) :: group_id

    integer :: hdf5_err ! HDF5 error code

    call h5gclose_f(group_id, hdf5_err)
    if (hdf5_err < 0) then
      call fatal_error("Unable to close HDF5 group.")
    end if
  end subroutine close_group

!===============================================================================
! OPEN_DATASET opens an existing HDF5 dataset
!===============================================================================

  function open_dataset(group_id, name) result(dataset_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of dataset
    integer(HID_T)             :: dataset_id

    logical :: exists   ! does the dataset exist
    integer :: hdf5_err ! HDF5 error code

    ! Check if group exists
    exists = object_exists(group_id, name)

    ! open group if it exists
    if (exists) then
      call h5dopen_f(group_id, trim(name), dataset_id, hdf5_err)
    else
      call fatal_error("The dataset '" // trim(name) // "' does not exist.")
    end if
  end function open_dataset

!===============================================================================
! CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine close_dataset(dataset_id)
    integer(HID_T), intent(inout) :: dataset_id

    integer :: hdf5_err ! HDF5 error code

    call h5dclose_f(dataset_id, hdf5_err)
    if (hdf5_err < 0) then
      call fatal_error("Unable to close HDF5 dataset.")
    end if
  end subroutine close_dataset

!===============================================================================
! WRITE_DOUBLE writes double precision scalar data
!===============================================================================

  subroutine write_double(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    real(8),      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up independentive vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double

!===============================================================================
! READ_DOUBLE reads double precision scalar data
!===============================================================================

  subroutine read_double(buffer, obj_id, name, indep)
    real(8),      target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer        :: hdf5_err
    integer        :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset_id
    type(c_ptr)    :: f_ptr

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_double

!===============================================================================
! WRITE_DOUBLE_1D writes double precision 1-D array data
!===============================================================================

  subroutine write_double_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_1D

  subroutine write_double_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_1D_explicit

!===============================================================================
! READ_DOUBLE_1D reads double precision 1-D array data
!===============================================================================

  subroutine read_double_1D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_double_1D_explicit(dset_id, dims, buffer, indep)
    else
      call read_double_1D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_double_1D

  subroutine read_double_1D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(1)
    real(8), target,   intent(inout) :: buffer(dims(1))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if
  end subroutine read_double_1D_explicit

!===============================================================================
! WRITE_DOUBLE_2D writes double precision 2-D array data
!===============================================================================

  subroutine write_double_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_2D

  subroutine write_double_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in)       :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(2, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_2D_explicit

!===============================================================================
! READ_DOUBLE_2D reads double precision 2-D array data
!===============================================================================

  subroutine read_double_2D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(2)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_double_2D_explicit(dset_id, dims, buffer, indep)
    else
      call read_double_2D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_double_2D

  subroutine read_double_2D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(2)
    real(8), target,   intent(inout) :: buffer(dims(1),dims(2))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if
  end subroutine read_double_2D_explicit

!===============================================================================
! WRITE_DOUBLE_3D writes double precision 3-D array data
!===============================================================================

  subroutine write_double_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_3D

  subroutine write_double_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(3, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_3D_explicit

!===============================================================================
! READ_DOUBLE_3D reads double precision 3-D array data
!===============================================================================

  subroutine read_double_3D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(3)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_double_3D_explicit(dset_id, dims, buffer, indep)
    else
      call read_double_3D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_double_3D

  subroutine read_double_3D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(3)
    real(8), target,   intent(inout) :: buffer(dims(1),dims(2),dims(3))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if
  end subroutine read_double_3D_explicit

!===============================================================================
! WRITE_DOUBLE_4D writes double precision 4-D array data
!===============================================================================

  subroutine write_double_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_4D

  subroutine write_double_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(4, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_4D_explicit

!===============================================================================
! READ_DOUBLE_4D reads double precision 4-D array data
!===============================================================================

  subroutine read_double_4D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(4)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_double_4D_explicit(dset_id, dims, buffer, indep)
    else
      call read_double_4D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_double_4D

  subroutine read_double_4D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(4)
    real(8), target,   intent(inout) :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if
  end subroutine read_double_4D_explicit

!===============================================================================
! WRITE_INTEGER writes integer precision scalar data
!===============================================================================

  subroutine write_integer(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer,      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer

!===============================================================================
! READ_INTEGER reads integer precision scalar data
!===============================================================================

  subroutine read_integer(buffer, obj_id, name, indep)
    integer,      target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer        :: hdf5_err
    integer        :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset_id
    type(c_ptr)    :: f_ptr

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_integer

!===============================================================================
! WRITE_INTEGER_1D writes integer precision 1-D array data
!===============================================================================

  subroutine write_integer_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_1D

  subroutine write_integer_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_1D_explicit

!===============================================================================
! READ_INTEGER_1D reads integer precision 1-D array data
!===============================================================================

  subroutine read_integer_1D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_integer_1D_explicit(dset_id, dims, buffer, indep)
    else
      call read_integer_1D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_integer_1D

  subroutine read_integer_1D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(1)
    integer, target,   intent(inout) :: buffer(dims(1))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if
  end subroutine read_integer_1D_explicit

!===============================================================================
! WRITE_INTEGER_2D writes integer precision 2-D array data
!===============================================================================

  subroutine write_integer_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_2D

  subroutine write_integer_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(2, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_2D_explicit

!===============================================================================
! READ_INTEGER_2D reads integer precision 2-D array data
!===============================================================================

  subroutine read_integer_2D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(2)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_integer_2D_explicit(dset_id, dims, buffer, indep)
    else
      call read_integer_2D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_integer_2D

  subroutine read_integer_2D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(2)
    integer, target,   intent(inout) :: buffer(dims(1),dims(2))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if
  end subroutine read_integer_2D_explicit

!===============================================================================
! WRITE_INTEGER_3D writes integer precision 3-D array data
!===============================================================================

  subroutine write_integer_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_3D

  subroutine write_integer_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(3, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_3D_explicit

!===============================================================================
! READ_INTEGER_3D reads integer precision 3-D array data
!===============================================================================

  subroutine read_integer_3D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(3)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_integer_3D_explicit(dset_id, dims, buffer, indep)
    else
      call read_integer_3D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_integer_3D

  subroutine read_integer_3D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(3)
    integer, target,   intent(inout) :: buffer(dims(1),dims(2),dims(3))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if
  end subroutine read_integer_3D_explicit

!===============================================================================
! WRITE_INTEGER_4D writes integer precision 4-D array data
!===============================================================================

  subroutine write_integer_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_4D

  subroutine write_integer_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(4, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_4D_explicit

!===============================================================================
! READ_INTEGER_4D reads integer precision 4-D array data
!===============================================================================

  subroutine read_integer_4D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(4)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_integer_4D_explicit(dset_id, dims, buffer, indep)
    else
      call read_integer_4D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_integer_4D

  subroutine read_integer_4D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),    intent(in)    :: dset_id
    integer(HSIZE_T),  intent(in)    :: dims(4)
    integer, target,   intent(inout) :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical, optional, intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if
  end subroutine read_integer_4D_explicit

!===============================================================================
! WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine write_long(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer(8),   intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), hdf5_integer8_t, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_long

!===============================================================================
! READ_LONG reads long integer scalar data
!===============================================================================

  subroutine read_long(buffer, obj_id, name, indep)
    integer(8),   target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer        :: hdf5_err
    integer        :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset_id
    type(c_ptr)    :: f_ptr

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr = c_loc(buffer)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, hdf5_integer8_t, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, hdf5_integer8_t, f_ptr, hdf5_err)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_long

!===============================================================================
! WRITE_STRING writes string data
!===============================================================================

  subroutine write_string(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), target   :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    integer(HID_T) :: filetype
    integer(SIZE_T) :: i, n
    character(kind=C_CHAR), allocatable, target :: temp_buffer(:)
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create datatype for HDF5 file based on C char
    n = len_trim(buffer)
    if (n > 0) then
      call h5tcopy_f(H5T_C_S1, filetype, hdf5_err)
      call h5tset_size_f(filetype, n, hdf5_err)

      ! Create dataspace/dataset
      call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
      call h5dcreate_f(group_id, trim(name), filetype, dspace, dset, hdf5_err)

      ! Copy string to temporary buffer
      allocate(temp_buffer(n))
      do i = 1, n
        temp_buffer(i) = buffer(i:i)
      end do

      ! Get pointer to start of string
      f_ptr = c_loc(temp_buffer(1))

      if (using_mpio_device(group_id)) then
#ifdef PHDF5
        call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
        call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
        call h5dwrite_f(dset, filetype, f_ptr, hdf5_err, xfer_prp=plist)
        call h5pclose_f(plist, hdf5_err)
#endif
      else
        call h5dwrite_f(dset, filetype, f_ptr, hdf5_err)
      end if

      call h5dclose_f(dset, hdf5_err)
      call h5sclose_f(dspace, hdf5_err)
      call h5tclose_f(filetype, hdf5_err)
    end if
  end subroutine write_string

!===============================================================================
! READ_STRING reads string data
!===============================================================================

  subroutine read_string(buffer, obj_id, name, indep)
    character(*), target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset_id
    integer(HID_T) :: filetype
    integer(HID_T) :: memtype
    integer(SIZE_T) :: i, n
    character(kind=C_CHAR), allocatable, target :: temp_buffer(:)
    type(c_ptr) :: f_ptr

    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Make sure buffer is large enough
    call h5dget_type_f(dset_id, filetype, hdf5_err)
    call h5tget_size_f(filetype, n, hdf5_err)
    if (n > len(buffer)) then
      call fatal_error("Character buffer is not long enough to &
           &read HDF5 string.")
    end if

    ! Get datatype in memory based on Fortran character
    call h5tcopy_f(H5T_C_S1, memtype, hdf5_err)
    call h5tset_size_f(memtype, n + 1, hdf5_err)

    ! Get pointer to start of string
    allocate(temp_buffer(n))
    f_ptr = c_loc(temp_buffer(1))

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, memtype, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, memtype, f_ptr, hdf5_err)
    end if

    buffer = ''
    do i = 1, n
      if (temp_buffer(i) == C_NULL_CHAR) cycle
      buffer(i:i) = temp_buffer(i)
    end do

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)

    call h5tclose_f(filetype, hdf5_err)
    call h5tclose_f(memtype, hdf5_err)
  end subroutine read_string

!===============================================================================
! WRITE_STRING_1D writes string 1-D array data
!===============================================================================

  subroutine write_string_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), target   :: buffer(:)  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_string_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_string_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_string_1D

  subroutine write_string_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name
    character(*), intent(in), target   :: buffer(dims(1))
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    integer(HID_T) :: filetype
    integer(HID_T) :: memtype
    integer(SIZE_T) :: n
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create datatype for HDF5 file based on C char
    n = maxval(len_trim(buffer))
    call h5tcopy_f(H5T_C_S1, filetype, hdf5_err)
    call h5tset_size_f(filetype, n + 1, hdf5_err)

    ! Create datatype in memory based on Fortran character
    call h5tcopy_f(H5T_FORTRAN_S1, memtype, hdf5_err)
    call h5tset_size_f(memtype, int(len(buffer(1)), SIZE_T), hdf5_err)

    ! Create dataspace/dataset
    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), filetype, dspace, dset, hdf5_err)

    ! Get pointer to start of string
    f_ptr = c_loc(buffer(1)(1:1))

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      if (n > 0) call h5dwrite_f(dset, memtype, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      if (n > 0) call h5dwrite_f(dset, memtype, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5tclose_f(memtype, hdf5_err)
    call h5tclose_f(filetype, hdf5_err)
  end subroutine write_string_1D_explicit

!===============================================================================
! READ_STRING_1D reads string 1-D array data
!===============================================================================

  subroutine read_string_1D(buffer, obj_id, name, indep)
    character(*), target,   intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_string_1D_explicit(dset_id, dims, buffer, indep)
    else
      call read_string_1D_explicit(dset_id, dims, buffer)
    end if
  end subroutine read_string_1D

  subroutine read_string_1D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),       intent(in)    :: dset_id
    integer(HSIZE_T),     intent(in)    :: dims(1)
    character(*), target, intent(inout) :: buffer(dims(1))
    logical, optional,    intent(in)    :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    integer(HID_T) :: space_id
    integer(HID_T) :: filetype
    integer(HID_T) :: memtype
    integer(SIZE_T) :: size
    integer(SIZE_T) :: n
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Get dataset and dataspace
    call h5dget_space_f(dset_id, space_id, hdf5_err)

    ! Make sure buffer is large enough
    call h5dget_type_f(dset_id, filetype, hdf5_err)
    call h5tget_size_f(filetype, size, hdf5_err)
    if (size > len(buffer(1)) + 1) then
      call fatal_error("Character buffer is not long enough to &
           &read HDF5 string array.")
    end if

    ! Get datatype in memory based on Fortran character
    n = len(buffer(1))
    call h5tcopy_f(H5T_FORTRAN_S1, memtype, hdf5_err)
    call h5tset_size_f(memtype, n, hdf5_err)

    ! Get pointer to start of string
    f_ptr = c_loc(buffer(1)(1:1))

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, memtype, f_ptr, hdf5_err, mem_space_id=space_id, &
           xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, memtype, f_ptr, hdf5_err, mem_space_id=space_id)
    end if

    call h5sclose_f(space_id, hdf5_err)
    call h5tclose_f(filetype, hdf5_err)
    call h5tclose_f(memtype, hdf5_err)
  end subroutine read_string_1D_explicit

!===============================================================================
! WRITE_ATTRIBUTE_STRING
!===============================================================================

  subroutine write_attribute_string(obj_id, name, buffer)
    integer(HID_T), intent(in)       :: obj_id    ! object to write attribute to
    character(*), intent(in)         :: name      ! name of attribute
    character(*), intent(in), target :: buffer    ! string to write

    integer        :: hdf5_err
    integer(HID_T) :: dspace_id
    integer(HID_T) :: attr_id
    integer(HID_T) :: filetype
    integer(SIZE_T) :: i
    integer(SIZE_T) :: n
    character(kind=C_CHAR), allocatable, target :: temp_buffer(:)
    type(c_ptr) :: f_ptr

    ! Create datatype for HDF5 file based on C char
    n = len_trim(buffer)
    if (n > 0) then
      call h5tcopy_f(H5T_C_S1, filetype, hdf5_err)
      call h5tset_size_f(filetype, n, hdf5_err)

      ! Create memory space and attribute
      call h5screate_f(H5S_SCALAR_F, dspace_id, hdf5_err)
      call h5acreate_f(obj_id, trim(name), filetype, dspace_id, &
           attr_id, hdf5_err)

      ! Copy string to temporary buffer
      allocate(temp_buffer(n))
      do i = 1, n
        temp_buffer(i) = buffer(i:i)
      end do

      ! Write attribute
      f_ptr = c_loc(buffer(1:1))
      call h5awrite_f(attr_id, filetype, f_ptr, hdf5_err)

      ! Close attribute
      call h5aclose_f(attr_id, hdf5_err)
      call h5sclose_f(dspace_id, hdf5_err)
      call h5tclose_f(filetype, hdf5_err)
    end if
  end subroutine write_attribute_string

  subroutine read_attribute_double(buffer, obj_id, name)
    real(8),        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    integer        :: hdf5_err
    integer(HID_T) :: attr_id
    type(c_ptr)    :: f_ptr

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_double

  subroutine write_attribute_double(obj_id, name, buffer)
    integer(HID_T), intent(in)  :: obj_id
    character(*),   intent(in)  :: name
    real(8), intent(in), target :: buffer

    integer        :: hdf5_err
    integer(HID_T) :: dspace_id
    integer(HID_T) :: attr_id
    type(C_PTR)    :: f_ptr

    call h5screate_f(H5S_SCALAR_F, dspace_id, hdf5_err)
    call h5acreate_f(obj_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, &
         attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
    call h5sclose_f(dspace_id, hdf5_err)
  end subroutine write_attribute_double

  subroutine read_attribute_double_1D(buffer, obj_id, name)
    real(8), target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer          :: hdf5_err
    integer(HID_T)   :: space_id
    integer(HID_T)   :: attr_id
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: maxdims(1)

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    if (allocated(buffer)) then
      dims(:) = shape(buffer)
    else
      call h5aget_space_f(attr_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      allocate(buffer(dims(1)))
      call h5sclose_f(space_id, hdf5_err)
    end if

    call read_attribute_double_1D_explicit(attr_id, dims, buffer)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_double_1D

  subroutine read_attribute_double_1D_explicit(attr_id, dims, buffer)
    integer(HID_T),   intent(in)    :: attr_id
    integer(HSIZE_T), intent(in)    :: dims(1)
    real(8), target,  intent(inout) :: buffer(dims(1))

    integer        :: hdf5_err
    type(c_ptr)    :: f_ptr

    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
  end subroutine read_attribute_double_1D_explicit

  subroutine write_attribute_double_1D(obj_id, name, buffer)
    integer(HID_T),  intent(in) :: obj_id
    character(*),    intent(in) :: name
    real(8), target, intent(in) :: buffer(:)

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    call write_attribute_double_1D_explicit(obj_id, dims, name, buffer)
  end subroutine write_attribute_double_1D

  subroutine write_attribute_double_1D_explicit(obj_id, dims, name, buffer)
    integer(HID_T),   intent(in) :: obj_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*),     intent(in) :: name
    real(8), target,  intent(in) :: buffer(dims(1))

    integer        :: hdf5_err
    integer(HID_T) :: dspace_id
    integer(HID_T) :: attr_id
    type(C_PTR)    :: f_ptr

    call h5screate_simple_f(1, dims, dspace_id, hdf5_err)
    call h5acreate_f(obj_id, trim(name), H5T_NATIVE_DOUBLE, dspace_id, &
         attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
    call h5sclose_f(dspace_id, hdf5_err)
  end subroutine write_attribute_double_1D_explicit

  subroutine read_attribute_double_2D(buffer, obj_id, name)
    real(8), target, allocatable, intent(inout) :: buffer(:,:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer          :: hdf5_err
    integer(HID_T)   :: space_id
    integer(HID_T)   :: attr_id
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T) :: maxdims(2)

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    if (allocated(buffer)) then
      dims(:) = shape(buffer)
    else
      call h5aget_space_f(attr_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      allocate(buffer(dims(1), dims(2)))
      call h5sclose_f(space_id, hdf5_err)
    end if

    call read_attribute_double_2D_explicit(attr_id, dims, buffer)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_double_2D

  subroutine read_attribute_double_2D_explicit(attr_id, dims, buffer)
    integer(HID_T),   intent(in)    :: attr_id
    integer(HSIZE_T), intent(in)    :: dims(2)
    real(8), target,  intent(inout) :: buffer(dims(1),dims(2))

    integer        :: hdf5_err
    type(c_ptr)    :: f_ptr

    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
  end subroutine read_attribute_double_2D_explicit

  subroutine read_attribute_integer(buffer, obj_id, name)
    integer,        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    integer        :: hdf5_err
    integer(HID_T) :: attr_id
    type(c_ptr)    :: f_ptr

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_integer

  subroutine write_attribute_integer(obj_id, name, buffer)
    integer(HID_T), intent(in)  :: obj_id
    character(*),   intent(in)  :: name
    integer, intent(in), target :: buffer

    integer        :: hdf5_err
    integer(HID_T) :: dspace_id
    integer(HID_T) :: attr_id
    type(C_PTR)    :: f_ptr

    call h5screate_f(H5S_SCALAR_F, dspace_id, hdf5_err)
    call h5acreate_f(obj_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, &
         attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
    call h5sclose_f(dspace_id, hdf5_err)
  end subroutine write_attribute_integer

  subroutine read_attribute_integer_1D(buffer, obj_id, name)
    integer, target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer          :: hdf5_err
    integer(HID_T)   :: space_id
    integer(HID_T)   :: attr_id
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: maxdims(1)

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    if (allocated(buffer)) then
      dims(:) = shape(buffer)
    else
      call h5aget_space_f(attr_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      allocate(buffer(dims(1)))
      call h5sclose_f(space_id, hdf5_err)
    end if

    call read_attribute_integer_1D_explicit(attr_id, dims, buffer)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_integer_1D

  subroutine read_attribute_integer_1D_explicit(attr_id, dims, buffer)
    integer(HID_T),   intent(in)    :: attr_id
    integer(HSIZE_T), intent(in)    :: dims(1)
    integer, target,  intent(inout) :: buffer(dims(1))

    integer        :: hdf5_err
    type(c_ptr)    :: f_ptr

    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
  end subroutine read_attribute_integer_1D_explicit

  subroutine write_attribute_integer_1D(obj_id, name, buffer)
    integer(HID_T),  intent(in) :: obj_id
    character(*),    intent(in) :: name
    integer, target, intent(in) :: buffer(:)

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    call write_attribute_integer_1D_explicit(obj_id, dims, name, buffer)
  end subroutine write_attribute_integer_1D

  subroutine write_attribute_integer_1D_explicit(obj_id, dims, name, buffer)
    integer(HID_T),   intent(in) :: obj_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*),     intent(in) :: name
    integer, target,  intent(in) :: buffer(dims(1))

    integer        :: hdf5_err
    integer(HID_T) :: dspace_id
    integer(HID_T) :: attr_id
    type(C_PTR)    :: f_ptr

    call h5screate_simple_f(1, dims, dspace_id, hdf5_err)
    call h5acreate_f(obj_id, trim(name), H5T_NATIVE_INTEGER, dspace_id, &
         attr_id, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)
    call h5sclose_f(dspace_id, hdf5_err)
  end subroutine write_attribute_integer_1D_explicit

  subroutine read_attribute_integer_2D(buffer, obj_id, name)
    integer, target, allocatable, intent(inout) :: buffer(:,:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer          :: hdf5_err
    integer(HID_T)   :: space_id
    integer(HID_T)   :: attr_id
    integer(HSIZE_T) :: dims(2)
    integer(HSIZE_T) :: maxdims(2)

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    if (allocated(buffer)) then
      dims(:) = shape(buffer)
    else
      call h5aget_space_f(attr_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      allocate(buffer(dims(1), dims(2)))
      call h5sclose_f(space_id, hdf5_err)
    end if

    call read_attribute_integer_2D_explicit(attr_id, dims, buffer)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_integer_2D

  subroutine read_attribute_integer_2D_explicit(attr_id, dims, buffer)
    integer(HID_T),   intent(in)    :: attr_id
    integer(HSIZE_T), intent(in)    :: dims(2)
    integer, target,  intent(inout) :: buffer(dims(1),dims(2))

    integer        :: hdf5_err
    type(c_ptr)    :: f_ptr

    f_ptr = c_loc(buffer)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
  end subroutine read_attribute_integer_2D_explicit

  subroutine read_attribute_string(buffer, obj_id, name)
    character(*),   intent(inout) :: buffer  ! read data to here
    integer(HID_T), intent(in)    :: obj_id
    character(*),   intent(in)    :: name    ! name for data

    integer         :: hdf5_err
    integer(HID_T)  :: attr_id ! data set handle
    integer(HID_T)  :: filetype
    integer(HID_T)  :: memtype
    integer(SIZE_T) :: i
    integer(SIZE_T) :: size
    character(kind=C_CHAR), allocatable, target :: temp_buffer(:)
    type(c_ptr)     :: f_ptr

    ! Get dataset and dataspace
    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    ! Make sure buffer is large enough
    call h5aget_type_f(attr_id, filetype, hdf5_err)
    call h5tget_size_f(filetype, size, hdf5_err)
    allocate(temp_buffer(size))
    if (size > len(buffer)) then
      call fatal_error("Character buffer is not long enough to &
           &read HDF5 string.")
    end if

    ! Get datatype in memory based on Fortran character
    call h5tcopy_f(H5T_C_S1, memtype, hdf5_err)
    call h5tset_size_f(memtype, size + 1, hdf5_err)

    ! Get pointer to start of string
    f_ptr = c_loc(temp_buffer(1))

    call h5aread_f(attr_id, memtype, f_ptr, hdf5_err)
    buffer = ''
    do i = 1, size
      if (temp_buffer(i) == C_NULL_CHAR) cycle
      buffer(i:i) = temp_buffer(i)
    end do

    call h5aclose_f(attr_id, hdf5_err)
    call h5tclose_f(filetype, hdf5_err)
    call h5tclose_f(memtype, hdf5_err)
  end subroutine read_attribute_string

  subroutine read_attribute_string_1D(buffer, obj_id, name)
    character(*), target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T), intent(in) :: obj_id
    character(*),   intent(in) :: name

    integer          :: hdf5_err
    integer(HID_T)   :: space_id
    integer(HID_T)   :: attr_id
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: maxdims(1)

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)

    if (allocated(buffer)) then
      dims(:) = shape(buffer)
    else
      call h5aget_space_f(attr_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      allocate(buffer(dims(1)))
      call h5sclose_f(space_id, hdf5_err)
    end if

    call read_attribute_string_1D_explicit(attr_id, dims, buffer)
    call h5aclose_f(attr_id, hdf5_err)
  end subroutine read_attribute_string_1D

  subroutine read_attribute_string_1D_explicit(attr_id, dims, buffer)
    integer(HID_T),       intent(in)    :: attr_id
    integer(HSIZE_T),     intent(in)    :: dims(1)
    character(*), target, intent(inout) :: buffer(dims(1))

    integer :: hdf5_err
    integer(HID_T) :: filetype
    integer(HID_T) :: memtype
    integer(SIZE_T) :: size
    integer(SIZE_T) :: n
    type(c_ptr) :: f_ptr

    ! Make sure buffer is large enough
    call h5aget_type_f(attr_id, filetype, hdf5_err)
    call h5tget_size_f(filetype, size, hdf5_err)
    if (size > len(buffer(1)) + 1) then
      call fatal_error("Character buffer is not long enough to &
           &read HDF5 string array.")
    end if

    ! Get datatype in memory based on Fortran character
    n = len(buffer(1))
    call h5tcopy_f(H5T_FORTRAN_S1, memtype, hdf5_err)
    call h5tset_size_f(memtype, n, hdf5_err)

    ! Get pointer to start of string
    f_ptr = c_loc(buffer(1)(1:1))

    call h5aread_f(attr_id, memtype, f_ptr, hdf5_err)

    call h5tclose_f(filetype, hdf5_err)
    call h5tclose_f(memtype, hdf5_err)
  end subroutine read_attribute_string_1D_explicit

  subroutine read_attribute_logical(buffer, obj_id, name)
    logical,        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    integer, target :: int_buffer
    integer         :: hdf5_err
    integer(HID_T)  :: attr_id
    type(c_ptr)     :: f_ptr

    call h5aopen_f(obj_id, trim(name), attr_id, hdf5_err)
    f_ptr = c_loc(int_buffer)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    call h5aclose_f(attr_id, hdf5_err)

    ! Convert to Fortran logical
    if (int_buffer == 0) then
      buffer = .false.
    else
      buffer = .true.
    end if
  end subroutine read_attribute_logical

  subroutine get_shape(obj_id, dims)
    integer(HID_T),   intent(in)  :: obj_id
    integer(HSIZE_T), intent(out) :: dims(:)

    integer          :: hdf5_err
    integer          :: type
    integer(HID_T)   :: space_id
    integer(HSIZE_T) :: maxdims(size(dims))

    call h5iget_type_f(obj_id, type, hdf5_err)
    if (type == H5I_DATASET_F) then
      call h5dget_space_f(obj_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      call h5sclose_f(space_id, hdf5_err)
    elseif (type == H5I_ATTR_F) then
      call h5aget_space_f(obj_id, space_id, hdf5_err)
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdf5_err)
      call h5sclose_f(space_id, hdf5_err)
    end if
  end subroutine get_shape

  subroutine get_ndims(obj_id, ndims)
    integer(HID_T), intent(in)  :: obj_id
    integer,        intent(out) :: ndims

    integer          :: hdf5_err
    integer          :: type
    integer(HID_T)   :: space_id

    call h5iget_type_f(obj_id, type, hdf5_err)
    if (type == H5I_DATASET_F) then
      call h5dget_space_f(obj_id, space_id, hdf5_err)
      call h5sget_simple_extent_ndims_f(space_id, ndims, hdf5_err)
      call h5sclose_f(space_id, hdf5_err)
    elseif (type == H5I_ATTR_F) then
      call h5aget_space_f(obj_id, space_id, hdf5_err)
      call h5sget_simple_extent_ndims_f(space_id, ndims, hdf5_err)
      call h5sclose_f(space_id, hdf5_err)
    end if
  end subroutine get_ndims

  function using_mpio_device(obj_id) result(mpio)
    integer(HID_T), intent(in) :: obj_id
    logical :: mpio

    integer :: hdf5_err
    integer(HID_T) :: driver
    integer(HID_T) :: file_id
    integer(HID_T) :: fapl_id

    ! Determine file that this object is part of
    call h5iget_file_id_f(obj_id, file_id, hdf5_err)

    ! Get file access property list
    call h5fget_access_plist_f(file_id, fapl_id, hdf5_err)

    ! Get low-level driver identifier
    call h5pget_driver_f(fapl_id, driver, hdf5_err)

    ! Close file access property list access
    call h5pclose_f(fapl_id, hdf5_err)

    ! Close file access -- note that this only decreases the reference count so
    ! that the file is not actually closed
    call h5fclose_f(file_id, hdf5_err)

    mpio = (driver == H5FD_MPIO_F)
  end function using_mpio_device

!===============================================================================
! READ_COMPLEX_2D reads double precision complex 2-D array data as output by
! the h5py HDF5 python module.
!===============================================================================

  subroutine read_complex_2D(buffer, obj_id, name, indep)
    complex(8),   target,   intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer          :: hdf5_err
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(2)

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      call h5dopen_f(obj_id, trim(name), dset_id, hdf5_err)
    else
      dset_id = obj_id
    end if

    dims(:) = shape(buffer)

    if (present(indep)) then
      call read_complex_2D_explicit(dset_id, dims, buffer, indep)
    else
      call read_complex_2D_explicit(dset_id, dims, buffer)
    end if

    if (present(name)) call h5dclose_f(dset_id, hdf5_err)
  end subroutine read_complex_2D

  subroutine read_complex_2D_explicit(dset_id, dims, buffer, indep)
    integer(HID_T),     intent(in)    :: dset_id
    integer(HSIZE_T),   intent(in)    :: dims(2)
    complex(8), target, intent(inout) :: buffer(dims(1), dims(2))
    logical, optional, intent(in)     :: indep   ! independent I/O

    real(8), target :: buffer_r(dims(1), dims(2))
    real(8), target :: buffer_i(dims(1), dims(2))

    integer(HSIZE_T) :: i, j

    integer :: hdf5_err
    integer :: data_xfer_mode
#ifdef PHDF5
    integer(HID_T) :: plist   ! property list
#endif
    type(c_ptr) :: f_ptr_r, f_ptr_i

    ! Components needed for complex type support
    integer(HID_T)   :: dtype_real
    integer(HID_T)   :: dtype_imag
    integer(SIZE_T)  :: size_double

    ! Create the complex type
    call h5tget_size_f(H5T_NATIVE_DOUBLE, size_double, hdf5_err)

    ! Insert the 'r' and 'i' identifiers
    call h5tcreate_f(H5T_COMPOUND_F, size_double, dtype_real, hdf5_err)
    call h5tcreate_f(H5T_COMPOUND_F, size_double, dtype_imag, hdf5_err)
    call h5tinsert_f(dtype_real, "r", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdf5_err)
    call h5tinsert_f(dtype_imag, "i", 0_SIZE_T, H5T_NATIVE_DOUBLE, hdf5_err)

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    f_ptr_r = c_loc(buffer_r)
    f_ptr_i = c_loc(buffer_i)

    if (using_mpio_device(dset_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset_id, dtype_real, f_ptr_r, hdf5_err, xfer_prp=plist)
      call h5dread_f(dset_id, dtype_imag, f_ptr_i, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset_id, dtype_real, f_ptr_r, hdf5_err)
      call h5dread_f(dset_id, dtype_imag, f_ptr_i, hdf5_err)
    end if

    ! Reconstitute the complex numbers
    do i = 1, dims(1)
      do j = 1, dims(2)
        buffer(i, j) = cmplx(buffer_r(i,j), buffer_i(i,j), kind=8)
      end do
    end do
  end subroutine read_complex_2D_explicit

end module hdf5_interface
