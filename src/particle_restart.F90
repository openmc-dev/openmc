module particle_restart 

  use bank_header,     only: Bank
  use constants
  use geometry_header, only: BASE_UNIVERSE
  use global
  use particle_header, only: deallocate_coord

#ifdef HDF5
  use hdf5
  use h5lt
#endif

  implicit none
  private
  public ::  write_particle_restart, run_particle_restart

  integer(HID_T) :: hdf5_particle_file

contains

!===============================================================================
! WRITE_PARTICLE_RESTART
!===============================================================================

  subroutine write_particle_restart()

    character(MAX_FILE_LEN) :: filename
    integer(HSIZE_T)        :: dims1(1)
    type(Bank), pointer     :: src => null()

    ! set up file name
    filename = 'particle_'//trim(int4_to_str(rank))//'.h5'

    ! create hdf5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_particle_file, hdf5_err)

    ! get information about source particle
    src => source_bank(current_work)

    ! write data to file
    call hdf5_write_integer(hdf5_particle_file, 'current_batch', current_batch)
    call hdf5_write_integer(hdf5_particle_file, 'gen_per_batch', gen_per_batch)
    call hdf5_write_integer(hdf5_particle_file, 'current_gen', current_gen)
    call hdf5_write_long(hdf5_particle_file, 'n_particles', n_particles)
    call hdf5_write_long(hdf5_particle_file, 'id', p % id)
    call hdf5_write_double(hdf5_particle_file, 'weight', src % wgt)
    call hdf5_write_double(hdf5_particle_file, 'energy', src % E)
    dims1 = (/3/)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'xyz', 1, dims1, &
         src % xyz, hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'uvw', 1, dims1, &
         src % uvw, hdf5_err)

    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)

  end subroutine write_particle_restart

!===============================================================================
! READ_PARTICLE_RESTART
!===============================================================================

  subroutine read_particle_restart()

    integer(HSIZE_T)        :: dims1(1)

    ! open hdf5 file
    call h5fopen_f(path_particle_restart, H5F_ACC_RDONLY_F, hdf5_particle_file,&
                   hdf5_err)

    ! read data from file
    call hdf5_read_integer(hdf5_particle_file, 'current_batch', current_batch)
    call hdf5_read_integer(hdf5_particle_file, 'gen_per_batch', gen_per_batch)
    call hdf5_read_integer(hdf5_particle_file, 'current_gen', current_gen)
    call hdf5_read_long(hdf5_particle_file, 'n_particles', n_particles)
    call hdf5_read_long(hdf5_particle_file, 'id', p % id)
    call hdf5_read_double(hdf5_particle_file, 'weight', p % wgt)
    call hdf5_read_double(hdf5_particle_file, 'energy', p % E)
    dims1 = (/3/)
    call h5ltread_dataset_double_f(hdf5_particle_file, 'xyz', p % coord % xyz, &
         dims1, hdf5_err)
    call h5ltread_dataset_double_f(hdf5_particle_file, 'uvw', p % coord % uvw, &
         dims1, hdf5_err)

    ! set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord % xyz
    p % last_E   = p % E

    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)

  end subroutine read_particle_restart

!===============================================================================
! RUN_PARTICLE_RESTART
!===============================================================================

  subroutine run_particle_restart()

    ! initialize the particle to be tracked
    call initialize_particle()

    ! read in the restart information
    call read_particle_restart()

    ! set all tallies to 0 for now (just tracking errors)
    n_tallies = 0

    ! compute seed

  end subroutine run_particle_restart

!===============================================================================
! INITIALIZE_PARTICLE
!===============================================================================

  subroutine initialize_particle()

    ! Allocate particle
    allocate(p)

    ! Set particle to neutron that's alive
    p % type  = NEUTRON
    p % alive = .true.

    ! clear attributes
    p % surface       = NONE
    p % cell_born     = NONE
    p % material      = NONE
    p % last_material = NONE
    p % wgt           = ONE
    p % last_wgt      = ONE
    p % absorb_wgt    = ZERO
    p % n_bank        = 0
    p % wgt_bank      = ZERO
    p % n_collision   = 0

    ! remove any original coordinates
    call deallocate_coord(p % coord0)

    ! Set up base level coordinates
    allocate(p % coord0)
    p % coord0 % universe = BASE_UNIVERSE
    p % coord             => p % coord0

  end subroutine initialize_particle

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
! HDF5_WRITE_LONG
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer)

    integer(HID_T),     intent(in) :: group
    character(*),       intent(in) :: name
    integer(8), target, intent(in) :: buffer

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr

    ! Create dataspace and dataset
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)
    call h5dcreate_f(group, name, hdf5_integer8_t, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

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
! HDF5_READ_LONG
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer)

    integer(HID_T),     intent(in)  :: group
    character(*),       intent(in)  :: name
    integer(8), target, intent(out) :: buffer

    integer(HID_T) :: dset
    type(c_ptr)    :: f_ptr

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get pointer to buffer
    f_ptr = c_loc(buffer)

    ! Read data from dataset
    call h5dread_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)

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
! INT4_TO_STR converts an integer(4) to a string.
!===============================================================================

  function int4_to_str(num) result(str)

    integer, intent(in) :: num
    character(11) :: str

    write (str, '(I11)') num
    str = adjustl(str)

  end function int4_to_str

end module particle_restart
