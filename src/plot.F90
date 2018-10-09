module plot

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error,           only: fatal_error, write_message
  use geometry,        only: find_cell, check_cell_overlap
  use geometry_header, only: Cell, root_universe, cells
  use hdf5_interface
  use output,          only: time_stamp
  use material_header, only: materials
  use mesh_header,     only: meshes, RegularMesh
  use particle_header
  use plot_header
  use progress_header, only: ProgressBar
  use settings,        only: check_overlaps
  use string,          only: to_str

  implicit none
  private

  public :: openmc_plot_geometry

  integer, parameter :: RED = 1
  integer, parameter :: GREEN = 2
  integer, parameter :: BLUE = 3

  interface
    subroutine position_rgb(p, pl, rgb, id) bind(C)
      import Particle, ObjectPlot, C_INT
      type(Particle), intent(inout) :: p
      type(ObjectPlot), intent(in) :: pl
      integer(C_INT),   intent(out) :: rgb(3)
      integer(C_INT), intent(out) :: id
    end subroutine position_rgb


    subroutine create_ppm(pl) bind(C)
      import ObjectPlot
      type(ObjectPlot), intent(in) :: pl
    end subroutine create_ppm
    
  end interface
contains

!===============================================================================
! RUN_PLOT controls the logic for making one or many plots
!===============================================================================

  function openmc_plot_geometry() result(err) bind(C)
    integer(C_INT) :: err

    integer :: i ! loop index for plots

    do i = 1, n_plots
      associate (pl => plots(i))
        ! Display output message
        call write_message("Processing plot " // trim(to_str(pl % id)) &
             // ": " // trim(pl % path_plot) // " ...", 5)

        if (pl % type == PLOT_TYPE_SLICE) then
          ! create 2d image
          call create_ppm(pl)
        else if (pl % type == PLOT_TYPE_VOXEL) then
          ! create dump for 3D silomesh utility script
          call create_voxel(pl)
        end if
      end associate
    end do

    err = 0
  end function openmc_plot_geometry

!===============================================================================
! CREATE_VOXEL outputs a binary file that can be input into silomesh for 3D
! geometry visualization.  It works the same way as create_ppm by dragging a
! particle across the geometry for the specified number of voxels. The first
! 3 int(4)'s in the binary are the number of x, y, and z voxels.  The next 3
! real(8)'s are the widths of the voxels in the x, y, and z directions. The next
! 3 real(8)'s are the x, y, and z coordinates of the lower left point. Finally
! the binary is filled with entries of four int(4)'s each. Each 'row' in the
! binary contains four int(4)'s: 3 for x,y,z position and 1 for cell or material
! id.  For 1 million voxels this produces a file of approximately 15MB.
!===============================================================================

  subroutine create_voxel(pl)
    type(ObjectPlot), intent(in) :: pl

    integer(C_INT) :: x, y, z      ! voxel location indices
    integer :: rgb(3)       ! colors (red, green, blue) from 0-255
    integer :: id           ! id of cell or material
    integer(C_INT), target :: data(pl%pixels(3),pl%pixels(2))
    integer(HID_T) :: file_id
    integer(HID_T) :: dspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dims(3)
    real(8) :: vox(3)       ! x, y, and z voxel widths
    real(8) :: ll(3)        ! lower left starting point for each sweep direction
    type(Particle)    :: p
    type(ProgressBar) :: progress

    interface
      subroutine voxel_init(file_id, dims, dspace, dset, memspace) bind(C)
        import HID_T, HSIZE_T
        integer(HID_T), value :: file_id
        integer(HSIZE_T), intent(in) :: dims(*)
        integer(HID_T), intent(out) :: dspace
        integer(HID_T), intent(out) :: dset
        integer(HID_T), intent(out) :: memspace
      end subroutine voxel_init

      subroutine voxel_write_slice(x, dspace, dset, memspace, buf) bind(C)
        import C_INT, HID_T, C_PTR
        integer(C_INT), value :: x
        integer(HID_T), value :: dspace
        integer(HID_T), value :: dset
        integer(HID_T), value :: memspace
        type(C_PTR), value :: buf
      end subroutine voxel_write_slice

      subroutine voxel_finalize(dspace, dset, memspace) bind(C)
        import HID_T
        integer(HID_T), value :: dspace
        integer(HID_T), value :: dset
        integer(HID_T), value :: memspace
      end subroutine voxel_finalize
    end interface

    ! compute voxel widths in each direction
    vox = pl % width/dble(pl % pixels)

    ! initial particle position
    ll = pl % origin - pl % width / TWO

    ! allocate and initialize particle
    call particle_initialize(p)
    p % coord(1) % xyz = ll
    p % coord(1) % uvw = [ HALF, HALF, HALF ]
    p % coord(1) % universe = root_universe

    ! Open binary plot file for writing
    file_id = file_open(pl%path_plot, 'w')

    ! write header info
    call write_attribute(file_id, "filetype", 'voxel')
    call write_attribute(file_id, "version", VERSION_VOXEL)
    call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
    call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

    ! Write current date and time
    call write_attribute(file_id, "date_and_time", time_stamp())

    call write_attribute(file_id, "num_voxels", pl%pixels)
    call write_attribute(file_id, "voxel_width", vox)
    call write_attribute(file_id, "lower_left", ll)

    ! Create dataset for voxel data -- note that the dimensions are reversed
    ! since we want the order in the file to be z, y, x
    dims(:) = pl % pixels
    call voxel_init(file_id, dims, dspace, dset, memspace)

    ! move to center of voxels
    ll = ll + vox / TWO

    do x = 1, pl % pixels(1)
      call progress % set_value(dble(x)/dble(pl % pixels(1))*100)
      do y = 1, pl % pixels(2)
        do z = 1, pl % pixels(3)
          ! get voxel color
          call position_rgb(p, pl, rgb, id)
          ! write to plot file
          data(z,y) = id

          ! advance particle in z direction
          p % coord(1) % xyz(3) = p % coord(1) % xyz(3) + vox(3)
        end do

        ! advance particle in y direction
        p % coord(1) % xyz(2) = p % coord(1) % xyz(2) + vox(2)
        p % coord(1) % xyz(3) = ll(3)
      end do

      ! advance particle in y direction
      p % coord(1) % xyz(1) = p % coord(1) % xyz(1) + vox(1)
      p % coord(1) % xyz(2) = ll(2)
      p % coord(1) % xyz(3) = ll(3)

      ! Write to HDF5 dataset
      call voxel_write_slice(x, dspace, dset, memspace, c_loc(data))
    end do

    call voxel_finalize(dspace, dset, memspace)
    call file_close(file_id)

  end subroutine create_voxel

end module plot
