module plot

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants
  use error,           only: fatal_error, write_message
  use geometry,        only: find_cell, check_cell_overlap
  use geometry_header, only: Cell, root_universe, cells
  use hdf5_interface
  use output,          only: time_stamp
  use material_header, only: materials
  use particle_header, only: LocalCoord, Particle
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

contains

!===============================================================================
! RUN_PLOT controls the logic for making one or many plots
!===============================================================================

  subroutine openmc_plot_geometry() bind(C)

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

  end subroutine openmc_plot_geometry

!===============================================================================
! POSITION_RGB computes the red/green/blue values for a given plot with the
! current particle's position
!===============================================================================

  subroutine position_rgb(p, pl, rgb, id)
    type(Particle), intent(inout) :: p
    type(ObjectPlot), intent(in)  :: pl
    integer, intent(out)          :: rgb(3)
    integer, intent(out)          :: id

    integer :: j
    logical :: found_cell

    p % n_coord = 1

    call find_cell(p, found_cell)
    j = p % n_coord
    if (check_overlaps) call check_cell_overlap(p)

    ! Set coordinate level if specified
    if (pl % level >= 0) j = pl % level + 1

    if (.not. found_cell) then
      ! If no cell, revert to default color
      rgb = pl % not_found % rgb
      id = -1
    else
      if (pl % color_by == PLOT_COLOR_MATS) then
        ! Assign color based on material
        associate (c => cells(p % coord(j) % cell))
          if (c % type == FILL_UNIVERSE) then
            ! If we stopped on a middle universe level, treat as if not found
            rgb = pl % not_found % rgb
            id = -1
          else if (p % material == MATERIAL_VOID) then
            ! By default, color void cells white
            rgb = 255
            id = -1
          else
            rgb = pl % colors(p % material) % rgb
            id = materials(p % material) % id
          end if
        end associate
      else if (pl % color_by == PLOT_COLOR_CELLS) then
        ! Assign color based on cell
        rgb = pl % colors(p % coord(j) % cell) % rgb
        id = cells(p % coord(j) % cell) % id
      else
        rgb = 0
        id = -1
      end if
    end if

  end subroutine position_rgb

!===============================================================================
! CREATE_PPM creates an image based on user input from a plots.xml <plot>
! specification in the portable pixmap format (PPM)
!===============================================================================

  subroutine create_ppm(pl)
    type(ObjectPlot), intent(in) :: pl

    integer :: in_i
    integer :: out_i
    integer :: x, y      ! pixel location
    integer :: rgb(3)    ! colors (red, green, blue) from 0-255
    integer :: id
    integer :: height, width
    real(8) :: in_pixel
    real(8) :: out_pixel
    real(8) :: xyz(3)
    integer, allocatable :: data(:,:,:)
    type(Particle)    :: p

    width = pl % pixels(1)
    height = pl % pixels(2)

    in_pixel  = pl % width(1)/dble(width)
    out_pixel = pl % width(2)/dble(height)

    ! Allocate and initialize results array
    allocate(data(3, width, height))
    data(:,:,:) = 0

    if (pl % basis == PLOT_BASIS_XY) then
      in_i  = 1
      out_i = 2
      xyz(1) = pl % origin(1) - pl % width(1) / TWO
      xyz(2) = pl % origin(2) + pl % width(2) / TWO
      xyz(3) = pl % origin(3)
    else if (pl % basis == PLOT_BASIS_XZ) then
      in_i  = 1
      out_i = 3
      xyz(1) = pl % origin(1) - pl % width(1) / TWO
      xyz(2) = pl % origin(2)
      xyz(3) = pl % origin(3) + pl % width(2) / TWO
    else if (pl % basis == PLOT_BASIS_YZ) then
      in_i  = 2
      out_i = 3
      xyz(1) = pl % origin(1)
      xyz(2) = pl % origin(2) - pl % width(1) / TWO
      xyz(3) = pl % origin(3) + pl % width(2) / TWO
    end if

    ! allocate and initialize particle
    call p % initialize()
    p % coord(1) % xyz = xyz
    p % coord(1) % uvw = [ HALF, HALF, HALF ]
    p % coord(1) % universe = root_universe

!$omp parallel do firstprivate(p) private(x, rgb, id) reduction(+ : data)
    do y = 1, height
      ! Set y coordinate
      p % coord(1) % xyz(out_i) = xyz(out_i) - out_pixel*(y - 1)
      do x = 1, width
        ! Set x coordinate
        p % coord(1) % xyz(in_i) = xyz(in_i) + in_pixel*(x - 1)

        ! get pixel color
        call position_rgb(p, pl, rgb, id)

        ! Create a pixel at (x,y) with color (r,g,b)
        data(:, x, y) = rgb
      end do
    end do
!$omp end parallel do

    ! Draw tally mesh boundaries on the image if requested
    if (associated(pl % meshlines_mesh)) call draw_mesh_lines(pl, data)

    ! Write out the ppm to a file
    call output_ppm(pl, data)

  end subroutine create_ppm

!===============================================================================
! DRAW_MESH_LINES draws mesh line boundaries on an image
!===============================================================================

  subroutine draw_mesh_lines(pl, data)
    type(ObjectPlot), intent(in)    :: pl
    integer,          intent(inout) :: data(:,:,:)

    logical :: in_mesh
    integer :: out_, in_  ! pixel location
    integer :: rgb(3)     ! RGB color for meshlines pixels
    integer :: outrange(2), inrange(2) ! range of pixel locations
    integer :: i, j       ! loop indices
    integer :: plus
    integer :: ijk_ll(3)  ! mesh bin ijk indicies of plot lower left
    integer :: ijk_ur(3)  ! mesh bin ijk indicies of plot upper right
    integer :: outer, inner
    real(8) :: frac
    real(8) :: width(3)   ! real widths of the plot
    real(8) :: xyz_ll_plot(3)  ! lower left xyz of plot image
    real(8) :: xyz_ur_plot(3)  ! upper right xyz of plot image
    real(8) :: xyz_ll(3)  ! lower left xyz
    real(8) :: xyz_ur(3)  ! upper right xyz

    rgb(:) = pl % meshlines_color % rgb

    select case (pl % basis)
      case(PLOT_BASIS_XY)
        outer = 1
        inner = 2
      case(PLOT_BASIS_XZ)
        outer = 1
        inner = 3
      case(PLOT_BASIS_YZ)
        outer = 2
        inner = 3
    end select

    xyz_ll_plot = pl % origin
    xyz_ur_plot = pl % origin

    xyz_ll_plot(outer) = pl % origin(outer) - pl % width(1) / TWO
    xyz_ll_plot(inner) = pl % origin(inner) - pl % width(2) / TWO
    xyz_ur_plot(outer) = pl % origin(outer) + pl % width(1) / TWO
    xyz_ur_plot(inner) = pl % origin(inner) + pl % width(2) / TWO

    width = xyz_ur_plot - xyz_ll_plot

    associate (m => pl % meshlines_mesh)
      call m % get_indices(xyz_ll_plot, ijk_ll(:m % n_dimension), in_mesh)
      call m % get_indices(xyz_ur_plot, ijk_ur(:m % n_dimension), in_mesh)

      ! sweep through all meshbins on this plane and draw borders
      do i = ijk_ll(outer), ijk_ur(outer)
        do j = ijk_ll(inner), ijk_ur(inner)
          ! check if we're in the mesh for this ijk
          if (i > 0 .and. i <= m % dimension(outer) .and. &
               j > 0 .and. j <= m % dimension(inner)) then

            ! get xyz's of lower left and upper right of this mesh cell
            xyz_ll(outer) = m % lower_left(outer) + m % width(outer) * (i - 1)
            xyz_ll(inner) = m % lower_left(inner) + m % width(inner) * (j - 1)
            xyz_ur(outer) = m % lower_left(outer) + m % width(outer) * i
            xyz_ur(inner) = m % lower_left(inner) + m % width(inner) * j

            ! map the xyz ranges to pixel ranges

            frac = (xyz_ll(outer) - xyz_ll_plot(outer)) / width(outer)
            outrange(1) = int(frac * real(pl % pixels(1), 8))
            frac = (xyz_ur(outer) - xyz_ll_plot(outer)) / width(outer)
            outrange(2) = int(frac * real(pl % pixels(1), 8))

            frac = (xyz_ur(inner) - xyz_ll_plot(inner)) / width(inner)
            inrange(1) = int((ONE - frac) * real(pl % pixels(2), 8))
            frac = (xyz_ll(inner) - xyz_ll_plot(inner)) / width(inner)
            inrange(2) = int((ONE - frac) * real(pl % pixels(2), 8))

            ! draw lines
            do out_ = outrange(1), outrange(2)
              do plus = 0, pl % meshlines_width
                data(:, out_ + 1, inrange(1) + plus + 1) = rgb
                data(:, out_ + 1, inrange(2) + plus + 1) = rgb
                data(:, out_ + 1, inrange(1) - plus + 1) = rgb
                data(:, out_ + 1, inrange(2) - plus + 1) = rgb
              end do
            end do
            do in_ = inrange(1), inrange(2)
              do plus = 0, pl % meshlines_width
                data(:, outrange(1) + plus + 1, in_ + 1) = rgb
                data(:, outrange(2) + plus + 1, in_ + 1) = rgb
                data(:, outrange(1) - plus + 1, in_ + 1) = rgb
                data(:, outrange(2) - plus + 1, in_ + 1) = rgb
              end do
            end do

          end if
        end do
      end do
    end associate

  end subroutine draw_mesh_lines

!===============================================================================
! OUTPUT_PPM writes out a previously generated image to a PPM file
!===============================================================================

  subroutine output_ppm(pl, data)
    type(ObjectPlot), intent(in) :: pl
    integer,          intent(in) :: data(:,:,:)

    integer :: y ! loop index for height
    integer :: x ! loop index for width
    integer :: unit_plot

    ! Open PPM file for writing
    open(NEWUNIT=unit_plot, FILE=pl % path_plot)

    ! Write header
    write(unit_plot, '(A2)') 'P6'
    write(unit_plot, '(I0,'' '',I0)') pl % pixels(1), pl % pixels(2)
    write(unit_plot, '(A)') '255'

    ! Write color for each pixel
    do y = 1, pl % pixels(2)
      do x = 1, pl % pixels(1)
        write(unit_plot, '(3A1)', advance='no') achar(data(RED, x, y)), &
             achar(data(GREEN, x, y)), achar(data(BLUE, x, y))
      end do
    end do

    ! Close plot file
    close(UNIT=unit_plot)
  end subroutine output_ppm

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

    integer :: x, y, z      ! voxel location indices
    integer :: rgb(3)       ! colors (red, green, blue) from 0-255
    integer :: id           ! id of cell or material
    integer :: hdf5_err
    integer, target :: data(pl%pixels(3),pl%pixels(2))
    integer(HID_T) :: file_id
    integer(HID_T) :: dspace
    integer(HID_T) :: memspace
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T) :: dims_slab(3)
    integer(HSIZE_T) :: offset(3)
    real(8) :: vox(3)       ! x, y, and z voxel widths
    real(8) :: ll(3)        ! lower left starting point for each sweep direction
    type(Particle)    :: p
    type(ProgressBar) :: progress
    type(c_ptr)       :: f_ptr

    ! compute voxel widths in each direction
    vox = pl % width/dble(pl % pixels)

    ! initial particle position
    ll = pl % origin - pl % width / TWO

    ! allocate and initialize particle
    call p % initialize()
    p % coord(1) % xyz = ll
    p % coord(1) % uvw = [ HALF, HALF, HALF ]
    p % coord(1) % universe = root_universe

    ! Open binary plot file for writing
    file_id = file_create(pl%path_plot)

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
    dims(:) = [pl%pixels(3), pl%pixels(2), pl%pixels(1)]
    call h5screate_simple_f(3, dims, dspace, hdf5_err)
    call h5dcreate_f(file_id, "data", H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Create another dataspace for 2D array in memory
    dims_slab(1) = pl%pixels(3)
    dims_slab(2) = pl%pixels(2)
    dims_slab(3) = 1
    call h5screate_simple_f(2, dims_slab(1:2), memspace, hdf5_err)

    ! Initialize offset and get pointer to data
    offset(:) = 0
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims_slab, hdf5_err)
    f_ptr = c_loc(data)

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
      offset(3) = x - 1
      call h5soffset_simple_f(dspace, offset, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, &
           mem_space_id=memspace, file_space_id=dspace)
    end do

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call file_close(file_id)

  end subroutine create_voxel

end module plot
