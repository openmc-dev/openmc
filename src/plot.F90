module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, check_cell_overlap
  use geometry_header, only: Cell, BASE_UNIVERSE
  use global
  use mesh,            only: get_mesh_indices
  use output,          only: write_message
  use particle_header, only: deallocate_coord, Particle, LocalCoord
  use plot_header
  use ppmlib,          only: Image, init_image, allocate_image, &
                             deallocate_image, set_pixel
  use progress_header, only: ProgressBar
  use string,          only: to_str

  implicit none

contains

!===============================================================================
! RUN_PLOT controls the logic for making one or many plots
!===============================================================================

  subroutine run_plot()

    integer :: i ! loop index for plots
    type(ObjectPlot), pointer :: pl => null()

    do i = 1, n_plots
      pl => plots(i)

      ! Display output message
      call write_message("Processing plot " // trim(to_str(pl % id)) &
           &// ": " // trim(pl % path_plot) // " ...", 5)

      if (pl % type == PLOT_TYPE_SLICE) then
        ! create 2d image
        call create_ppm(pl)
      else if (pl % type == PLOT_TYPE_VOXEL) then
        ! create dump for 3D silomesh utility script
        call create_3d_dump(pl)
      end if
    end do

  end subroutine run_plot

!===============================================================================
! POSITION_RGB computes the red/green/blue values for a given plot with the
! current particle's position
!===============================================================================

  subroutine position_rgb(p, pl, rgb, id)

    type(Particle), intent(inout)         :: p
    type(ObjectPlot), pointer, intent(in) :: pl
    integer, intent(out)                  :: rgb(3)
    integer, intent(out)                  :: id

    logical :: found_cell
    integer :: level
    type(Cell),       pointer :: c => null()
    type(LocalCoord), pointer :: coord => null()

    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0

    call find_cell(p, found_cell)
    if (check_overlaps) call check_cell_overlap(p)

    ! Loop through universes and stop on any specified level
    level = 0
    coord => p % coord0
    do
      if (level == pl % level) exit
      if (.not. associated(coord % next)) exit
      coord => coord % next
      level = level + 1
    end do

    if (.not. found_cell) then
      ! If no cell, revert to default color
      rgb = pl % not_found % rgb
      id = -1
    else
      if (pl % color_by == PLOT_COLOR_MATS) then
        ! Assign color based on material
        c => cells(coord % cell)
        if (c % material == MATERIAL_VOID) then
          ! By default, color void cells white
          rgb = 255
          id = -1
        else if (c % type == CELL_FILL) then
          ! If we stopped on a middle universe level, treat as if not found
          rgb = pl % not_found % rgb
          id = -1
        else
          rgb = pl % colors(c % material) % rgb
          id = materials(c % material) % id
        end if
      else if (pl % color_by == PLOT_COLOR_CELLS) then
        ! Assign color based on cell
        rgb = pl % colors(coord % cell) % rgb
        id = cells(coord % cell) % id
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

    type(ObjectPlot), pointer :: pl

    integer :: in_i
    integer :: out_i
    integer :: x, y      ! pixel location
    integer :: rgb(3)    ! colors (red, green, blue) from 0-255
    integer :: id
    real(8) :: in_pixel
    real(8) :: out_pixel
    real(8) :: xyz(3)
    type(Image)       :: img
    type(Particle)    :: p
    type(ProgressBar) :: progress

    ! Initialize and allocate space for image
    call init_image(img)
    call allocate_image(img, pl % pixels(1), pl % pixels(2))

    in_pixel  = pl % width(1)/dble(pl % pixels(1))
    out_pixel = pl % width(2)/dble(pl % pixels(2))

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
    p % coord % xyz = xyz
    p % coord % uvw = [ HALF, HALF, HALF ]
    p % coord % universe = BASE_UNIVERSE

    do y = 1, img % height
      call progress % set_value(dble(y)/dble(img % height)*100)
      do x = 1, img % width

        ! get pixel color
        call position_rgb(p, pl, rgb, id)

        ! Create a pixel at (x,y) with color (r,g,b)
        call set_pixel(img, x-1, y-1, rgb(1), rgb(2), rgb(3))

        ! Advance pixel in first direction
        p % coord0 % xyz(in_i) = p % coord0 % xyz(in_i) + in_pixel
      end do

      ! Advance pixel in second direction
      p % coord0 % xyz(in_i)  = xyz(in_i)
      p % coord0 % xyz(out_i) = p % coord0 % xyz(out_i) - out_pixel
    end do

    ! Draw tally mesh boundaries on the image if requested
    if (associated(pl % meshlines_mesh)) call draw_mesh_lines(pl, img)

    ! Write out the ppm to a file
    call output_ppm(pl,img)

    ! Free up space
    call deallocate_image(img)

    ! Clear particle
    call p % clear()

  end subroutine create_ppm

!===============================================================================
! DRAW_MESH_LINES draws mesh line boundaries on an image
!===============================================================================
  subroutine draw_mesh_lines(pl, img)

    type(ObjectPlot), pointer, intent(in)    :: pl
    type(Image),               intent(inout) :: img

    logical :: in_mesh
    integer :: out_, in_  ! pixel location
    integer :: r, g, b    ! RGB color for meshlines pixels
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
    type(StructuredMesh), pointer :: m => null()

    m => pl % meshlines_mesh

    r = pl % meshlines_color % rgb(1)
    g = pl % meshlines_color % rgb(2)
    b = pl % meshlines_color % rgb(3)

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

    xyz_ll_plot(outer) = pl % origin(1) - pl % width(1) / TWO
    xyz_ll_plot(inner) = pl % origin(2) - pl % width(2) / TWO
    xyz_ur_plot(outer) = pl % origin(1) + pl % width(1) / TWO
    xyz_ur_plot(inner) = pl % origin(2) + pl % width(2) / TWO

    width = xyz_ur_plot - xyz_ll_plot

    call get_mesh_indices(m, xyz_ll_plot, ijk_ll(:m % n_dimension), in_mesh)
    call get_mesh_indices(m, xyz_ur_plot, ijk_ur(:m % n_dimension), in_mesh)

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
          outrange(1) = int(frac * real(img % width, 8))
          frac = (xyz_ur(outer) - xyz_ll_plot(outer)) / width(outer)
          outrange(2) = int(frac * real(img % width, 8))

          frac = (xyz_ur(inner) - xyz_ll_plot(inner)) / width(inner)
          inrange(1) = int((ONE - frac) * real(img % height, 8))
          frac = (xyz_ll(inner) - xyz_ll_plot(inner)) / width(inner)
          inrange(2) = int((ONE - frac) * real(img % height, 8))

          ! draw lines
          do out_ = outrange(1), outrange(2)
            do plus = 0, pl % meshlines_width
              call set_pixel(img, out_, inrange(1) + plus, r, g, b)
              call set_pixel(img, out_, inrange(2) + plus, r, g, b)
              call set_pixel(img, out_, inrange(1) - plus, r, g, b)
              call set_pixel(img, out_, inrange(2) - plus, r, g, b)
            end do
          end do
          do in_ = inrange(1), inrange(2)
            do plus = 0, pl % meshlines_width
              call set_pixel(img, outrange(1) + plus, in_, r, g, b)
              call set_pixel(img, outrange(2) + plus, in_, r, g, b)
              call set_pixel(img, outrange(1) - plus, in_, r, g, b)
              call set_pixel(img, outrange(2) - plus, in_, r, g, b)
            end do
          end do

        end if
      end do
    end do

  end subroutine draw_mesh_lines

!===============================================================================
! OUTPUT_PPM writes out a previously generated image to a PPM file
!===============================================================================

  subroutine output_ppm(pl, img)

    type(ObjectPlot), pointer :: pl
    type(Image), intent(in)  :: img

    integer :: i ! loop index for height
    integer :: j ! loop index for width

    ! Open PPM file for writing
    open(UNIT=UNIT_PLOT, FILE=pl % path_plot)

    ! Write header
    write(UNIT_PLOT, '(A2)') 'P6'
    write(UNIT_PLOT, '(I0,'' '',I0)') img%width, img%height
    write(UNIT_PLOT, '(A)') '255'

    ! Write color for each pixel
    do j = 1, img % height
      do i = 1, img % width
        write(UNIT_PLOT, '(3A1)', advance='no') achar(img%red(i,j)), &
             achar(img%green(i,j)), achar(img%blue(i,j))
      end do
    end do

    ! Close plot file
    close(UNIT=UNIT_PLOT)

  end subroutine output_ppm

!===============================================================================
! CREATE_3D_DUMP outputs a binary file that can be input into silomesh for 3D
! geometry visualization.  It works the same way as create_ppm by dragging a
! particle across the geometry for the specified number of voxels. The first
! 3 int(4)'s in the binary are the number of x, y, and z voxels.  The next 3
! real(8)'s are the widths of the voxels in the x, y, and z directions. The next
! 3 real(8)'s are the x, y, and z coordinates of the lower left point. Finally
! the binary is filled with entries of four int(4)'s each. Each 'row' in the
! binary contains four int(4)'s: 3 for x,y,z position and 1 for cell or material
! id.  For 1 million voxels this produces a file of approximately 15MB.
!===============================================================================

  subroutine create_3d_dump(pl)

    type(ObjectPlot), pointer :: pl

    integer :: x, y, z      ! voxel location indices
    integer :: rgb(3)       ! colors (red, green, blue) from 0-255
    integer :: id           ! id of cell or material
    real(8) :: vox(3)       ! x, y, and z voxel widths
    real(8) :: ll(3)        ! lower left starting point for each sweep direction
    type(Particle)    :: p
    type(ProgressBar) :: progress

    ! compute voxel widths in each direction
    vox = pl % width/dble(pl % pixels)

    ! initial particle position
    ll = pl % origin - pl % width / TWO

    ! allocate and initialize particle
    call p % initialize()
    p % coord0 % xyz = ll
    p % coord0 % uvw = [ HALF, HALF, HALF ]
    p % coord0 % universe = BASE_UNIVERSE

    ! Open binary plot file for writing
    open(UNIT=UNIT_PLOT, FILE=pl % path_plot, STATUS='replace', &
         ACCESS='stream')

    ! write plot header info
    write(UNIT_PLOT) pl % pixels, vox, ll

    ! move to center of voxels
    ll = ll + vox / TWO

    do x = 1, pl % pixels(1)
      call progress % set_value(dble(x)/dble(pl % pixels(1))*100)
      do y = 1, pl % pixels(2)
        do z = 1, pl % pixels(3)

          ! get voxel color
          call position_rgb(p, pl, rgb, id)

          ! write to plot file
          write(UNIT_PLOT) id

          ! advance particle in z direction
          p % coord0 % xyz(3) = p % coord0 % xyz(3) + vox(3)

        end do

        ! advance particle in y direction
        p % coord0 % xyz(2) = p % coord0 % xyz(2) + vox(2)
        p % coord0 % xyz(3) = ll(3)

      end do

      ! advance particle in y direction
      p % coord0 % xyz(1) = p % coord0 % xyz(1) + vox(1)
      p % coord0 % xyz(2) = ll(2)
      p % coord0 % xyz(3) = ll(3)

    end do

    close(UNIT_PLOT)

  end subroutine create_3d_dump

end module plot
