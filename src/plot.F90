module plotter

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, cell_contains
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use output,          only: write_message
  use particle_header, only: LocalCoord, deallocate_coord
  use plot_header
  use ppmlib,          only: Image, init_image, allocate_image, &
                             deallocate_image, set_pixel
  use source,          only: initialize_particle
  use string,          only: to_str

  implicit none

contains


!===============================================================================
! RUN_PLOT
!===============================================================================

  subroutine run_plot()

    integer                  :: i
    type(Plot),    pointer   :: pl => null()

    do i=1,n_plots
      pl => plots(i)

      ! Display output message
      message = "Processing plot " // trim(to_str(pl % id)) // "..."
      call write_message(5)

      if (pl % type == PLOT_TYPE_SLICE) then

        ! create 2d image
        call create_ppm(pl)

      end if
    end do

  end subroutine run_plot

!===============================================================================
! create_ppm
!===============================================================================
  subroutine create_ppm(pl)

    type(Plot),         pointer :: pl

    type(Image)         :: img
    integer             :: in_i, out_i
    integer             :: x, y
    integer             :: r, g, b
    real(8)             :: in_pixel, out_pixel
    real(8)             :: xyz(3)
    logical             :: found_cell
    type(Cell), pointer :: c    => null()

    call init_image(img)
    call allocate_image(img, pl % pixels(1), pl % pixels(2))

    if (pl % basis == PLOT_BASIS_XY) then
      in_i  = 1
      out_i = 2
      in_pixel  = pl % width(1)/dble(pl % pixels(1))
      out_pixel = pl % width(2)/dble(pl % pixels(2))
      xyz(1) = pl % origin(1) - pl % width(1) / 2.0
      xyz(2) = pl % origin(2) + pl % width(2) / 2.0
      xyz(3) = pl % origin(3)
    else if (pl % basis == PLOT_BASIS_XZ) then
      in_i  = 1
      out_i = 3
      in_pixel  = pl % width(1)/dble(pl % pixels(1))
      out_pixel = pl % width(2)/dble(pl % pixels(2))
      xyz(1) = pl % origin(1) - pl % width(1) / 2.0
      xyz(2) = pl % origin(2)
      xyz(3) = pl % origin(3) + pl % width(2) / 2.0
    else if (pl % basis == PLOT_BASIS_YZ) then
      in_i  = 2
      out_i = 3
      in_pixel  = pl % width(1)/dble(pl % pixels(1))
      out_pixel = pl % width(2)/dble(pl % pixels(2))
      xyz(1) = pl % origin(1)
      xyz(2) = pl % origin(2) - pl % width(1) / 2.0
      xyz(3) = pl % origin(3) + pl % width(2) / 2.0
    end if

    ! allocate and initialize particle
    allocate(p)
    call initialize_particle()
    p % coord % xyz = xyz
    p % coord % uvw = (/ 1, 0, 0 /)
    p % coord % universe = BASE_UNIVERSE

    do y=1, img%height
      do x=1, img%width

        call deallocate_coord(p % coord0 % next)
        p % coord => p % coord0

        call find_cell(found_cell)

        if (.not. found_cell) then
          r = pl % not_found % rgb(1)
          g = pl % not_found % rgb(2)
          b = pl % not_found % rgb(3)
        else
          if (pl % color_by == PLOT_COLOR_MATS) then
            c => cells(p % coord % cell)
            r = pl % colors(c % material) % rgb(1)
            g = pl % colors(c % material) % rgb(2)
            b = pl % colors(c % material) % rgb(3)
          else if (pl % color_by == PLOT_COLOR_CELLS) then
            r = pl % colors(p % coord % cell) % rgb(1)
            g = pl % colors(p % coord % cell) % rgb(2)
            b = pl % colors(p % coord % cell) % rgb(3)
          else
            r = 0
            g = 0
            b = 0
          end if
        end if

        call set_pixel(img, x, y, r, g, b)

        p % coord0 % xyz(in_i) = p % coord0 % xyz(in_i) + in_pixel
      end do

      p % coord0 % xyz(in_i)  = xyz(in_i)
      p % coord0 % xyz(out_i) = p % coord0 % xyz(out_i) - out_pixel
    end do

    call output_ppm(pl,img)

    call deallocate_image(img)

  end subroutine create_ppm


!===============================================================================
! output_ppm
!===============================================================================
  subroutine output_ppm(pl,img)

    type(Plot),        pointer :: pl
    type(Image),    intent(in) :: img

    integer :: i, j
    character(MAX_LINE_LEN)    :: path_plot ! unit for binary plot file

    path_plot = trim(path_input) // "slice" // trim(to_str(pl % id)) // ".ppm"
    open(UNIT=UNIT_PLOT, FILE=path_plot)

    write(UNIT_PLOT, '(A2)') 'P6'
    write(UNIT_PLOT, '(I0,'' '',I0)') img%width, img%height
    write(UNIT_PLOT, '(A)') '255'
 
    do j=1, img%height
       do i=1, img%width
          write(UNIT_PLOT, '(3A1)', advance='no') achar(img%red(i,j)), &
                                                  achar(img%green(i,j)), &
                                                  achar(img%blue(i,j))
       end do
    end do

    ! Close plot file
    close(UNIT=UNIT_PLOT)
 
  end subroutine output_ppm


end module plotter
