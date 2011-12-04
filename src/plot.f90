module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, cell_contains
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use particle_header, only: Particle, initialize_particle, LocalCoord,      &
                             deallocate_coord

  implicit none

contains

!===============================================================================
! RUN_PLOT generates a binary stream file containing a list of surface/lattice
! crossings and what cell was traveled through. A Python script can then be used
! to generate a plot based on the recorded crossings and cells
!===============================================================================

  subroutine run_plot()

    integer :: i               ! loop index
    integer :: surface_crossed ! surface which particle is on
    integer :: last_cell       ! most recent cell particle was in
    real(8) :: xyz(3)        ! starting coordinates
    real(8) :: last_x_coord    ! bounding x coordinate
    real(8) :: last_y_coord    ! bounding y coordinate
    real(8) :: d               ! distance to boundary
    real(8) :: distance        ! distance particle travels
    logical :: found_cell      ! found cell which particle is in?
    logical :: lattice_crossed ! is surface crossing in lattice?
    character(MAX_LINE_LEN) :: path_plot ! unit for binary plot file
    type(Cell),       pointer :: c    => null()
    type(Universe),   pointer :: univ => null()
    type(Particle),   pointer :: p    => null()
    type(LocalCoord), pointer :: coord => null()

    ! Open plot file for binary writing
    path_plot = trim(path_input) // "plot.out"
    open(UNIT=UNIT_PLOT, FILE=path_plot, STATUS="replace", ACCESS="stream")

    ! Write origin, width, basis, and pixel width to file
    write(UNIT=UNIT_PLOT) plot_origin
    write(UNIT=UNIT_PLOT) plot_width
    write(UNIT=UNIT_PLOT) plot_basis
    write(UNIT=UNIT_PLOT) pixel

    ! Determine coordinates of the upper-left corner of the plot
    xyz(1) = plot_origin(1) - plot_width(1) / 2.0
    xyz(2) = plot_origin(2) + (plot_width(2) - pixel) / 2.0
    xyz(3) = plot_origin(3)

    ! Determine bounding x and y coordinates for plot
    last_x_coord = plot_origin(1) + plot_width(1) / 2.0
    last_y_coord = plot_origin(2) - plot_width(2) / 2.0

    ! allocate and initialize particle
    allocate(p)

    ! loop over horizontal rays
    do while(xyz(2) > last_y_coord)

       ! initialize the particle and set starting coordinate and direction
       call initialize_particle(p)

       p % coord % xyz = xyz
       p % coord % uvw = (/ 1, 0, 0 /)

       ! write starting coordinate to file
       write(UNIT=UNIT_PLOT) p % coord % xyz

       ! Find cell that particle is currently in
       call find_cell(p, found_cell)

       ! =======================================================================
       ! MOVE PARTICLE FORWARD TO NEXT CELL

       if (.not. found_cell) then
          ! Clear any coordinates beyond first level
          call deallocate_coord(p % coord0 % next)
          p % coord => p % coord0

          univ => universes(BASE_UNIVERSE)
          do i = 1, univ % n_cells
             p % coord0 % xyz = xyz
             p % coord0 % cell = univ % cells(i)

             distance = INFINITY
             call distance_to_boundary(p, d, surface_crossed, lattice_crossed)
             if (d < distance) then
                ! Move particle forward to next surface
                ! Advance particle
                p % coord0 % xyz = p % coord0 % xyz + d * p % coord0 % uvw

                ! Check to make sure particle is actually going into this cell
                ! by moving it slightly forward and seeing if the cell contains
                ! that coordinate

                p % coord0 % xyz = p % coord0 % xyz + 1e-4 * p % coord0 % uvw

                c => cells(p % coord0 % cell)
                if (.not. cell_contains(c, p)) cycle

                ! Reset coordinate to surface crossing
                p % coord0 % xyz = p % coord0 % xyz - 1e-4 * p % coord0 % uvw

                ! Set new distance and retain pointer to this cell
                distance = d
                last_cell = p % coord0 % cell
             end if
          end do

          ! No cell was found on this horizontal ray
          if (distance == INFINITY) then
             p % coord0 % xyz(1) = last_x_coord
             write(UNIT_PLOT) p % coord0 % xyz, 0

             ! Move to next horizontal ray
             xyz(2) = xyz(2) - pixel
             cycle
          end if

          ! Write coordinate where next cell begins
          write(UNIT=UNIT_PLOT) p % coord0 % xyz, 0

          ! Process surface crossing for next cell
          p % coord0 % cell = NONE
          p % surface = -surface_crossed
          call cross_surface(p, last_cell)
       end if

       ! =======================================================================
       ! MOVE PARTICLE ACROSS HORIZONTAL TRACK

       do while (p % alive)
          ! save particle's current cell
          last_cell = p % coord % cell

          ! Calculate distance to next boundary
          call distance_to_boundary(p, distance, surface_crossed, lattice_crossed)

          ! Advance particle
          coord => p % coord0
          do while (associated(coord))
             coord % xyz = coord % xyz + distance * coord % uvw
             coord => coord % next
          end do

          ! If next boundary crossing is out of range of the plot, only include
          ! the visible portion and move to next horizontal ray
          if (p % coord0 % xyz(1) >= last_x_coord) then
             p % alive = .false.
             p % coord0 % xyz(1) = last_x_coord

             ! If there is no cell beyond this boundary, mark it as cell 0
             if (distance == INFINITY) p % coord % cell = 0

             ! Write ending coordinates to file
             write(UNIT=UNIT_PLOT) p % coord0 % xyz, last_cell
             cycle
          end if

          ! Write boundary crossing coordinates to file
          write(UNIT=UNIT_PLOT) p % coord0 % xyz, last_cell

          p % coord % cell = 0
          if (lattice_crossed) then
             p % surface = NONE
             call cross_lattice(p)
          else
             p % surface = surface_crossed
             call cross_surface(p, last_cell)

             ! Since boundary conditions are disabled in plotting mode, we need
             ! to manually add the last segment
             if (surfaces(abs(surface_crossed)) % bc == BC_VACUUM) then
                p % coord0 % xyz(1) = last_x_coord
                write(UNIT=UNIT_PLOT) p % coord0 % xyz, 0
                exit
             end if
          end if

       end do

       ! Move y-coordinate to next position
       xyz(2) = xyz(2) - pixel
    end do

    ! Close plot file
    close(UNIT=UNIT_PLOT)

  end subroutine run_plot

end module plot
