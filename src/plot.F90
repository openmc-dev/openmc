module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use mesh,            only: get_mesh_bin, mesh_indices_to_bin
  use output,          only: write_message
  use particle_header, only: deallocate_coord
  use plot_header
  use ppmlib,          only: Image, init_image, allocate_image, &
                             deallocate_image, set_pixel
  use progress_header
  use string,          only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! RUN_PLOTs controls the logic for making one or many plots in MODE_PLOTTING
!===============================================================================

  subroutine run_plots()

    integer :: i ! loop index for plots
    type(ObjectPlot), pointer :: pl => null()

    do i = 1, n_plots
      pl => plots(i)

      ! Display output message
      message = "Processing plot " // trim(to_str(pl % id)) // "..."
      call write_message(5)

      if (pl % type == PLOT_TYPE_SLICE) then
        ! create 2d image
        call create_ppm(pl)
      else if (pl % type == PLOT_TYPE_VOXEL) then
        ! create dump for 3D silomesh utility script
        call create_3d_dump(pl)
      end if
    end do

  end subroutine run_plots

!===============================================================================
! SCORE_RXN_RATE_PLOTS
!===============================================================================

  subroutine score_rxn_rate_plots()

    integer :: i                             ! loop index for plots
    integer :: n
    integer :: bin
    type(ObjectPlot), pointer     :: pl => null()
    type(StructuredMesh), pointer :: m => null()

    do i = 1, n_plots
      pl => plots(i)

      if (.not. pl % type == PLOT_TYPE_RXNRATE) cycle

      if (.not. nuclides(p % event_nuclide) % fissionable) then
        select case(pl % rrtype)
          case (PLOT_RXN_FLUX_THERMAL)
            if (.not. p % E <= 0.625e-6) cycle
          case (PLOT_RXN_FLUX_FAST)
            if (.not. p % E >= 0.625e-6) cycle
        end select
      end if

      m => pl % pixmesh
      call get_mesh_bin(m, p % coord0 % xyz, bin)

      if (bin == NO_BIN_FOUND) cycle

      if (nuclides(p % event_nuclide) % fissionable) then
        n = nuclides(p % event_nuclide) % index_fission(1)
        pl % fisswgt(bin) = pl% fisswgt(bin) + p % last_wgt * &
                nuclides(p % event_nuclide) % reactions(n) % Q_value
      else
        pl % fluxwgt(bin) = pl % fluxwgt(bin) + p % last_wgt / material_xs % total
      end if
      
    end do

  end subroutine score_rxn_rate_plots

!===============================================================================
! FINALIZE_RXN_PLOTS
!===============================================================================

  subroutine finalize_rxn_plots()

    integer :: i                             ! loop index for plots
    type(ObjectPlot), pointer     :: pl => null()

    do i = 1, n_plots
      pl => plots(i)

      if (.not. pl % type == PLOT_TYPE_RXNRATE) cycle

#ifdef MPI
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, pl % fisswgt, size(pl % fisswgt), &
             MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(MPI_IN_PLACE, pl % fluxwgt, size(pl % fluxwgt), &
             MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        call MPI_REDUCE(pl % fisswgt, pl % fisswgt, size(pl % fisswgt), &
             MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_REDUCE(pl % fluxwgt, pl % fluxwgt, size(pl % fluxwgt), &
             MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
        cycle
      end if
#endif

      pl % fisswgt = pl % fisswgt / maxval(pl % fisswgt)
      pl % fluxwgt = pl % fluxwgt / maxval(pl % fluxwgt)

      ! Display output message
      message = "Processing rxn plot " // trim(to_str(pl % id)) // "..."
      call write_message(5)

      call create_ppm(pl, rxn_overlay=.true.)

    end do

  end subroutine finalize_rxn_plots

!===============================================================================
! POSITION_RGB computes the red/green/blue values for a given plot with the 
! current particle's position
!===============================================================================

  subroutine position_rgb(pl, rgb, id, fissionable)
  
    type(ObjectPlot), pointer, intent(in) :: pl
    integer, intent(out)                  :: rgb(3)
    integer, intent(out)                  :: id
    logical, intent(out), optional        :: fissionable
    
    logical :: found_cell
    type(Cell), pointer :: c => null()
    
    if (present(fissionable)) fissionable = .false.

    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0

    call find_cell(found_cell)

    if (.not. found_cell) then
      ! If no cell, revert to default color
      rgb = pl % not_found % rgb
      id = -1
    else
      c => cells(p % coord % cell)
      if (pl % color_by == PLOT_COLOR_MATS) then
        ! Assign color based on material
        id = materials(c % material) % id
        if (c % material == MATERIAL_VOID) then
          ! By default, color void cells white
          rgb = 255
        else
          rgb = pl % colors(c % material) % rgb
        end if
      else if (pl % color_by == PLOT_COLOR_CELLS) then
        ! Assign color based on cell
        rgb = pl % colors(p % coord % cell) % rgb
        id = cells(p % coord % cell) % id
      else
        rgb = 0
        id = -1
      end if
      if (present(fissionable)) then
        fissionable = materials(c % material) % fissionable
      end if
    end if
    
  end subroutine position_rgb

!===============================================================================
! CREATE_PPM creates an image based on user input from a plots.xml <plot>
! specification in the portable pixmap format (PPM)
!===============================================================================

  subroutine create_ppm(pl, rxn_overlay)

    type(ObjectPlot), pointer :: pl
    logical, optional         :: rxn_overlay

    integer :: in_i
    integer :: out_i
    integer :: x, y      ! pixel location
    integer :: ijk(2)    ! pixmesh indices
    integer :: rgb(3)    ! colors (red, green, blue) from 0-255
    integer :: id
    integer :: bin       ! pixmesh results bin
    real(8) :: rxnval    ! pixmesh results value
    real(8) :: in_pixel
    real(8) :: out_pixel
    real(8) :: xyz(3)
    logical :: rxn = .false.
    logical :: fissionable
    type(Image) :: img
    type(ProgressBar) :: progress
    type(StructuredMesh), pointer :: m => null()

    if (present(rxn_overlay) .and. rxn_overlay) then
      rxn = .true.
      m => pl % pixmesh
    end if

    ! Initialize and allocate space for image
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
    ! We can't depend on source.F90 for initialize_particle due to
    ! circular dependencies (since physics.F90 depends on this file)
    ! To fix this, transport and physics could be 
    allocate(p)
    allocate(p % coord0)
    p % coord0 % universe = BASE_UNIVERSE
    p % coord => p % coord0
    p % coord % xyz = xyz
    p % coord % uvw = (/ 0.5, 0.5, 0.5 /)
    p % coord % universe = BASE_UNIVERSE

    do y = 1, img % height
      call progress % set_value(dble(y)/dble(img % height)*100.)
      do x = 1, img % width

        ! Get pixel color
        call position_rgb(pl, rgb, id, fissionable)

        ! Add rxn rate overlay if appropriate
        if (rxn) then
          ijk(1) = x
          ijk(2) = y
          bin = mesh_indices_to_bin(m, ijk)
          if (fissionable) then
            rxnval = pl % fisswgt(bin)
            rgb(1) = min(int(1010.*rxnval),255)
            rgb(2) = int(255.*rxnval)
            rgb(3) = int(255.*rxnval)
          else
            rxnval = pl % fluxwgt(bin)
            rgb(1) = int(255.*rxnval)
            rgb(2) = int(255.*rxnval)
            rgb(3) = min(int(1010.*rxnval),255)
          end if
        end if

        ! Create a pixel at (x,y) with color (r,g,b)
        call set_pixel(img, x, y, rgb(1), rgb(2), rgb(3))

        ! Advance pixel in first direction
        p % coord0 % xyz(in_i) = p % coord0 % xyz(in_i) + in_pixel
      end do

      ! Advance pixel in second direction
      p % coord0 % xyz(in_i)  = xyz(in_i)
      p % coord0 % xyz(out_i) = p % coord0 % xyz(out_i) - out_pixel
    end do

    ! Write out the ppm to a file
    call output_ppm(pl,img)

    ! Free up space
    call deallocate_image(img)

  end subroutine create_ppm

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
    type(ProgressBar) :: progress

    ! compute voxel widths in each direction
    vox = pl % width/dble(pl % pixels)
    
    ! initial particle position
    ll = pl % origin - pl % width / 2.0

    ! allocate and initialize particle
    ! We can't depend on source.F90 for initialize_particle due to
    ! circular dependencies (since physics.F90 depends on this file)
    ! To fix this, transport and physics could be 
    allocate(p)
    allocate(p % coord0)
    p % coord0 % universe = BASE_UNIVERSE
    p % coord => p % coord0
    p % coord0 % xyz = ll
    p % coord0 % uvw = (/ 0.5, 0.5, 0.5 /)
    p % coord0 % universe = BASE_UNIVERSE

    ! Open binary plot file for writing
    open(UNIT=UNIT_PLOT, FILE=pl % path_plot, STATUS='replace', &
         ACCESS='stream')

    ! write plot header info
    write(UNIT_PLOT) pl % pixels, vox, ll

    ! move to center of voxels    
    ll = ll + vox / 2.0

    do x = 1, pl % pixels(1)
      call progress % set_value(dble(x)/dble(pl % pixels(1))*100.)
      do y = 1, pl % pixels(2)
        do z = 1, pl % pixels(3)

          ! get voxel color
          call position_rgb(pl, rgb, id)

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
