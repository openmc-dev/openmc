module plot

  use constants
  use error,           only: fatal_error
  use geometry
  use geometry_header
  use global
  use output,          only: write_message
  use particle_header, only: deallocate_coord, Particle
  use plot_header
  use ppmlib,          only: Image, init_image, allocate_image, &
                             deallocate_image, set_pixel
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
    type(Cell), pointer :: c => null()
    
    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0

    call find_cell(p, found_cell)
    if (check_overlaps) call check_cell_overlap(p)

    if (.not. found_cell) then
      ! If no cell, revert to default color
      rgb = pl % not_found % rgb
      id = -1
    else
      if (pl % color_by == PLOT_COLOR_MATS) then
        ! Assign color based on material
        c => cells(p % coord % cell)
        if (c % material == MATERIAL_VOID) then
          ! By default, color void cells white
          rgb = 255
          id = -1
        else
          rgb = pl % colors(c % material) % rgb
          id = materials(c % material) % id
        end if
      else if (pl % color_by == PLOT_COLOR_CELLS) then
        ! Assign color based on cell
        rgb = pl % colors(p % coord % cell) % rgb
        id = cells(p % coord % cell) % id
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
    logical :: found_cell
    integer :: cellid
    type(Image) :: img
    type(Particle) :: p
    type(Universe), pointer :: univ => null()


    cellid = -1


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
    call p % initialize()
    p % coord % xyz = xyz
    p % coord % uvw = (/ 0.5, 0.5, 0.5 /)
    p % coord % universe = BASE_UNIVERSE

    univ => universes(BASE_UNIVERSE)
    call deallocate_coord(p % coord % next)
    
    call calc_offsets(univ,cellid)
    call print_all(univ)
    do y = 1, img % height
      do x = 1, img % width

        ! get pixel color
        call position_rgb(p, pl, rgb, id)

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
    type(Particle) :: p

    ! compute voxel widths in each direction
    vox = pl % width/dble(pl % pixels)
    
    ! initial particle position
    ll = pl % origin - pl % width / 2.0

    ! allocate and initialize particle
    call p % initialize()
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

!===============================================================================
! CALC_OFFSETS calculates and stores the offsets in all cells
!===============================================================================

  recursive subroutine calc_offsets(univ, cellid)

    type(Universe), intent(in) :: univ
    integer, intent(in) :: cellid

    integer :: i                    ! index over cells
    integer :: i_x, i_y, i_z        ! indices in lattice
    integer :: n_x, n_y, n_z        ! size of lattice
    integer :: n                    ! number of cells to search
    integer :: prevoffset           ! total offset for the previous cell
    integer :: tempoffset           ! total offset for a given cell
    integer :: index_cell           ! index in cells array
    integer :: index_univ           ! index to next universe in universes array
    real(8) :: xyz(3)               ! temporary location
    real(8) :: upper_right(3)       ! lattice upper_right
    type(Cell),     pointer :: c => null()    ! pointer to cell
    type(Lattice),  pointer :: lat => null()  ! pointer to lattice
    type(Universe), pointer :: univ_next => null() ! next universe to loop through

    n = univ % n_cells
    
    print *, 'Uni ', univ % id , ' has ', n , ' cells'
    
    do i = 1, n
    
      
      index_cell = univ % cells(i)

      ! get pointer to cell
      c => cells(index_cell)
      tempoffset = 0
      call count_target_cell(c,cellid,tempoffset)
           
      
      if (i /= 1) then
      
        prevoffset = c % offset
      
      else
      
        prevoffset = 0
        c % offset = 0
        
      end if
       
      if (i /= n) then
      
        c => cells(index_cell + 1)
        c % offset = prevoffset + tempoffset
        c => cells(index_cell)
        
        print *, 'Cell ', c % id , ' has offset ' , c % offset
      
      else 
        
        print *, 'Cell ', c % id , ' has offset ' , c % offset
      
      end if
      
      
      
      !print *, 'Cell ', c % id
      
      if (c % type == CELL_NORMAL) then
        ! ====================================================================
        ! AT LOWEST UNIVERSE, TERMINATE SEARCH

      elseif (c % type == CELL_FILL) then
        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

        univ_next => universes(c % fill)
        call calc_offsets(univ_next,cellid)

      elseif (c % type == CELL_LATTICE) then
        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill)
        print *, 'Lattice ', lat % id


        n_x = lat % dimension(1)
        n_y = lat % dimension(2)
        if (lat % n_dimension == 3) then
          n_z = lat % dimension(3)
        else
          n_z = 1
        end if
        print *, 'Dimensions:',n_x,',',n_y,',',n_z
        ! Loop over lattice coordinates
        do i_x = 1, n_x
          do i_y = 1, n_y
            do i_z = 1, n_z

              univ_next => universes(lat % universes(i_x,i_y,i_z))
              
              tempoffset = 0
              if (i_x == 1 .AND. i_y == 1 .AND. i_z == 1) then
           
                lat % offset(1,1,1) = 0
                prevoffset = 0
                call count_target_univ(univ_next,cellid,tempoffset)

              else
              
                prevoffset = lat % offset(i_x,i_y,i_z)
                call count_target_univ(univ_next,cellid,tempoffset)
               
              end if
              

              if (i_z + 1 <= n_z) then
                lat % offset(i_x,i_y,i_z+1) = tempoffset + prevoffset              
              elseif (i_y + 1 <= n_y) then
                lat % offset(i_x,i_y+1,1) = tempoffset + prevoffset              
              elseif (i_x + 1 <= n_x) then
                lat % offset(i_x+1,1,1) = tempoffset + prevoffset
              end if
              
              
              print *, 'Lat',lat % id,' (',i_x,',',i_y,',',i_z,') has offset ' , lat % offset(i_x,i_y,i_z)
!              call calc_offsets(univ_next,cellid)

            end do
          end do
        end do

      end if
    end do
    
    return
    
  end subroutine calc_offsets
  
!===============================================================================
! COUNT_TARGET_CELL recursively totals the numbers of occurances of a given cell id
! beginning with the cell given. Using -1 for the cell id will total all 
! normal cells.
!===============================================================================

  recursive subroutine count_target_cell(c, cellid, kount)

    type(Cell), intent(in) :: c
    integer, intent(in) :: cellid
    integer, intent(inout) :: kount

    integer :: i                    ! index over cells
    integer :: i_x, i_y, i_z        ! indices in lattice
    integer :: n_x, n_y, n_z        ! size of lattice
    integer :: n                    ! number of cells to search
    integer :: index_cell           ! index in cells array
    integer :: index_univ           ! index to next universe in universes array
    real(8) :: xyz(3)               ! temporary location
    real(8) :: upper_right(3)       ! lattice upper_right
    type(Cell),     pointer, save :: c_next => null()    ! pointer to cell
    type(Lattice),  pointer, save :: lat => null()  ! pointer to lattice
    type(Universe), pointer, save :: univ_next => null() ! next universe to loop through

    if (c % type == CELL_NORMAL) then
      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (cellid == c % id .OR. cellid == -1) then
        kount = kount + 1
      endif

    elseif (c % type == CELL_FILL) then
      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

      univ_next => universes(c % fill)
      
      n = univ_next % n_cells

      !print *, 'Uni ', univ % id , ' has ', n , ' cells'
      
      do i = 1, n
      
        index_cell = univ_next % cells(i)

        ! get pointer to cell
        c_next => cells(index_cell)
        !print *, 'Cell ', c % id
        
        call count_target_cell(c_next, cellid, kount)
        
      end do

    elseif (c % type == CELL_LATTICE) then
      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

      ! Set current lattice
      lat => lattices(c % fill)


      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if
      
      ! Loop over lattice coordinates
      do i_x = 1, n_x
        do i_y = 1, n_y
          do i_z = 1, n_z

           ! print *, 'Lattice (',i_x,',',i_y,',',i_z,')'
!              index_univ = universe_dict % get_key()
            univ_next => universes(lat % universes(i_x,i_y,i_z))
            
            n = univ_next % n_cells
    
            !print *, 'Uni ', univ % id , ' has ', n , ' cells'
            
            do i = 1, n
            
              index_cell = univ_next % cells(i)

              ! get pointer to cell
              c_next => cells(index_cell)
              !print *, 'Cell ', c % id
              
              call count_target_cell(c_next, cellid, kount)
              
            end do
            

          end do
        end do
      end do

    end if
    
  end subroutine count_target_cell

!===============================================================================
! COUNT_TARGET_UNIV recursively totals the numbers of occurances of a given cell id
! beginning with the universe given. Using -1 for the cell id will total all 
! normal cells.
!===============================================================================

  recursive subroutine count_target_univ(univ, cellid, kount)

    type(Universe), intent(in) :: univ
    integer, intent(in) :: cellid
    integer, intent(inout) :: kount

    integer :: i                    ! index over cells
    integer :: i_x, i_y, i_z        ! indices in lattice
    integer :: n_x, n_y, n_z        ! size of lattice
    integer :: n                    ! number of cells to search
    integer :: index_cell           ! index in cells array
    integer :: index_univ           ! index to next universe in universes array
    real(8) :: xyz(3)               ! temporary location
    real(8) :: upper_right(3)       ! lattice upper_right
    type(Cell),     pointer, save :: c => null()    ! pointer to cell
    type(Lattice),  pointer, save :: lat => null()  ! pointer to lattice
    type(Universe), pointer, save :: univ_next => null() ! next universe to loop through

    n = univ % n_cells
    
    do i = 1, n
    
      
      index_cell = univ % cells(i)

      ! get pointer to cell
      c => cells(index_cell)
      
      !print *, 'Cell ', c % id
      
      if (c % type == CELL_NORMAL) then
        ! ====================================================================
        ! AT LOWEST UNIVERSE, TERMINATE SEARCH
        if (cellid == c % id .OR. cellid == -1) then
          kount = kount + 1
        endif
        
      elseif (c % type == CELL_FILL) then
        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

        univ_next => universes(c % fill)
        call count_target_univ(univ_next,cellid,kount)

      elseif (c % type == CELL_LATTICE) then
        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill)
        
        n_x = lat % dimension(1)
        n_y = lat % dimension(2)
        if (lat % n_dimension == 3) then
          n_z = lat % dimension(3)
        else
          n_z = 1
        end if
        
        ! Loop over lattice coordinates
        do i_x = 1, n_x
          do i_y = 1, n_y
            do i_z = 1, n_z

              univ_next => universes(lat % universes(i_x,i_y,i_z))
              call count_target_univ(univ_next,cellid,kount)

            end do
          end do
        end do

      end if
    end do
             
  end subroutine count_target_univ
  
!===============================================================================
! 
!===============================================================================

  recursive subroutine populate_cell_list(univ, celllist, ind, cellid)

    type(Universe), intent(in) :: univ
    integer, intent(inout) :: celllist(:)
    integer, intent(inout) :: ind
    integer, intent(inout) :: cellid

    integer :: i                    ! index over cells
    integer :: i_x, i_y, i_z        ! indices in lattice
    integer :: n_x, n_y, n_z        ! size of lattice
    integer :: n                    ! number of cells to search
    integer :: index_cell           ! index in cells array
    integer :: index_univ           ! index to next universe in universes array
    real(8) :: xyz(3)               ! temporary location
    real(8) :: upper_right(3)       ! lattice upper_right
    type(Cell),     pointer, save :: c => null()    ! pointer to cell
    type(Lattice),  pointer, save :: lat => null()  ! pointer to lattice
    type(Universe), pointer, save :: univ_next => null() ! next universe to loop through

    n = univ % n_cells
    
    do i = 1, n
    
      
      index_cell = univ % cells(i)

      ! get pointer to cell
      c => cells(index_cell)
      
      !print *, 'Cell ', c % id
      
      if (c % type == CELL_NORMAL) then
        ! ====================================================================
        ! AT LOWEST UNIVERSE, TERMINATE SEARCH
        if (cellid == c % id .OR. cellid == -1) then
          print  *,'c%id:',c%id
          print  *,'ind:',ind
          celllist(ind) = c % id
          ind = ind + 1          
        endif
        
      elseif (c % type == CELL_FILL) then
        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

        univ_next => universes(c % fill)
        call populate_cell_list(univ_next,celllist,ind,cellid)

      elseif (c % type == CELL_LATTICE) then
        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill)
        
        n_x = lat % dimension(1)
        n_y = lat % dimension(2)
        if (lat % n_dimension == 3) then
          n_z = lat % dimension(3)
        else
          n_z = 1
        end if
        
        ! Loop over lattice coordinates
        do i_x = 1, n_x
          do i_y = 1, n_y
            do i_z = 1, n_z

              univ_next => universes(lat % universes(i_x,i_y,i_z))
              call populate_cell_list(univ_next,celllist,ind,cellid)

            end do
          end do
        end do

      end if
    end do
             
  end subroutine populate_cell_list

!===============================================================================
! PRINT_ALL prints the tally array index and coordinate chain for 
! every normal cell.
!===============================================================================

  subroutine print_all(univ)

    type(Universe), intent(in) :: univ

    integer :: i                     ! index over cells
    integer :: i_x, i_y, i_z         ! indices in lattice
    integer :: n_x, n_y, n_z         ! size of lattice
    integer :: n                     ! number of cells to search
    integer :: index_cell            ! index in cells array
    integer :: index_univ            ! index to next universe in universes array
    integer :: kount                 ! number of normal cells in the system
    integer :: cellid                ! cellid to use. will be set to -1 to signal all cells
    integer, allocatable :: celllist(:) ! vector of normal cells in system
    real(8) :: xyz(3)                ! temporary location
    real(8) :: upper_right(3)        ! lattice upper_right
    type(Cell),     pointer, save :: c => null()    ! pointer to cell
    type(Lattice),  pointer, save :: lat => null()  ! pointer to lattice
    type(Universe), pointer, save :: univ_next => null() ! next universe to loop through


    cellid = -1
    i = 1
    call count_target_univ(univ, cellid, kount)
    allocate(celllist(kount))
    call populate_cell_list(univ,celllist,i,cellid)
    print *,kount
    do i = 1, kount
    
      !print  *,celllist(i)
      
    end do
             
  deallocate(celllist)
  end subroutine print_all

  
end module plot
