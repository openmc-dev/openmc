module fileio

  use global
  use types,  only: Cell, Surface, ExtSource
  use string, only: split_string_wl, lower_case, str_to_real, split_string, &
       &            concatenate
  use output, only: message, warning, error
  use data_structures, only: dict_create, dict_add_key, dict_get_key, &
       &                     DICT_NULL

  implicit none

!=====================================================================
! READ_DATA interface allows data to be read with one function
! regardless of whether it is integer or real data. E.g. NXS and JXS
! can be read with the integer version and XSS can be read with the
! real version
!=====================================================================

  interface read_data
     module procedure read_data_int, read_data_real
  end interface

contains

!=====================================================================
! READ_COMMAND_LINE reads all parameters from the command line
!=====================================================================

  subroutine read_command_line()

    integer :: argc
    character(250) :: msg
    logical :: file_exists
    character(7) :: readable

    argc = COMMAND_ARGUMENT_COUNT()
    if ( argc > 0 ) then
       call GET_COMMAND_ARGUMENT(1, path_input)
    else
       msg = "No input file specified!"
       call error(msg)
    end if

    ! Check if input file exists and is readable
    inquire(FILE=path_input, EXIST=file_exists, READ=readable)
    if (.not. file_exists) then
       msg = "Input file '" // trim(path_input) // "' does not exist!"
       call error(msg)
    elseif (readable(1:3) == 'NO') then
       ! Need to explicitly check for a NO status -- Intel compiler
       ! looks at file attributes only if the file is open on a
       ! unit. Therefore, it will always return UNKNOWN
       msg = "Input file '" // trim(path_input) // "' is not readable! &
            &Change file permissions with chmod command."
       call error(msg)
    end if

  end subroutine read_command_line

!=====================================================================
! READ_COUNT makes a first pass through the input file to determine
! how many cells, surfaces, materials, sources, etc there are in the
! problem.
!=====================================================================

  subroutine read_count(filename)
    
    character(*), intent(in) :: filename

    character(250) :: line, msg
    character(32) :: words(max_words)
    integer :: in = 7
    integer :: ioError
    integer :: n

    msg = "First pass through input file..."
    call message(msg, 5)

    n_cells     = 0
    n_surfaces  = 0
    n_materials = 0

    open(FILE=filename, UNIT=in, STATUS='old', & 
         & ACTION='read', IOSTAT=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call error(msg)
    end if

    do
       call get_next_line(in, words, n, ioError)
       if (ioError /= 0) exit
       if (n == 0) cycle

       select case (trim(words(1)))
       case ('cell')
          n_cells = n_cells + 1
       case ('surface')
          n_surfaces = n_surfaces + 1
       case ('material') 
          n_materials = n_materials + 1
       end select
    end do

    ! Check to make sure there are cells, surface, and materials
    ! defined for the problem
    if (n_cells == 0) then
       msg = "No cells specified!"
       close(UNIT=in)
       call error(msg)
    elseif (n_surfaces == 0) then
       msg = "No surfaces specified!"
       close(UNIT=in)
       call error(msg)
    elseif (n_materials == 0) then
       msg = "No materials specified!"
       close(UNIT=in)
       call error(msg)
    end if

    close(UNIT=in)

    ! Allocate arrays for cells, surfaces, etc.
    allocate(cells(n_cells))
    allocate(surfaces(n_surfaces))
    allocate(materials(n_materials))

    ! Set up dictionaries
    call dict_create(cell_dict)
    call dict_create(surface_dict)
    call dict_create(material_dict)

  end subroutine read_count

!=====================================================================
! READ_INPUT takes a second pass through the input file and parses all
! the data on each line, calling the appropriate subroutine for each
! type of data item.
!=====================================================================

  subroutine read_input(filename)

    character(*), intent(in) :: filename

    character(250) :: line, msg
    character(32) :: words(max_words)
    integer :: in = 7
    integer :: ioError
    integer :: n, i
    integer :: index_cell
    integer :: index_surface
    integer :: index_material
    integer :: index_source

    msg = "Second pass through input file..."
    call message(msg, 5)

    index_cell     = 0
    index_surface  = 0
    index_material = 0
    index_source   = 0

    open(FILE=filename, UNIT=in, STATUS='old', & 
         & ACTION='read', IOSTAT=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call error(msg)
    end if

    do
       call get_next_line(in, words, n, ioError)
       if ( ioError /= 0 ) exit

       ! skip comments
       if (n==0 .or. words(1)(1:1) == '#') cycle

       select case (trim(words(1)))
       case ('cell')
          ! Read data from cell entry
          index_cell = index_cell + 1
          call read_cell(index_cell, words, n)
          
       case ('surface')
          ! Read data
          index_surface = index_surface + 1
          call read_surface(index_surface, words, n)

       case ('material')
          index_material = index_material + 1
          call read_material(index_material, words, n)

       case ('source')
          call read_source(words, n)

       case ('criticality')
          call read_criticality(words, n)

       case ('xs_data')
          path_xsdata = words(2)

       case default
          ! do nothing
       end select

    end do

    close(UNIT=in)

    call adjust_indices()

  end subroutine read_input

!=====================================================================
! ADJUST_INDICES changes the boundary_list values for each cell and
! the material index assigned to each to the indices in the surfaces
! and material array rather than the unique IDs assigned to each
! surface and material
!=====================================================================

  subroutine adjust_indices()

    type(Cell), pointer :: c

    integer :: i, j
    integer :: index
    integer :: surf_num
    character(250) :: msg
    
    do i = 1, n_cells
       ! adjust boundary list
       c => cells(i)
       do j = 1, c%n_items
          surf_num = c%boundary_list(j)
          if (surf_num < OP_DIFFERENCE) then
             index = dict_get_key(surface_dict, abs(surf_num))
             if (index == DICT_NULL) then
                surf_num = abs(surf_num)
                msg = "Could not find surface " // trim(int_to_str(surf_num)) // &
                     & " specified on cell " // trim(int_to_str(c%uid))
                call error(msg)
             end if
             c%boundary_list(j) = sign(index, surf_num)
          end if
       end do

       ! adjust material indices
       index = dict_get_key(material_dict, c%material)
       if (index == DICT_NULL) then
          msg = "Could not find material " // trim(int_to_str(c%material)) // &
               & " specified on cell " // trim(int_to_str(c%uid))
          call error(msg)
       end if
       c%material = index
    end do

  end subroutine adjust_indices

!=====================================================================
! READ_CELL parses the data on a cell card.
!=====================================================================

  subroutine read_cell(index, words, n_words)

    integer,      intent(in) :: index
    character(*), intent(in) :: words(max_words)
    integer,      intent(in) :: n_words

    integer :: ioError
    integer :: i
    integer :: n_items
    character(250) :: msg
    character(32) :: word
    type(Cell), pointer :: this_cell => null()

    this_cell => cells(index)

    ! Read cell identifier
    read(words(2), FMT='(I8)', IOSTAT=ioError) this_cell%uid
    if (ioError > 0) then
       msg = "Invalid cell name: " // words(2)
       call error(msg)
    end if
    call dict_add_key(cell_dict, this_cell%uid, index)

    ! Read cell material
    read(words(3), FMT='(I8)', IOSTAT=ioError) this_cell%material

    ! Assign number of items
    n_items = n_words - 3
    this_cell%n_items = n_items

    ! Read list of surfaces
    allocate(this_cell%boundary_list(n_items))
    do i = 1, n_items
       word = words(i+3)
       if (word(1:1) == '(') then
          this_cell%boundary_list(i) = OP_LEFT_PAREN
       elseif (word(1:1) == ')') then
          this_cell%boundary_list(i) = OP_RIGHT_PAREN
       elseif (word(1:1) == ':') then
          this_cell%boundary_list(i) = OP_UNION
       elseif (word(1:1) == '#') then
          this_cell%boundary_list(i) = OP_DIFFERENCE
       else
          read(word, FMT='(I8)', IOSTAT=ioError) this_cell%boundary_list(i)
       end if
    end do

  end subroutine read_cell

!=====================================================================
! READ_SURFACE parses the data on a surface card.
!=====================================================================

  subroutine read_surface(index, words, n_words)

    integer, intent(in) :: index
    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

    integer :: ioError
    integer :: i
    integer :: coeffs_reqd
    character(250) :: msg
    character(32) :: word
    type(Surface), pointer :: this_surface => null()

    this_surface => surfaces(index)

    ! Read surface identifier
    read(words(2), FMT='(I8)', IOSTAT=ioError) this_surface%uid
    if (ioError > 0) then
       msg = "Invalid surface name: " // words(2)
       call error(msg)
    end if
    call dict_add_key(surface_dict, this_surface%uid, index)

    ! Read surface type
    word = words(3)
    call lower_case(word)
    select case (trim(word))
       case ('px')
          this_surface%type = SURF_PX
          coeffs_reqd  = 1
       case ('py')
          this_surface%type = SURF_PY
          coeffs_reqd  = 1
       case ('pz')
          this_surface%type = SURF_PZ
          coeffs_reqd  = 1
       case ('plane')
          this_surface%type = SURF_PLANE
          coeffs_reqd  = 4
       case ('cylx')
          this_surface%type = SURF_CYL_X
          coeffs_reqd  = 3
       case ('cyly')
          this_surface%type = SURF_CYL_Y
          coeffs_reqd  = 3
       case ('cylz')
          this_surface%type = SURF_CYL_Z
          coeffs_reqd  = 3
       case ('sph')
          this_surface%type = SURF_SPHERE
          coeffs_reqd  = 4
       case ('boxx')
          this_surface%type = SURF_BOX_X
          coeffs_reqd  = 4
       case ('boxy')
          this_surface%type = SURF_BOX_Y
          coeffs_reqd  = 4
       case ('boxz')
          this_surface%type = SURF_BOX_Z
          coeffs_reqd  = 4
       case ('box') 
          this_surface%type = SURF_BOX
          coeffs_reqd  = 6
       case ('gq')
          this_surface%type = SURF_GQ
          coeffs_reqd  = 10
       case default
          msg = "Invalid surface type: " // words(3)
          call error(msg)
    end select

    ! Make sure there are enough coefficients for surface type
    if (n_words-3 < coeffs_reqd) then
       msg = "Not enough coefficients for surface: " // words(2)
       call error(msg)
    end if
    
    ! Read list of surfaces
    allocate(this_surface%coeffs(n_words-3))
    do i = 1, n_words-3
       this_surface%coeffs(i) = str_to_real(words(i+3))
    end do

  end subroutine read_surface

!=====================================================================
! READ_SOURCE parses the data on a source card.
!=====================================================================

  subroutine read_source(words, n_words)

    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

    character(250) :: msg
    character(32) :: word
    integer :: ioError
    integer :: values_reqd
    integer :: i

    ! Read source type
    word = words(2)
    call lower_case(word)
    select case (trim(word))
    case ('box')
       external_source%type = SRC_BOX
       values_reqd = 6
    case default
       msg = "Invalid source type: " // words(2)
       call error(msg)
    end select

    ! Make sure there are enough values for this source type
    if (n_words-2 < values_reqd) then
       msg = "Not enough values for source of type: " // words(2)
       call error(msg)
    end if
    
    ! Read list of surfaces
    allocate(external_source%values(n_words-2))
    do i = 1, n_words-2
       external_source%values(i) = str_to_real(words(i+2))
    end do
    
  end subroutine read_source

!=====================================================================
! READ_MATERIAL parses a material card. Note that atom percents and
! densities are normalized in a separate routine
!=====================================================================

  subroutine read_material(index, words, n_words)

    integer, intent(in) :: index
    character(*), intent(in) :: words(n_words)
    integer, intent(in) :: n_words

    character(100) :: msg
    integer :: i
    integer :: ioError
    integer :: n_isotopes
    type(Material), pointer :: mat

    if (mod(n_words,2) == 0 .or. n_words < 5) then
       msg = "Invalid number of arguments for material: " // words(2)
       call error(msg)
    end if

    n_isotopes = (n_words-3)/2
    mat => materials(index)
    mat%n_isotopes = n_isotopes

    ! Read surface identifier
    read(words(2), FMT='(I8)', IOSTAT=ioError) mat%uid
    if (ioError > 0) then
       msg = "Invalid surface name: " // words(2)
       call error(msg)
    end if
    call dict_add_key(material_dict, mat%uid, index)

    ! Read atom density
    mat%atom_density = str_to_real(words(3))

    ! allocate isotope and density list
    allocate(mat%names(n_isotopes))
    allocate(mat%isotopes(n_isotopes))
    allocate(mat%table(n_isotopes))
    allocate(mat%atom_percent(n_isotopes))

    ! read names and percentage
    do i = 1, n_isotopes
       mat%names(i) = words(2*i+2)
       mat%atom_percent(i) = str_to_real(words(2*i+3))
    end do

  end subroutine read_material

!=====================================================================
! NORMALIZE_AO normalizes the atom or weight percentages for each
! material
!=====================================================================

  subroutine normalize_ao()

    type(xsData), pointer :: iso
    type(Material), pointer :: mat
    integer :: index
    integer :: i, j
    real(8) :: sum_percent
    real(8) :: awr              ! atomic weight ratio
    real(8) :: w                ! weight percent
    real(8) :: x                ! atom percent
    logical :: percent_in_atom
    logical :: density_in_atom
    character(10) :: key
    character(100) :: msg
    
    ! first find the index in the xsdata array for each isotope in
    ! each material
    do i = 1, n_materials
       mat => materials(i)

       ! Check to make sure either all atom percents or all weight
       ! percents are given
       if (.not. (all(mat%atom_percent > 0) .or. & 
            & all(mat%atom_percent < 0))) then
          msg = "Cannot mix atom and weight percents in material " // &
               & int_to_str(mat%uid)
          call error(msg)
       end if

       percent_in_atom = (mat%atom_percent(1) > 0)
       density_in_atom = (mat%atom_density > 0)

       sum_percent = 0
       do j = 1, mat%n_isotopes
          ! Set indices for isotopes
          key = mat%names(j)
          index = dict_get_key(xsdata_dict, key)
          mat%isotopes(j) = index

          ! determine atomic weight ratio
          awr = xsdatas(index)%awr

          ! if given weight percent, convert all values so that they
          ! are divided by awr. thus, when a sum is done over the
          ! values, it's actually sum(w/awr)
          if (.not. percent_in_atom) then
             mat%atom_percent(j) = -mat%atom_percent(j) / awr
          end if
       end do

       ! determine normalized atom percents. if given atom percents,
       ! this is straightforward. if given weight percents, the value
       ! is w/awr and is divided by sum(w/awr)
       sum_percent = sum(mat%atom_percent)
       mat%atom_percent = mat%atom_percent / sum_percent

       ! Change density in g/cm^3 to atom/b-cm. Since all values are
       ! now in atom percent, the sum needs to be re-evaluated as
       ! 1/sum(x*awr)
       if (.not. density_in_atom) then
          sum_percent = 0
          do j = 1, mat%n_isotopes
             index = mat%isotopes(j)
             awr = xsdatas(index)%awr
             x = mat%atom_percent(j)
             sum_percent = sum_percent + x*awr
          end do
          sum_percent = 1.0_8 / sum_percent
          mat%atom_density = -mat%atom_density * N_AVOGADRO & 
               & / MASS_NEUTRON * sum_percent
       end if
    end do

  end subroutine normalize_ao

!=====================================================================
! READ_XS_LIBRARY parses the data on a xs_library card. This card
! specifies what cross section library should be used by default
!=====================================================================

  subroutine read_xs_library( words, n_words )

    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

  end subroutine read_xs_library

!=====================================================================
! READ_CRITICALITY parses the data on a criticality card. This card
! specifies that the problem at hand is a criticality calculation and
! gives the number of cycles, inactive cycles, and particles per
! cycle.
!=====================================================================

  subroutine read_criticality( words, n_words )

    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

    character(250) :: msg
    integer :: ioError

    ! Set problem type to criticality
    problem_type = PROB_CRITICALITY

    ! Read number of cycles
    read(words(2), FMT='(I8)', IOSTAT=ioError) n_cycles
    if (ioError > 0) then
       msg = "Invalid number of cycles: " // words(2)
       call error(msg)
    end if

    ! Read number of inactive cycles
    read(words(3), FMT='(I8)', IOSTAT=ioError) n_inactive
    if (ioError > 0) then
       msg = "Invalid number of inactive cycles: " // words(2)
       call error(msg)
    end if

    ! Read number of particles
    read(words(4), FMT='(I8)', IOSTAT=ioError) n_particles
    if (ioError > 0) then
       msg = "Invalid number of particles: " // words(2)
       call error(msg)
    end if

  end subroutine read_criticality

!=====================================================================
! READ_LINE reads a line from a file open on a unit
!=====================================================================

  subroutine read_line(unit, line, ioError)

    integer,             intent(in)  :: unit
    character(max_line), intent(out) :: line
    integer,             intent(out) :: ioError

    read(UNIT=unit, FMT='(A250)', IOSTAT=ioError) line

  end subroutine read_line

!=====================================================================
! GET_NEXT_LINE reads the next line to the file connected on the
! specified unit including any continuation lines. If a line ends in
! an ampersand, the next line is read and its words are appended to
! the final array
!=====================================================================

  subroutine get_next_line(unit, words, n, ioError)

    integer,      intent(in)  :: unit
    character(*), intent(out) :: words(max_words)
    integer,      intent(out) :: n
    integer,      intent(out) :: ioError

    character(250) :: line 
    character(32)  :: local_words(max_words)
    integer        :: index

    index = 0
    do     
       read(UNIT=unit, FMT='(A250)', IOSTAT=ioError) line
       if (ioError /= 0) return
       call split_string_wl(line, local_words, n)
       if (n == 0) exit
       if (local_words(n) == '&') then
          words(index+1:index+n-1) = local_words(1:n-1)
          index = index + n-1
       else
          words(index+1:index+n) = local_words(1:n)
          index = index + n
          exit
       end if
    end do
    n = index

  end subroutine get_next_line

!=====================================================================
! SKIP_LINES skips 'n_lines' lines from a file open on a unit
!=====================================================================

  subroutine skip_lines(unit, n_lines, ioError)

    integer, intent(in) :: unit
    integer, intent(in) :: n_lines
    integer, intent(out) :: ioError

    integer :: i
    character(max_line) :: tmp

    do i = 1, n_lines
       read(UNIT=unit, FMT='(A250)', IOSTAT=ioError) tmp
    end do

  end subroutine skip_lines

!=====================================================================
! READ_DATA_INT reads integer data into an arrya from a file open
!=====================================================================

  subroutine read_data_int(in, array, n, lines, words_per_line)

    integer, intent(in) :: in
    integer, intent(out) :: array(n)
    integer, intent(in) :: n
    integer, intent(in) :: lines
    integer, intent(in) :: words_per_line

    integer :: i, j
    character(250) :: line
    character(32) :: words(max_words)
    integer :: ioError
    integer :: n_words
    integer :: index
    integer :: n_last

    index = 1
    do i = 1, lines
       if (i == lines) then
          read(in,*) array(index:n)
       else
          read(in,*) array(index:index+words_per_line-1)
          index = index + words_per_line
       end if
    end do

  end subroutine read_data_int

!=====================================================================
! READ_DATA_REAL reads real(8) data into an arrya from a file open
!=====================================================================

  subroutine read_data_real(in, array, n, lines, words_per_line)

    integer, intent(in)  :: in
    real(8), intent(out) :: array(n)
    integer, intent(in)  :: n
    integer, intent(in)  :: lines
    integer, intent(in)  :: words_per_line

    integer :: i, j
    character(250) :: line
    character(32) :: words(max_words)
    integer :: ioError
    integer :: n_words
    integer :: index

    index = 1
    do i = 1, lines
       if (i == lines) then
          read(in,*) array(index:n)
       else
          read(in,*) array(index:index+words_per_line-1)
          index = index + words_per_line
       end if
    end do

  end subroutine read_data_real

end module fileio
