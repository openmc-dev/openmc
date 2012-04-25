module input_xml

  use cmfd_input,      only: read_cmfd_xml 
  use constants
  use datatypes,       only: dict_add_key, dict_has_key, dict_get_key
  use error,           only: fatal_error, warning
  use geometry_header, only: Cell, Surface, Lattice
  use global
  use mesh_header,     only: StructuredMesh
  use output,          only: write_message
  use plot_header
  use random_lcg,      only: prn
  use string,          only: lower_case, to_str, str_to_int, str_to_real, &
                             split_string, starts_with, ends_with
  use tally_header,    only: TallyObject

  implicit none

  type(DictionaryII), pointer :: &  ! used to count how many cells each
       cells_in_univ_dict => null() ! universe contains

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    call read_settings_xml()
    call read_geometry_xml()
    call read_materials_xml()
    call read_tallies_xml()
    if (cmfd_on) call read_cmfd_xml()

  end subroutine read_input_xml

!===============================================================================
! READ_SETTINGS_XML reads data from a settings.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_settings_xml()

    use xml_data_settings_t

    integer :: n
    integer :: coeffs_reqd
    logical :: file_exists
    character(MAX_FILE_LEN) :: env_variable
    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: filename

    ! Display output message
    message = "Reading settings XML file..."
    call write_message(5)

    ! Check if settings.xml exists
    filename = trim(path_input) // "settings.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Settings XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if

    ! Initialize XML scalar variables
    cross_sections_ = ""
    verbosity_ = 0
    energy_grid_ = "union"
    seed_ = 0_8
    write_source_ = ""
    no_reduce_ = ""
    source_ % type = ""

    ! Parse settings.xml file
    call read_xml_file_settings_t(filename)

    ! Find cross_sections.xml file -- the first place to look is the
    ! settings.xml file. If no file is found there, then we check the
    ! CROSS_SECTIONS environment variable

    if (len_trim(cross_sections_) == 0) then
       ! No cross_sections.xml file specified in settings.xml, check environment
       ! variable
       call get_environment_variable("CROSS_SECTIONS", env_variable)
       if (len_trim(env_variable) == 0) then
          message = "No cross_sections.xml file was specified in settings.xml &
               &or in the CROSS_SECTIONS environment variable."
          call fatal_error()
       else
          path_cross_sections = trim(env_variable)
       end if
    else
       path_cross_sections = trim(cross_sections_)
    end if

    ! Criticality information
    if (criticality % batches > 0) then
       if (run_mode /= MODE_PLOTTING) run_mode = MODE_CRITICALITY

       ! Check number of particles
       if (len_trim(criticality % particles) == 0) then
          message = "Need to specify number of particles per cycles."
          call fatal_error()
       end if
       n_particles = str_to_int(criticality % particles)

       ! Copy cycle information
       n_batches     = criticality % batches
       n_inactive    = criticality % inactive
       n_active      = n_batches - n_inactive
       gen_per_batch = criticality % generations_per_batch

       ! Check number of active batches
       if (n_active <= 0) then
          message = "Number of active batches must be greater than 0."
          call fatal_error()
       end if
    else
       message = "Need to specify number of batches with <batches> tag."
       call fatal_error()
    end if

    ! Copy random number seed if specified
    if (seed_ > 0) seed = seed_

    ! Energy grid methods
    select case (energy_grid_)
    case ('nuclide')
       grid_method = GRID_NUCLIDE
    case ('union')
       grid_method = GRID_UNION
    case ('lethargy')
       message = "Lethargy mapped energy grid not yet supported."
       call fatal_error()
    case default
       message = "Unknown energy grid method: " // energy_grid_
       call fatal_error()
    end select

    ! Verbosity
    if (verbosity_ > 0) verbosity = verbosity_

    if (len(source_ % type) > 0) then
       ! Determine external source type
       type = source_ % type
       call lower_case(type)
       select case (type)
       case ('box')
          external_source % type = SRC_BOX
          coeffs_reqd = 6
       case ('point')
          external_source % type = SRC_POINT
          coeffs_reqd = 3
       case ('file')
          external_source % type = SRC_FILE
          external_source % path = trim(source_ % path)
       case default
          message = "Invalid source type: " // trim(source_ % type)
          call fatal_error()
       end select

       ! Coefficients for external surface
       if (type /= 'file') then
          n = size(source_ % coeffs)
          if (n < coeffs_reqd) then
             message = "Not enough coefficients specified for external source."
             call fatal_error()
          elseif (n > coeffs_reqd) then
             message = "Too many coefficients specified for external source."
             call fatal_error()
          else
             allocate(external_source % values(n))
             external_source % values = source_ % coeffs
          end if
       end if
    end if

    ! Survival biasing
    if (trim(survival_) == 'on') survival_biasing = .true.

    ! Probability tables
    if (ptables_ == 'off') urr_ptables_on = .false.

    ! Cutoffs
    if (size(cutoff_) > 0) then
       weight_cutoff = cutoff_(1) % weight
       weight_survive = cutoff_(1) % weight_avg
    end if

    ! Particle trace
    if (associated(trace_)) then
       trace_batch    = trace_(1)
       trace_gen      = trace_(2)
       trace_particle = trace_(3)
    end if

    ! Shannon Entropy mesh
    if (size(entropy_) > 0) then
       ! Check to make sure enough values were supplied
       if (size(entropy_(1) % lower_left) /= 3) then
          message = "Need to specify (x,y,z) coordinates of lower-left corner &
               &of Shannon entropy mesh."
       elseif (size(entropy_(1) % upper_right) /= 3) then
          message = "Need to specify (x,y,z) coordinates of upper-right corner &
               &of Shannon entropy mesh."
       end if

       ! Allocate mesh object and coordinates on mesh
       allocate(entropy_mesh)
       allocate(entropy_mesh % lower_left(3))
       allocate(entropy_mesh % upper_right(3))

       ! Copy values
       entropy_mesh % lower_left  = entropy_(1) % lower_left
       entropy_mesh % upper_right = entropy_(1) % upper_right

       ! Check on values provided
       if (.not. all(entropy_mesh % upper_right > entropy_mesh % lower_left)) then
          message = "Upper-right coordinate must be greater than lower-left &
               &coordinate for Shannon entropy mesh."
          call fatal_error()
       end if

       ! Check if dimensions were specified -- if not, they will be calculated
       ! automatically upon first entry into shannon_entropy

       if (associated(entropy_(1) % dimension)) then
          ! If so, make sure proper number of values were given
          if (size(entropy_(1) % dimension) /= 3) then
             message = "Dimension of entropy mesh must be given as three &
                  &integers."
             call fatal_error()
          end if

          ! Allocate dimensions
          entropy_mesh % n_dimension = 3
          allocate(entropy_mesh % dimension(3))

          ! Copy dimensions
          entropy_mesh % dimension = entropy_(1) % dimension
       end if

       ! Turn on Shannon entropy calculation
       entropy_on = .true.
    end if

    ! Uniform fission source weighting mesh
    if (size(uniform_fs_) > 0) then
       ! Check to make sure enough values were supplied
       if (size(uniform_fs_(1) % lower_left) /= 3) then
          message = "Need to specify (x,y,z) coordinates of lower-left corner &
               &of UFS mesh."
       elseif (size(uniform_fs_(1) % upper_right) /= 3) then
          message = "Need to specify (x,y,z) coordinates of upper-right corner &
               &of UFS mesh."
       elseif (size(uniform_fs_(1) % dimension) /= 3) then
          message = "Dimension of UFS mesh must be given as three &
               &integers."
          call fatal_error()
       end if

       ! Allocate mesh object and coordinates on mesh
       allocate(ufs_mesh)
       allocate(ufs_mesh % lower_left(3))
       allocate(ufs_mesh % upper_right(3))
       allocate(ufs_mesh % width(3))

       ! Allocate dimensions
       ufs_mesh % n_dimension = 3
       allocate(ufs_mesh % dimension(3))

       ! Copy dimensions
       ufs_mesh % dimension = uniform_fs_(1) % dimension

       ! Copy values
       ufs_mesh % lower_left  = uniform_fs_(1) % lower_left
       ufs_mesh % upper_right = uniform_fs_(1) % upper_right

       ! Check on values provided
       if (.not. all(ufs_mesh % upper_right > ufs_mesh % lower_left)) then
          message = "Upper-right coordinate must be greater than lower-left &
               &coordinate for UFS mesh."
          call fatal_error()
       end if

       ! Calculate width
       ufs_mesh % width = (ufs_mesh % upper_right - &
            ufs_mesh % lower_left) / ufs_mesh % dimension

       ! Calculate volume fraction of each cell
       ufs_mesh % volume_frac = ONE/real(product(ufs_mesh % dimension),8)

       ! Turn on uniform fission source weighting
       ufs = .true.

       ! Allocate source_frac
       allocate(source_frac(1, ufs_mesh % dimension(1), &
            ufs_mesh % dimension(2), ufs_mesh % dimension(3)))
    end if

    ! Check if the user has specified to write binary source file
    if (trim(write_source_) == 'on') write_source = .true.

    ! Check if the user has specified to not reduce tallies at the end of every
    ! batch
    if (trim(no_reduce_) == 'on') no_reduce = .true.

    ! Determine number of realizations
    if (no_reduce) then
       n_realizations = n_active * n_procs
    else
       n_realizations = n_active
    end if
       
    ! check for cmfd run
    if (run_cmfd_) cmfd_on = .true.

  end subroutine read_settings_xml

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    use xml_data_geometry_t

    integer :: i, j, k
    integer :: n
    integer :: n_x, n_y
    integer :: universe_num
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    logical :: file_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    type(Cell),    pointer :: c => null()
    type(Surface), pointer :: s => null()
    type(Lattice), pointer :: l => null()

    ! Display output message
    message = "Reading geometry XML file..."
    call write_message(5)

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Geometry XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if

    ! Parse geometry.xml file
    call read_xml_file_geometry_t(filename)

    ! Allocate cells array
    n_cells = size(cell_)
    allocate(cells(n_cells))

    n_universes = 0
    do i = 1, n_cells
       c => cells(i)
       
       ! Copy data into cells
       c % id       = cell_(i) % id
       c % universe = cell_(i) % universe
       c % material = cell_(i) % material
       c % fill     = cell_(i) % fill

       ! Check to make sure that either material or fill was specified
       if (c % material == 0 .and. c % fill == 0) then
          message = "Neither material nor fill was specified for cell " // & 
               trim(to_str(c % id))
          call fatal_error()
       end if

       ! Check to make sure that both material and fill haven't been
       ! specified simultaneously
       if (c % material /= 0 .and. c % fill /= 0) then
          message = "Cannot specify material and fill simultaneously"
          call fatal_error()
       end if

       ! Check to make sure that surfaces were specified
       if (.not. associated(cell_(i) % surfaces)) then
          message = "No surfaces specified for cell " // &
               trim(to_str(c % id))
          call fatal_error()
       end if

       ! Allocate array for surfaces and copy
       n = size(cell_(i) % surfaces)
       c % n_surfaces = n
       allocate(c % surfaces(n))
       c % surfaces = cell_(i) % surfaces

       ! Add cell to dictionary
       call dict_add_key(cell_dict, c % id, i)

       ! For cells, we also need to check if there's a new universe --
       ! also for every cell add 1 to the count of cells for the
       ! specified universe
       universe_num = cell_(i) % universe
       if (.not. dict_has_key(cells_in_univ_dict, universe_num)) then
          n_universes = n_universes + 1
          n_cells_in_univ = 1
          call dict_add_key(universe_dict, universe_num, n_universes)
       else
          n_cells_in_univ = 1 + dict_get_key(cells_in_univ_dict, universe_num)
       end if
       call dict_add_key(cells_in_univ_dict, universe_num, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! Allocate cells array
    n_surfaces = size(surface_)
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
       s => surfaces(i)
       
       ! Copy data into cells
       s % id = surface_(i) % id

       ! Copy and interpret surface type
       word = surface_(i) % type
       call lower_case(word)
       select case(trim(word))
       case ('x-plane')
          s % type = SURF_PX
          coeffs_reqd  = 1
       case ('y-plane')
          s % type = SURF_PY
          coeffs_reqd  = 1
       case ('z-plane')
          s % type = SURF_PZ
          coeffs_reqd  = 1
       case ('plane')
          s % type = SURF_PLANE
          coeffs_reqd  = 4
       case ('x-cylinder')
          s % type = SURF_CYL_X
          coeffs_reqd  = 3
       case ('y-cylinder')
          s % type = SURF_CYL_Y
          coeffs_reqd  = 3
       case ('z-cylinder')
          s % type = SURF_CYL_Z
          coeffs_reqd  = 3
       case ('sphere')
          s % type = SURF_SPHERE
          coeffs_reqd  = 4
       case ('box-x')
          s % type = SURF_BOX_X
          coeffs_reqd  = 4
       case ('box-y')
          s % type = SURF_BOX_Y
          coeffs_reqd  = 4
       case ('box-z')
          s % type = SURF_BOX_Z
          coeffs_reqd  = 4
       case ('box') 
          s % type = SURF_BOX
          coeffs_reqd  = 6
       case ('quadratic')
          s % type = SURF_GQ
          coeffs_reqd  = 10
       case default
          message = "Invalid surface type: " // trim(surface_(i) % type)
          call fatal_error()
       end select

       ! Check to make sure that the proper number of coefficients
       ! have been specified for the given type of surface. Then copy
       ! surface coordinates.

       n = size(surface_(i) % coeffs)
       if (n < coeffs_reqd) then
          message = "Not enough coefficients specified for surface: " // & 
               trim(to_str(s % id))
          call fatal_error()
       elseif (n > coeffs_reqd) then
          message = "Too many coefficients specified for surface: " // &
               trim(to_str(s % id))
          call fatal_error()
       else
          allocate(s % coeffs(n))
          s % coeffs = surface_(i) % coeffs
       end if

       ! Boundary conditions
       word = surface_(i) % boundary
       call lower_case(word)
       select case (trim(word))
       case ('transmission', 'transmit', '')
          s % bc = BC_TRANSMIT
       case ('vacuum')
          s % bc = BC_VACUUM
       case ('reflective', 'reflect', 'reflecting')
          s % bc = BC_REFLECT
       case ('periodic')
          s % bc = BC_PERIODIC
       case default
          message = "Unknown boundary condition '" // trim(word) // &
               "' specified on surface " // trim(to_str(s % id))
          call fatal_error()
       end select

       ! Add surface to dictionary
       call dict_add_key(surface_dict, s % id, i)

    end do

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    ! Allocate lattices array
    n_lattices = size(lattice_)
    allocate(lattices(n_lattices))

    do i = 1, n_lattices
       l => lattices(i)

       ! ID of lattice
       l % id = lattice_(i) % id

       ! Read lattice type
       word = lattice_(i) % type
       call lower_case(word)
       select case (trim(word))
       case ('rect', 'rectangle', 'rectangular')
          l % type = LATTICE_RECT
       case ('hex', 'hexagon', 'hexagonal')
          l % type = LATTICE_HEX
       case default
          message = "Invalid lattice type: " // trim(lattice_(i) % type)
          call fatal_error()
       end select

       ! Read number of lattice cells in each dimension
       n = size(lattice_(i) % dimension)
       if (n /= 2 .and. n /= 3) then
          message = "Lattice must be two or three dimensions."
          call fatal_error()
       end if
       n_x = lattice_(i) % dimension(1)
       n_y = lattice_(i) % dimension(2)
       l % n_x = n_x
       l % n_y = n_y

       ! Read lattice lower-left location
       if (size(lattice_(i) % dimension) /= size(lattice_(i) % lower_left)) then
          message = "Number of entries on <lower_left> must be the same as " // &
               "the number of entries on <dimension>."
          call fatal_error()
       end if
       l % x0 = lattice_(i) % lower_left(1)
       l % y0 = lattice_(i) % lower_left(2)

       ! Read lattice widths
       if (size(lattice_(i) % width) /= size(lattice_(i) % lower_left)) then
          message = "Number of entries on <width> must be the same as " // &
               "the number of entries on <lower_left>."
          call fatal_error()
       end if
       l % width_x = lattice_(i) % width(1)
       l % width_y = lattice_(i) % width(2)

       ! Read universes
       allocate(l % element(n_x, n_y))
       do k = 0, n_y - 1
          do j = 1, n_x
             l % element(j, n_y - k) = lattice_(i) % universes(j + k*n_x)
          end do
       end do

       ! Add lattice to dictionary
       call dict_add_key(lattice_dict, l % id, i)

    end do

  end subroutine read_geometry_xml

!===============================================================================
! READ_MATERIAL_XML reads data from a materials.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_materials_xml

    use xml_data_materials_t

    integer :: i           ! loop index for materials
    integer :: j           ! loop index for nuclides
    integer :: n           ! number of nuclides
    real(8) :: val         ! value entered for density
    logical :: file_exists ! does materials.xml exist?
    logical :: sum_density ! density is taken to be sum of nuclide densities
    character(3)            :: default_xs ! default xs identifier (e.g. 70c)
    character(12)           :: name       ! name of nuclide
    character(MAX_WORD_LEN) :: units      ! units on density
    character(MAX_LINE_LEN) :: filename   ! absolute path to materials.xml
    type(Material),    pointer :: m => null()
    type(nuclide_xml), pointer :: nuc => null()
    type(sab_xml),     pointer :: sab => null()

    ! Display output message
    message = "Reading materials XML file..."
    call write_message(5)

    ! Check is materials.xml exists
    filename = trim(path_input) // "materials.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Material XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if

    ! Initialize default cross section variable
    default_xs_ = ""

    ! Parse materials.xml file
    call read_xml_file_materials_t(filename)

    ! Copy default cross section if present
    default_xs = default_xs_

    ! Allocate cells array
    n_materials = size(material_)
    allocate(materials(n_materials))

    do i = 1, n_materials
       m => materials(i)

       ! Copy material id
       m % id = material_(i) % id

       ! Copy value and units
       val   = material_(i) % density % value
       units = material_(i) % density % units

       if (units == 'sum') then
          ! If the user gave the units as 'sum', then the total density of the
          ! material is taken to be the sum of the atom fractions listed on the
          ! nuclides

          sum_density = .true.

       else
          ! Check for erroneous density
          sum_density = .false.
          if (val <= ZERO) then
             message = "Need to specify a positive density on material " // &
                  trim(to_str(m % id)) // "."
             call fatal_error()
          end if

          ! Adjust material density based on specified units
          call lower_case(units)
          select case(trim(units))
          case ('g/cc', 'g/cm3')
             m % density = -val
          case ('kg/m3')
             m % density = -0.001 * val
          case ('atom/b-cm')
             m % density = val
          case ('atom/cm3', 'atom/cc')
             m % density = 1.0e-24 * val
          case default
             message = "Unkwown units '" // trim(material_(i) % density % units) &
                  // "' specified on material " // trim(to_str(m % id))
             call fatal_error()
          end select
       end if
       
       ! Check to ensure material has at least one nuclide
       if (.not. associated(material_(i) % nuclides)) then
          message = "No nuclides specified on material " // &
               trim(to_str(m % id))
          call fatal_error()
       end if

       ! allocate arrays in Material object
       n = size(material_(i) % nuclides)
       m % n_nuclides = n
       allocate(m % names(n))
       allocate(m % nuclide(n))
       allocate(m % xs_listing(n))
       allocate(m % atom_density(n))
       allocate(m % atom_percent(n))

       do j = 1, n
          ! Combine nuclide identifier and cross section and copy into names
          nuc => material_(i) % nuclides(j)

          ! Check for empty name on nuclide
          if (len_trim(nuc % name) == 0) then
             message = "No name specified on nuclide in material " // &
                  trim(to_str(m % id))
             call fatal_error()
          end if

          ! Check for cross section
          if (len_trim(nuc % xs) == 0) then
             if (default_xs == '') then
                message = "No cross section specified for nuclide in material " &
                     // trim(to_str(m % id))
                call fatal_error()
             else
                nuc % xs = default_xs
             end if
          end if

          ! copy full name
          name = trim(nuc % name) // "." // trim(nuc % xs)
          m % names(j) = name

          ! Check if no atom/weight percents were specified or if both atom and
          ! weight percents were specified
          if (nuc % ao == ZERO .and. nuc % wo == ZERO) then
             message = "No atom or weight percent specified for nuclide " // &
                  trim(name)
             call fatal_error()
          elseif (nuc % ao /= ZERO .and. nuc % wo /= ZERO) then
             message = "Cannot specify both atom and weight percents for a &
                  &nuclide: " // trim(name)
             call fatal_error()
          end if

          ! Copy atom/weight percents
          if (nuc % ao /= ZERO) then
             m % atom_percent(j) = nuc % ao
          else
             m % atom_percent(j) = -nuc % wo
          end if

          ! Read S(a,b) table information
          if (size(material_(i) % sab) == 1) then
             sab => material_(i) % sab(1)
             name = trim(sab % name) // "." // trim(sab % xs)
             m % sab_name = name
             m % has_sab_table = .true.
          elseif (size(material_(i) % sab) > 1) then
             message = "Cannot have multiple S(a,b) tables on a single material."
             call fatal_error()
          end if
       end do

       ! Determine density if it is a sum value
       if (sum_density) m % density = sum(m % atom_percent)

       ! Add material to dictionary
       call dict_add_key(material_dict, m % id, i)

    end do

  end subroutine read_materials_xml

!===============================================================================
! READ_TALLIES_XML reads data from a tallies.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_tallies_xml

    use xml_data_tallies_t

    integer :: i             ! loop over user-specified tallies
    integer :: j             ! loop over words
    integer :: i_analog      ! index in analog_tallies array
    integer :: i_tracklength ! index in tracklength_tallies array
    integer :: i_current     ! index in current_tallies array
    integer :: id            ! user-specified identifier
    integer :: i_mesh        ! index in meshes array
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filters     ! number of filters
    integer :: filters(N_FILTER_TYPES) ! temporary list of filters
    logical :: file_exists   ! does tallies.xml file exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: words(MAX_WORDS)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       ! Since a tallies.xml file is optional, no error is issued here
       n_tallies = 0
       return
    end if
    
    ! Display output message
    message = "Reading tallies XML file..."
    call write_message(5)

    ! Parse tallies.xml file
    call read_xml_file_tallies_t(filename)

    ! ==========================================================================
    ! DETERMINE SIZE OF ARRAYS AND ALLOCATE

    ! Check for meshes and allocate
    if (.not. associated(mesh_)) then
       n_meshes = 0
    else
       n_user_meshes = size(mesh_)
       if (cmfd_on) then
         n_meshes = n_user_meshes + n_cmfd_meshes
       else
         n_meshes = n_user_meshes
       end if
       allocate(meshes(n_meshes))
    end if

    ! Allocate tallies array
    if (.not. associated(tally_)) then
       n_tallies = 0
       message = "No tallies present in tallies.xml file!"
       call warning()
    else
       n_user_tallies = size(tally_)
       if (cmfd_on) then
         n_tallies = n_user_tallies + n_cmfd_tallies
       else
         n_tallies = n_user_tallies
       end if
       allocate(tallies(n_tallies))
    end if

    ! Check for <assume_separate> setting
    if (separate_ == 'yes') assume_separate = .true.

    ! ==========================================================================
    ! READ MESH DATA

    do i = 1, n_user_meshes
       m => meshes(i)

       ! copy mesh id
       m % id = mesh_(i) % id

       ! Read mesh type
       word = mesh_(i) % type
       call lower_case(word)
       select case (trim(word))
       case ('rect', 'rectangle', 'rectangular')
          m % type = LATTICE_RECT
       case ('hex', 'hexagon', 'hexagonal')
          m % type = LATTICE_HEX
       case default
          message = "Invalid mesh type: " // trim(mesh_(i) % type)
          call fatal_error()
       end select

       ! Determine number of dimensions for mesh
       n = size(mesh_(i) % dimension)
       if (n /= 2 .and. n /= 3) then
          message = "Mesh must be two or three dimensions."
          call fatal_error()
       end if
       m % n_dimension = n

       ! Allocate attribute arrays
       allocate(m % dimension(n))
       allocate(m % lower_left(n))
       allocate(m % width(n))
       allocate(m % upper_right(n))

       ! Read dimensions in each direction
       m % dimension = mesh_(i) % dimension

       ! Read mesh lower-left corner location
       if (m % n_dimension /= size(mesh_(i) % lower_left)) then
          message = "Number of entries on <lower_left> must be the same as &
               &the number of entries on <dimension>."
          call fatal_error()
       end if
       m % lower_left = mesh_(i) % lower_left

       ! Read mesh widths
       if (size(mesh_(i) % width) /= size(mesh_(i) % lower_left)) then
          message = "Number of entries on <width> must be the same as the &
               &number of entries on <lower_left>."
          call fatal_error()
       end if
       m % width = mesh_(i) % width

       ! Set upper right coordinate
       m % upper_right = m % lower_left + m % dimension * m % width

       ! Set volume fraction
       m % volume_frac = ONE/real(product(m % dimension),8)

       ! Add mesh to dictionary
       call dict_add_key(mesh_dict, m % id, i)
    end do

    ! ==========================================================================
    ! READ TALLY DATA

    do i = 1, n_user_tallies
       t => tallies(i)

       ! Allocate arrays for number of bins and stride in scores array
       allocate(t % n_filter_bins(N_FILTER_TYPES))
       allocate(t % stride(N_FILTER_TYPES))

       ! Initialize number of bins and stride
       t % n_filter_bins = 0
       t % stride = 0

       ! Initialize filters
       n_filters = 0
       filters = 0

       ! Set tally type to volume by default
       t % type = TALLY_VOLUME

       ! It's desirable to use a track-length esimator for tallies since
       ! generally more events will score to the tally, reducing the
       ! variance. However, for tallies that require information on
       ! post-collision parameters (e.g. tally with an energyout filter) the
       ! analog esimator must be used.

       if (trim(tally_(i) % estimator) == "analog") then
         t % estimator = ESTIMATOR_ANALOG
       else
         t % estimator = ESTIMATOR_TRACKLENGTH
       end if

       ! set tally reset property
       t % reset = tally_(i) % reset

       ! Copy material id
       t % id = tally_(i) % id

       ! Check to make sure that both cells and surfaces were not specified
       if (len_trim(tally_(i) % filters % cell) > 0 .and. &
            len_trim(tally_(i) % filters % surface) > 0) then
          message = "Cannot specify both cell and surface filters for tally " &
               // trim(to_str(t % id))
          call fatal_error()
       end if

       ! TODO: Parse logical expressions instead of just each word

       ! Read cell filter bins
       if (len_trim(tally_(i) % filters % cell) > 0) then
          call split_string(tally_(i) % filters % cell, words, n_words)
          allocate(t % cell_bins(n_words))
          do j = 1, n_words
             t % cell_bins(j) % scalar = int(str_to_int(words(j)),4)
          end do
          t % n_filter_bins(FILTER_CELL) = n_words

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_CELL
       end if

       ! Read surface filter bins
       if (len_trim(tally_(i) % filters % surface) > 0) then
          call split_string(tally_(i) % filters % surface, words, n_words)
          allocate(t % surface_bins(n_words))
          do j = 1, n_words
             t % surface_bins(j) % scalar = int(str_to_int(words(j)),4)
          end do
          t % n_filter_bins(FILTER_SURFACE) = n_words

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_SURFACE
       end if

       ! Read universe filter bins
       if (len_trim(tally_(i) % filters % universe) > 0) then
          call split_string(tally_(i) % filters % universe, words, n_words)
          allocate(t % universe_bins(n_words))
          do j = 1, n_words
             t % universe_bins(j) % scalar = int(str_to_int(words(j)),4)
          end do
          t % n_filter_bins(FILTER_UNIVERSE) = n_words

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_UNIVERSE
       end if

       ! Read material filter bins
       if (len_trim(tally_(i) % filters % material) > 0) then
          call split_string(tally_(i) % filters % material, words, n_words)
          allocate(t % material_bins(n_words))
          do j = 1, n_words
             t % material_bins(j) % scalar = int(str_to_int(words(j)),4)
          end do
          t % n_filter_bins(FILTER_MATERIAL) = n_words

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_MATERIAL
       end if

       ! Read mesh filter bins
       t % mesh = tally_(i) % filters % mesh
       if (t % mesh > 0) then
          ! Determine index in mesh array for this bin
          id = t % mesh
          if (dict_has_key(mesh_dict, id)) then
             i_mesh = dict_get_key(mesh_dict, id)
             m => meshes(i_mesh)
          else
             message = "Could not find mesh " // trim(to_str(id)) // &
                  " specified on tally " // trim(to_str(t % id))
             call fatal_error()
          end if

          t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) + product(m % dimension)

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_MESH
       end if

       ! Read birth region filter bins
       if (len_trim(tally_(i) % filters % cellborn) > 0) then
          call split_string(tally_(i) % filters % cellborn, words, n_words)
          allocate(t % cellborn_bins(n_words))
          do j = 1, n_words
             t % cellborn_bins(j) % scalar = int(str_to_int(words(j)),4)
          end do
          t % n_filter_bins(FILTER_CELLBORN) = n_words

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_CELLBORN
       end if

       ! Read incoming energy filter bins
       if (associated(tally_(i) % filters % energy)) then
          n = size(tally_(i) % filters % energy)
          allocate(t % energy_in(n))
          t % energy_in = tally_(i) % filters % energy
          t % n_filter_bins(FILTER_ENERGYIN) = n - 1

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_ENERGYIN
       end if

       ! Read outgoing energy filter bins
       if (associated(tally_(i) % filters % energyout)) then
          n = size(tally_(i) % filters % energyout)
          allocate(t % energy_out(n))
          t % energy_out = tally_(i) % filters % energyout
          t % n_filter_bins(FILTER_ENERGYOUT) = n - 1

          ! Set tally estimator to analog
          t % estimator = ESTIMATOR_ANALOG

          n_filters = n_filters + 1
          filters(n_filters) = FILTER_ENERGYOUT
       end if

       ! Allocate and set filters
       t % n_filters = n_filters
       allocate(t % filters(n_filters))
       t % filters = filters(1:n_filters)

       ! Read macro reactions
       if (len_trim(tally_(i) % scores) > 0) then
          call split_string(tally_(i) % scores, words, n_words)
          allocate(t % score_bins(n_words))
          do j = 1, n_words
             word = words(j)
             call lower_case(word)
             select case (trim(word))
             case ('flux')
                t % score_bins(j) % scalar = SCORE_FLUX
                if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
                   message = "Cannot tally flux with an outgoing energy filter."
                   call fatal_error()
                end if
             case ('total')
                t % score_bins(j) % scalar = SCORE_TOTAL
                if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
                   message = "Cannot tally total reaction rate with an &
                        &outgoing energy filter."
                   call fatal_error()
                end if
             case ('scatter')
                t % score_bins(j) % scalar = SCORE_SCATTER
             case ('nu-scatter')
                t % score_bins(j) % scalar = SCORE_NU_SCATTER

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('scatter-1')
                t % score_bins(j) % scalar = SCORE_SCATTER_1

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('scatter-2')
                t % score_bins(j) % scalar = SCORE_SCATTER_2

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('scatter-3')
                t % score_bins(j) % scalar = SCORE_SCATTER_3

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('diffusion')
                t % score_bins(j) % scalar = SCORE_DIFFUSION

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('n1n')
                t % score_bins(j) % scalar = SCORE_N_1N

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('n2n')
                t % score_bins(j) % scalar = SCORE_N_2N

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('n3n')
                t % score_bins(j) % scalar = SCORE_N_3N

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('n4n')
                t % score_bins(j) % scalar = SCORE_N_4N

                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
             case ('absorption')
                t % score_bins(j) % scalar = SCORE_ABSORPTION
                if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
                   message = "Cannot tally absorption rate with an outgoing &
                        &energy filter."
                   call fatal_error()
                end if
             case ('fission')
                t % score_bins(j) % scalar = SCORE_FISSION
                if (t % n_filter_bins(FILTER_ENERGYOUT) > 0) then
                   message = "Cannot tally fission rate with an outgoing &
                        &energy filter."
                   call fatal_error()
                end if
             case ('nu-fission')
                t % score_bins(j) % scalar = SCORE_NU_FISSION
             case ('current')
                t % score_bins(j) % scalar = SCORE_CURRENT
                t % type = TALLY_SURFACE_CURRENT

                ! Check to make sure that current is the only desired response
                ! for this tally
                if (n_words > 1) then
                   message = "Cannot tally other scoring functions in the same &
                        &tally as surface currents. Separate other scoring &
                        &functions into a distinct tally."
                   call fatal_error()
                end if

                ! Since the number of bins for the mesh filter was already set
                ! assuming it was a flux tally, we need to adjust the number of
                ! bins
                t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) &
                     - product(m % dimension)

                ! Get pointer to mesh
                id = t % mesh
                i_mesh = dict_get_key(mesh_dict, id)
                m => meshes(i_mesh)

                ! We need to increase the dimension by one since we also need
                ! currents coming into and out of the boundary mesh cells.
                if (size(m % dimension) == 2) then
                   t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) &
                        + product(m % dimension + 1) * 4
                elseif (size(m % dimension) == 3) then
                   t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) &
                        + product(m % dimension + 1) * 6
                end if

             case default
                message = "Unknown scoring function: " // trim(words(j))
                call fatal_error()
             end select
          end do
          t % n_score_bins = n_words
       end if

       ! Count number of tallies by type
       if (t % type == TALLY_VOLUME) then
          if (t % estimator == ESTIMATOR_ANALOG) then
             n_user_analog_tallies = n_user_analog_tallies + 1
          elseif (t % estimator == ESTIMATOR_TRACKLENGTH) then
             n_user_tracklength_tallies = n_user_tracklength_tallies + 1
          end if
       elseif (t % type == TALLY_SURFACE_CURRENT) then
          n_user_current_tallies = n_user_current_tallies + 1
       end if
    end do

    ! Determine number of types of tallies
    if (cmfd_on) then
      n_analog_tallies = n_user_analog_tallies + n_cmfd_analog_tallies
      n_tracklength_tallies = n_user_tracklength_tallies + n_cmfd_tracklength_tallies
      n_current_tallies = n_user_current_tallies + n_cmfd_current_tallies
    else
      n_analog_tallies = n_user_analog_tallies
      n_tracklength_tallies = n_user_tracklength_tallies
      n_current_tallies = n_user_current_tallies
    end if

    ! Allocate list of pointers for tallies by type
    allocate(analog_tallies(n_analog_tallies))
    allocate(tracklength_tallies(n_tracklength_tallies))
    allocate(current_tallies(n_current_tallies))

    ! Set indices for tally pointer lists to zero
    i_analog = 0
    i_tracklength = 0
    i_current = 0

    do i = 1, n_user_tallies
       t => tallies(i)

       ! Increment the appropriate index and set pointer
       if (t % type == TALLY_VOLUME) then
          if (t % estimator == ESTIMATOR_ANALOG) then
             i_analog = i_analog + 1
             analog_tallies(i_analog) = i
          elseif (t % estimator == ESTIMATOR_TRACKLENGTH) then
             i_tracklength = i_tracklength + 1
             tracklength_tallies(i_tracklength) = i
          end if
       elseif (t % type == TALLY_SURFACE_CURRENT) then
          i_current = i_current + 1
          current_tallies(i_current) = i
       end if
    end do

  end subroutine read_tallies_xml

!===============================================================================
! READ_PLOTS_XML reads data from a plots.xml file
!===============================================================================

  subroutine read_plots_xml

    use xml_data_plots_t

    integer i, j
    integer n_cols, col_id
    logical :: file_exists              ! does plots.xml file exist?
    character(MAX_LINE_LEN) :: filename ! absolute path to plots.xml
    type(Plot),         pointer :: pl => null()

    ! Check if plots.xml exists
    filename = trim(path_input) // "plots.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       message = "Plots XML file '" // trim(filename) // "' does not exist!"
       call fatal_error()
    end if
    
    ! Display output message
    message = "Reading plot XML file..."
    call write_message(5)

    ! Parse plots.xml file
    call read_xml_file_plots_t(filename)

    ! Allocate plots array
    n_plots = size(plot_)
    allocate(plots(n_plots))

    do i = 1, n_plots
      pl => plots(i)

      ! Copy data into plots
      pl % id       = plot_(i) % id

      ! Set output file path
      pl % path_plot = trim(path_input) // trim(to_str(pl % id)) // &
                       "_" // trim(plot_(i) % filename) // ".ppm"

      ! Copy plot pixel size
      if (size(plot_(i) % pixels) == 2) then
        pl % pixels = plot_(i) % pixels
      else
        message = "<pixels> must be length 2 in plot " // to_str(pl % id)
        call fatal_error()
      end if

      ! Copy plot background color
      if (associated(plot_(i) % background)) then
         if (size(plot_(i) % background) == 3) then
            pl % not_found % rgb = plot_(i) % background
         else
            message = "Bad background RGB " &
                 // "in plot " // trim(to_str(pl % id))
            call fatal_error()
         end if
      else
         pl % not_found % rgb = (/ 255, 255, 255 /)
      end if

      ! Copy plot type
      select case (plot_(i) % type)
        case ("slice")
          pl % type = PLOT_TYPE_SLICE
        !case ("points")
        !  pl % type = PLOT_TYPE_POINTS
        case default
          message = "Unsupported plot type '" // plot_(i) % type &
                    // "' in plot " // trim(to_str(pl % id))
          call fatal_error()
      end select

      ! Copy plot basis
      select case (plot_(i) % basis)
        case ("xy")
          pl % basis = PLOT_BASIS_XY
        case ("xz")
          pl % basis = PLOT_BASIS_XZ
        case ("yz")
          pl % basis = PLOT_BASIS_YZ
        case default
          message = "Unsupported plot basis '" // plot_(i) % basis &
                    // "' in plot " // trim(to_str(pl % id))
          call fatal_error()
      end select

      ! Copy plotting origin
      if (size(plot_(i) % origin) == 3) then
         pl % origin = plot_(i) % origin
      else
          message = "Origin must be length 3 " &
                    // "in plot " // trim(to_str(pl % id))
          call fatal_error()
      end if

      ! Copy plotting width
      if (size(plot_(i) % width) == 3) then
         pl % width = plot_(i) % width
      else if (size(plot_(i) % width) == 2) then
         pl % width(1) = plot_(i) % width(1)
         pl % width(2) = plot_(i) % width(2)
      else
          message = "Bad plot width " &
                    // "in plot " // trim(to_str(pl % id))
          call fatal_error()
      end if

      ! Copy plot color type and initialize all colors randomly
      select case (plot_(i) % color)
        case ("cell")

          pl % color_by = PLOT_COLOR_CELLS
          allocate(pl % colors(n_cells))
          do j = 1, n_cells
            pl % colors(j) % rgb(1) = prn()*255
            pl % colors(j) % rgb(2) = prn()*255
            pl % colors(j) % rgb(3) = prn()*255
          end do

        case ("mat", "material")

          pl % color_by = PLOT_COLOR_MATS
          allocate(pl % colors(n_materials))
          do j = 1, n_materials
            pl % colors(j) % rgb(1) = prn()*255
            pl % colors(j) % rgb(2) = prn()*255
            pl % colors(j) % rgb(3) = prn()*255
          end do

        case default
          message = "Unsupported plot color type '" // plot_(i) % color &
                    // "' in plot " // trim(to_str(pl % id))
          call fatal_error()
      end select

      ! Copy user specified colors
      if (associated(plot_(i) % col_spec_)) then
         n_cols = size(plot_(i) % col_spec_)
         do j = 1, n_cols
            if (size(plot_(i) % col_spec_(j) % rgb) /= 3) then
               message = "Bad RGB " &
                    // "in plot " // trim(to_str(pl % id))
               call fatal_error()          
            end if

            col_id = plot_(i) % col_spec_(j) % id

            if (pl % color_by == PLOT_COLOR_CELLS) then

               if (dict_has_key(cell_dict, col_id)) then
                  pl % colors(col_id) % rgb = plot_(i) % col_spec_(j) % rgb
               else
                  message = "Could not find cell " // trim(to_str(col_id)) // &
                       " specified in plot " // trim(to_str(pl % id))
                  call fatal_error()
               end if

            else if (pl % color_by == PLOT_COLOR_MATS) then

               if (dict_has_key(material_dict, col_id)) then
                  pl % colors(col_id) % rgb = plot_(i) % col_spec_(j) % rgb
               else
                  message = "Could not find material " // trim(to_str(col_id)) // &
                       " specified in plot " // trim(to_str(pl % id))
                  call fatal_error()
               end if

            end if
         end do
      end if

      ! Alter colors based on mask information
      if (associated(plot_(i) % mask_)) then
         if (size(plot_(i) % mask_) > 1) then
            message = "Mutliple masks" // &
                 " specified in plot " // trim(to_str(pl % id))
            call fatal_error()
         else if (.not. size(plot_(i) % mask_) == 0) then
            do j=1,size(pl % colors)
               if (.not. any(j .eq. plot_(i) % mask_(1) % components)) then
                  pl % colors(j) % rgb = plot_(i) % mask_(1) % background
               end if
            end do
         end if
      end if

    end do

  end subroutine read_plots_xml

!===============================================================================
! READ_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the ACE cross sections that may be used.
!===============================================================================

  subroutine read_cross_sections_xml()

    use xml_data_cross_sections_t

    integer :: i           ! loop index
    integer :: filetype    ! default file type
    integer :: recl        ! default record length
    integer :: entries     ! default number of entries
    logical :: file_exists ! does cross_sections.xml exist?
    character(MAX_WORD_LEN)  :: directory ! directory with cross sections
    type(XsListing), pointer :: listing => null()

    ! Check if cross_sections.xml exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
       ! Could not find cross_sections.xml file
       message = "Cross sections XML file '" // trim(path_cross_sections) // &
            "' does not exist!"
       call fatal_error()
    end if
    
    message = "Reading cross sections XML file..."
    call write_message(5)

    ! Initialize variables that may go unused
    directory_ = ""
    filetype_ = ""
    record_length_ = 0
    entries_ = 0

    ! Parse cross_sections.xml file
    call read_xml_file_cross_sections_t(path_cross_sections)

    if (len_trim(directory_) > 0) then
       ! Copy directory information if present
       directory = trim(directory_)
    else
       ! If no directory is listed in cross_sections.xml, by default select the
       ! directory in which the cross_sections.xml file resides
       i = index(path_cross_sections, "/", BACK=.true.)
       directory = path_cross_sections(1:i)
    end if

    ! determine whether binary/ascii
    if (filetype_ == 'ascii') then
       filetype = ASCII
    elseif (filetype_ == 'binary') then
       filetype = BINARY
    elseif (len_trim(filetype_) == 0) then
       filetype = ASCII
    else
       message = "Unknown filetype in cross_sections.xml: " // trim(filetype_)
       call fatal_error()
    end if

    ! copy default record length and entries for binary files
    recl = record_length_
    entries = entries_

    ! Allocate xs_listings array
    if (.not. associated(ace_tables_)) then
       message = "No ACE table listings present in cross_sections.xml file!"
       call fatal_error()
    else
       n_listings = size(ace_tables_)
       allocate(xs_listings(n_listings))
    end if
    

    do i = 1, n_listings
       listing => xs_listings(i)

       ! copy a number of attributes
       listing % name       = trim(ace_tables_(i) % name)
       listing % alias      = trim(ace_tables_(i) % alias)
       listing % zaid       = ace_tables_(i) % zaid
       listing % awr        = ace_tables_(i) % awr
       listing % kT         = ace_tables_(i) % temperature
       listing % location   = ace_tables_(i) % location

       ! determine type of cross section
       if (ends_with(listing % name, 'c')) then
          listing % type = ACE_NEUTRON
       elseif (ends_with(listing % name, 't')) then
          listing % type = ACE_THERMAL
       end if

       ! set filetype, record length, and number of entries
       listing % filetype = filetype
       listing % recl     = recl
       listing % entries  = entries

       ! determine metastable state
       if (ace_tables_(i) % metastable == 0) then
          listing % metastable = .false.
       else
          listing % metastable = .true.
       end if

       ! determine path of cross section table
       if (starts_with(ace_tables_(i) % path, '/')) then
          listing % path = ace_tables_(i) % path
       else
          if (ends_with(directory,'/')) then
             listing % path = trim(directory) // trim(ace_tables_(i) % path)
          else
             listing % path = trim(directory) // '/' // trim(ace_tables_(i) % path)
          end if
       end if

       ! create dictionary entry for both name and alias
       call dict_add_key(xs_listing_dict, listing % name, i)
       call dict_add_key(xs_listing_dict, listing % alias, i)
    end do

  end subroutine read_cross_sections_xml

end module input_xml
