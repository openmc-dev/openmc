module input_xml

  use cmfd_input,       only: configure_cmfd
  use constants
  use datatypes,        only: dict_add_key, dict_has_key, dict_get_key, &
                              dict_keys
  use datatypes_header, only: ListKeyValueCI
  use error,            only: fatal_error, warning
  use geometry_header,  only: Cell, Surface, Lattice
  use global
  use mesh_header,      only: StructuredMesh
  use output,           only: write_message
  use plot_header
  use random_lcg,       only: prn
  use string,           only: lower_case, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with
  use tally_header,     only: TallyObject, TallyFilter

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
    call read_cross_sections_xml()
    call read_geometry_xml()
    call read_materials_xml()
    call read_tallies_xml()
    if (cmfd_run) call configure_cmfd()

  end subroutine read_input_xml

!===============================================================================
! READ_SETTINGS_XML reads data from a settings.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_settings_xml()

    use xml_data_settings_t

    integer :: i ! loop index
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
    cross_sections_ = ''
    verbosity_ = 0
    energy_grid_ = 'union'
    seed_ = 0_8
    source_ % file = ''
    source_ % space % type = ''
    source_ % angle % type = ''
    source_ % energy % type = ''

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

    ! Make sure that either criticality or fixed source was specified
    if (eigenvalue_ % batches == 0 .and. fixed_source_ % batches == 0 &
         .and. criticality_ % batches == 0) then
      message = "Number of batches on <eigenvalue> or <fixed_source> &
           &tag was zero."
      call fatal_error()
    end if

    ! Check for old <criticality> tag
    if (criticality_ % batches > 0) then
      eigenvalue_ = criticality_
      message = "The <criticality> element has been deprecated and &
           &replaced by <eigenvalue>."
      call warning()
    end if

    ! Eigenvalue information
    if (eigenvalue_ % batches > 0) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_EIGENVALUE

      ! Check number of particles
      if (len_trim(eigenvalue_ % particles) == 0) then
        message = "Need to specify number of particles per cycles."
        call fatal_error()
      end if

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = str_to_int(eigenvalue_ % particles)

      ! Copy batch and generation information
      n_batches     = eigenvalue_ % batches
      n_inactive    = eigenvalue_ % inactive
      n_active      = n_batches - n_inactive
      gen_per_batch = eigenvalue_ % generations_per_batch

      ! Allocate array for batch keff and entropy
      allocate(k_batch(n_batches))
      allocate(entropy(n_batches))
      entropy = ZERO
    end if

    ! Fixed source calculation information
    if (fixed_source_ % batches > 0) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_FIXEDSOURCE

      ! Check number of particles
      if (len_trim(fixed_source_ % particles) == 0) then
        message = "Need to specify number of particles per cycles."
        call fatal_error()
      end if

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = str_to_int(fixed_source_ % particles)

      ! Copy batch information
      n_batches     = fixed_source_ % batches
      n_active      = fixed_source_ % batches
      n_inactive    = 0
      gen_per_batch = 1
    end if

    ! Check number of active batches, inactive batches, and particles
    if (n_active <= 0) then
      message = "Number of active batches must be greater than zero."
      call fatal_error()
    elseif (n_inactive < 0) then
      message = "Number of inactive batches must be non-negative."
      call fatal_error()
    elseif (n_particles <= 0) then
      message = "Number of particles must be greater than zero."
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

    ! ==========================================================================
    ! EXTERNAL SOURCE

    if (source_ % file /= '') then
      ! Copy path of source file
      path_source = source_ % file

      ! Check if source file exists
      inquire(FILE=path_source, EXIST=file_exists)
      if (.not. file_exists) then
        message = "Binary source file '" // trim(path_source) // &
             "' does not exist!"
        call fatal_error()
      end if

    else
      ! Spatial distribution for external source
      if (source_ % space % type /= '') then
        ! Read type of spatial distribution
        type = source_ % space % type
        call lower_case(type)
        select case (type)
        case ('box')
          external_source % type_space = SRC_SPACE_BOX
          coeffs_reqd = 6
        case ('point')
          external_source % type_space = SRC_SPACE_POINT
          coeffs_reqd = 3
        case default
          message = "Invalid spatial distribution for external source: " &
               // trim(source_ % space % type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (associated(source_ % space % parameters)) then
          n = size(source_ % space % parameters)
        else
          n = 0
        end if

        ! Read parameters for spatial distribution
        if (n < coeffs_reqd) then
          message = "Not enough parameters specified for spatial " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for spatial " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > 0) then
          allocate(external_source % params_space(n))
          external_source % params_space = source_ % space % parameters
        end if
      else
        message = "No spatial distribution specified for external source!"
        call fatal_error()
      end if

      ! Determine external source angular distribution
      if (source_ % angle % type /= '') then
        ! Read type of angular distribution
        type = source_ % angle % type
        call lower_case(type)
        select case (type)
        case ('isotropic')
          external_source % type_angle = SRC_ANGLE_ISOTROPIC
          coeffs_reqd = 0
        case ('monodirectional')
          external_source % type_angle = SRC_ANGLE_MONO
          coeffs_reqd = 3
        case ('tabular')
          external_source % type_angle = SRC_ANGLE_TABULAR
        case default
          message = "Invalid angular distribution for external source: " &
               // trim(source_ % angle % type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (associated(source_ % angle % parameters)) then
          n = size(source_ % angle % parameters)
        else
          n = 0
        end if

        ! Read parameters for angle distribution
        if (n < coeffs_reqd) then
          message = "Not enough parameters specified for angle " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for angle " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > 0) then
          allocate(external_source % params_angle(n))
          external_source % params_angle = source_ % angle % parameters
        end if
      else
        ! Set default angular distribution isotropic
        external_source % type_angle  = SRC_ANGLE_ISOTROPIC
      end if

      ! Determine external source energy distribution
      if (source_ % energy % type /= '') then
        ! Read type of energy distribution
        type = source_ % energy % type
        call lower_case(type)
        select case (type)
        case ('monoenergetic')
          external_source % type_energy = SRC_ENERGY_MONO
          coeffs_reqd = 1
        case ('maxwell')
          external_source % type_energy = SRC_ENERGY_MAXWELL
          coeffs_reqd = 1
        case ('watt')
          external_source % type_energy = SRC_ENERGY_WATT
          coeffs_reqd = 2
        case ('tabular')
          external_source % type_energy = SRC_ENERGY_TABULAR
        case default
          message = "Invalid energy distribution for external source: " &
               // trim(source_ % energy % type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (associated(source_ % energy % parameters)) then
          n = size(source_ % energy % parameters)
        else
          n = 0
        end if

        ! Read parameters for energy distribution
        if (n < coeffs_reqd) then
          message = "Not enough parameters specified for energy " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for energy " &
               // "distribution of external source."
          call fatal_error()
        elseif (n > 0) then
          allocate(external_source % params_energy(n))
          external_source % params_energy = source_ % energy % parameters
        end if
      else
        ! Set default energy distribution to Watt fission spectrum
        external_source % type_energy = SRC_ENERGY_WATT
        allocate(external_source % params_energy(2))
        external_source % params_energy = (/ 0.988_8, 2.249_8 /)
      end if
    end if

    ! Survival biasing
    call lower_case(survival_)
    if (survival_ == 'true' .or. survival_ == '1') survival_biasing = .true.

    ! Probability tables
    call lower_case(ptables_)
    if (ptables_ == 'false' .or. ptables_ == '0') urr_ptables_on = .false.

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

    ! Check if the user has specified to write state points
    if (size(state_point_) > 0) then
      ! Determine number of batches at which to store state points
      if (associated(state_point_(1) % batches)) then
        n_state_points = size(state_point_(1) % batches)
      else
        n_state_points = 0
      end if

      if (n_state_points > 0) then
        ! User gave specific batches to write state points
        allocate(statepoint_batch(n_state_points))
        statepoint_batch = state_point_(1) % batches

      elseif (state_point_(1) % interval /= 0) then
        ! User gave an interval for writing state points
        n_state_points = n_batches / state_point_(1) % interval
        allocate(statepoint_batch(n_state_points))
        statepoint_batch = (/ (state_point_(1) % interval * i, i = 1, &
             n_state_points) /)
      else
        ! If neither were specified, write state point at last batch
        n_state_points = 1
        allocate(statepoint_batch(n_state_points))
        statepoint_batch(1) = n_batches
      end if

      ! Check if the user has specified to write binary source file
      call lower_case(state_point_(1) % source_separate)
      if (state_point_(1) % source_separate == 'true' .or. &
           state_point_(1) % source_separate == '1') source_separate = .true.
    else
      ! If no <state_point> tag was present, by default write state point at
      ! last batch only
      n_state_points = 1
      allocate(statepoint_batch(n_state_points))
      statepoint_batch(1) = n_batches
    end if

    ! Check if the user has specified to not reduce tallies at the end of every
    ! batch
    call lower_case(no_reduce_)
    if (no_reduce_ == 'true' .or. no_reduce_ == '1') reduce_tallies = .false.

    ! Check if the user has specified to use confidence intervals for
    ! uncertainties rather than standard deviations
    call lower_case(confidence_intervals_)
    if (confidence_intervals_ == 'true' .or. &
         confidence_intervals_ == '1') confidence_intervals = .true.

    ! Check for output options
    if (associated(output_)) then
      do i = 1, size(output_)
        call lower_case(output_(i))
        select case (output_(i))
        case ('summary')
          output_summary = .true.
        case ('cross_sections')
          output_xs = .true.
        case ('tallies')
          output_tallies = .true.
        case ('none')
          output_summary = .false.
          output_xs = .false.
          output_tallies = .false.
        end select
      end do
    end if

    ! check for cmfd run
    call lower_case(run_cmfd_)
    if (run_cmfd_ == 'true' .or. run_cmfd_ == '1') then
      cmfd_run = .true.
#ifndef PETSC
      if (master) then
        message = 'CMFD is not available, compile OpenMC with PETSc'
        call fatal_error()
      end if
#endif
    end if

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
    real(8) :: phi, theta, psi
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

    ! Get number of <cell> tags
    n_cells = size(cell_)

    ! Check for no cells
    if (n_cells == 0) then
      message = "No cells found in geometry.xml!"
      call fatal_error()
    end if

    ! Allocate cells array
    allocate(cells(n_cells))

    n_universes = 0
    do i = 1, n_cells
      c => cells(i)

      ! Copy data into cells
      c % id       = cell_(i) % id
      c % universe = cell_(i) % universe
      c % fill     = cell_(i) % fill

      ! Read material
      word = cell_(i) % material
      call lower_case(word)
      select case(word)
      case ('void')
        c % material = MATERIAL_VOID

      case ('')
        ! This case is called if no material was specified
        c % material = 0

      case default
        c % material = int(str_to_int(word), 4)

        ! Check for error
        if (c % material == ERROR_INT) then
          message = "Invalid material specified on cell " // to_str(c % id)
          call fatal_error()
        end if
      end select

      ! Check to make sure that either material or fill was specified
      if (c % material == NONE .and. c % fill == NONE) then
        message = "Neither material nor fill was specified for cell " // & 
             trim(to_str(c % id))
        call fatal_error()
      end if

      ! Check to make sure that both material and fill haven't been
      ! specified simultaneously
      if (c % material /= NONE .and. c % fill /= NONE) then
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

      ! Rotation matrix
      if (associated(cell_(i) % rotation)) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a rotation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of rotation parameters
        n = size(cell_(i) % rotation)
        if (n /= 3) then
          message = "Incorrect number of rotation parameters on cell " // &
               to_str(c % id)
          call fatal_error()
        end if

        ! Copy rotation angles in x,y,z directions
        phi   = -cell_(i) % rotation(1) * PI/180.0
        theta = -cell_(i) % rotation(2) * PI/180.0
        psi   = -cell_(i) % rotation(3) * PI/180.0

        ! Calculate rotation matrix based on angles given
        allocate(c % rotation(3,3))
        c % rotation = reshape((/ &
             cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta), &
             -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), &
             cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), &
             sin(phi)*cos(theta), &
             sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi), &
             -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi), &
             cos(phi)*cos(theta) /), (/ 3,3 /))
      end if

      ! Translation vector
      if (associated(cell_(i) % translation)) then
        ! Translations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a translation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of translation parameters
        n = size(cell_(i) % translation)
        if (n /= 3) then
          message = "Incorrect number of translation parameters on cell " &
               // to_str(c % id)
          call fatal_error()
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        c % translation = cell_(i) % translation
      end if

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

    ! Get number of <surface> tags
    n_surfaces = size(surface_)

    ! Check for no surfaces
    if (n_surfaces == 0) then
      message = "No surfaces found in geometry.xml!"
      call fatal_error()
    end if

    ! Allocate cells array
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
      case ('x-cone')
        s % type = SURF_CONE_X
        coeffs_reqd  = 4
      case ('y-cone')
        s % type = SURF_CONE_Y
        coeffs_reqd  = 4
      case ('z-cone')
        s % type = SURF_CONE_Z
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

    integer :: i             ! loop index for materials
    integer :: j             ! loop index for nuclides
    integer :: n             ! number of nuclides
    integer :: index_list    ! index in xs_listings array
    integer :: index_nuclide ! index in nuclides
    integer :: index_sab     ! index in sab_tables
    real(8) :: val           ! value entered for density
    logical :: file_exists   ! does materials.xml exist?
    logical :: sum_density   ! density is taken to be sum of nuclide densities
    character(12) :: name       ! name of isotope, e.g. 92235.03c
    character(12) :: alias      ! alias of nuclide, e.g. U-235.03c
    character(MAX_WORD_LEN) :: units    ! units on density
    character(MAX_LINE_LEN) :: filename ! absolute path to materials.xml
    type(Material),    pointer :: mat => null()
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

    ! Initialize count for number of nuclides/S(a,b) tables
    index_nuclide = 0
    index_sab = 0

    do i = 1, n_materials
      mat => materials(i)

      ! Copy material id
      mat % id = material_(i) % id

      ! =======================================================================
      ! READ AND PARSE <density> TAG

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
               trim(to_str(mat % id)) // "."
          call fatal_error()
        end if

        ! Adjust material density based on specified units
        call lower_case(units)
        select case(trim(units))
        case ('g/cc', 'g/cm3')
          mat % density = -val
        case ('kg/m3')
          mat % density = -0.001 * val
        case ('atom/b-cm')
          mat % density = val
        case ('atom/cm3', 'atom/cc')
          mat % density = 1.0e-24 * val
        case default
          message = "Unkwown units '" // trim(material_(i) % density % units) &
               // "' specified on material " // trim(to_str(mat % id))
          call fatal_error()
        end select
      end if

      ! =======================================================================
      ! READ AND PARSE <nuclide> TAGS

      ! Check to ensure material has at least one nuclide
      if (.not. associated(material_(i) % nuclides)) then
        message = "No nuclides specified on material " // &
             trim(to_str(mat % id))
        call fatal_error()
      end if

      ! allocate arrays in Material object
      n = size(material_(i) % nuclides)
      mat % n_nuclides = n
      allocate(mat % names(n))
      allocate(mat % nuclide(n))
      allocate(mat % atom_density(n))

      do j = 1, mat % n_nuclides
        ! Combine nuclide identifier and cross section and copy into names
        nuc => material_(i) % nuclides(j)

        ! Check for empty name on nuclide
        if (len_trim(nuc % name) == 0) then
          message = "No name specified on nuclide in material " // &
               trim(to_str(mat % id))
          call fatal_error()
        end if

        ! Check for cross section
        if (len_trim(nuc % xs) == 0) then
          if (default_xs == '') then
            message = "No cross section specified for nuclide in material " &
                 // trim(to_str(mat % id))
            call fatal_error()
          else
            nuc % xs = default_xs
          end if
        end if

        ! copy full name
        name = trim(nuc % name) // "." // trim(nuc % xs)
        mat % names(j) = name

        ! Check that this nuclide is listed in the cross_sections.xml file
        if (.not. dict_has_key(xs_listing_dict, name)) then
          message = "Could not find nuclide " // trim(name) // &
               " in cross_sections.xml file!"
          call fatal_error()
        end if

        ! Check to make sure cross-section is continuous energy neutron table
        n = len_trim(name)
        if (name(n:n) /= 'c') then
          message = "Cross-section table " // trim(name) // & 
               " is not a continuous-energy neutron table."
          call fatal_error()
        end if

        ! Find xs_listing and set the name/alias according to the listing
        index_list = dict_get_key(xs_listing_dict, name)
        name       = xs_listings(index_list) % name
        alias      = xs_listings(index_list) % alias

        ! If this nuclide hasn't been encountered yet, we need to add its name
        ! and alias to the nuclide_dict
        if (.not. dict_has_key(nuclide_dict, name)) then
          index_nuclide    = index_nuclide + 1
          mat % nuclide(j) = index_nuclide

          call dict_add_key(nuclide_dict, name,  index_nuclide)
          call dict_add_key(nuclide_dict, alias, index_nuclide)
        else
          mat % nuclide(j) = dict_get_key(nuclide_dict, name)
        end if

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
          mat % atom_density(j) = nuc % ao
        else
          mat % atom_density(j) = -nuc % wo
        end if
      end do

      ! Check to make sure either all atom percents or all weight percents are
      ! given
      if (.not. (all(mat % atom_density > ZERO) .or. & 
           all(mat % atom_density < ZERO))) then
        message = "Cannot mix atom and weight percents in material " // &
             to_str(mat % id)
        call fatal_error()
      end if

      ! Determine density if it is a sum value
      if (sum_density) mat % density = sum(mat % atom_density)

      ! =======================================================================
      ! READ AND PARSE <sab> TAG FOR S(a,b) DATA

      if (size(material_(i) % sab) == 1) then
        ! Get pointer to S(a,b) table
        sab => material_(i) % sab(1)

        ! Determine name of S(a,b) table
        name = trim(sab % name) // "." // trim(sab % xs)
        mat % sab_name = name

        ! Check that this nuclide is listed in the cross_sections.xml file
        if (.not. dict_has_key(xs_listing_dict, name)) then
          message = "Could not find S(a,b) table " // trim(name) // &
               " in cross_sections.xml file!"
          call fatal_error()
        end if
        mat % has_sab_table = .true.

        ! Find index in xs_listing and set the name and alias according to the
        ! listing
        index_list = dict_get_key(xs_listing_dict, name)
        name       = xs_listings(index_list) % name
        alias      = xs_listings(index_list) % alias

        ! If this S(a,b) table hasn't been encountered yet, we need to add its
        ! name and alias to the sab_dict
        if (.not. dict_has_key(sab_dict, name)) then
          index_sab       = index_sab + 1
          mat % sab_table = index_sab

          call dict_add_key(sab_dict, name,  index_sab)
          call dict_add_key(sab_dict, alias, index_sab)
        else
          mat % sab_table = dict_get_key(sab_dict, name)
        end if

      elseif (size(material_(i) % sab) > 1) then
        message = "Cannot have multiple S(a,b) tables on a single material."
        call fatal_error()
      end if

      ! Add material to dictionary
      call dict_add_key(material_dict, mat % id, i)
    end do

    ! Set total number of nuclides and S(a,b) tables
    n_nuclides_total = index_nuclide
    n_sab_tables     = index_sab

  end subroutine read_materials_xml

!===============================================================================
! READ_TALLIES_XML reads data from a tallies.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_tallies_xml

    use xml_data_tallies_t

    integer :: i             ! loop over user-specified tallies
    integer :: j             ! loop over words
    integer :: k             ! another loop index
    integer :: i_analog      ! index in analog_tallies array
    integer :: i_tracklength ! index in tracklength_tallies array
    integer :: i_current     ! index in current_tallies array
    integer :: id            ! user-specified identifier
    integer :: i_mesh        ! index in meshes array
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filters     ! number of filters
    logical :: file_exists   ! does tallies.xml file exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    type(ListKeyValueCI), pointer :: key_list => null()
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter), allocatable :: filters(:) ! temporary filters

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! Since a tallies.xml file is optional, no error is issued here
      return
    end if

    ! Display output message
    message = "Reading tallies XML file..."
    call write_message(5)

    ! Parse tallies.xml file
    call read_xml_file_tallies_t(filename)

    ! ==========================================================================
    ! DETERMINE SIZE OF ARRAYS AND ALLOCATE

    ! Check for user meshes
    if (.not. associated(mesh_)) then
      n_user_meshes = 0
    else
      n_user_meshes = size(mesh_)
      if (cmfd_run) then
        n_meshes = n_user_meshes + n_cmfd_meshes
      else
        n_meshes = n_user_meshes
      end if
    end if

    ! Allocate mesh array
    if (n_meshes > 0) allocate(meshes(n_meshes))

    ! Check for user tallies
    if (.not. associated(tally_)) then
      n_user_tallies = 0
      message = "No tallies present in tallies.xml file!"
      call warning()
    else
      n_user_tallies = size(tally_)
      if (cmfd_run) then
        n_tallies = n_user_tallies + n_cmfd_tallies
      else
        n_tallies = n_user_tallies
      end if
    end if

    ! Allocate tally array
    if (n_tallies > 0) allocate(tallies(n_tallies))

    ! Check for <assume_separate> setting
    if (separate_ == 'yes') assume_separate = .true.

    ! ==========================================================================
    ! READ MESH DATA

    do i = 1, n_user_meshes
      m => meshes(i)

      ! copy mesh id
      m % id = mesh_(i) % id

      ! Read mesh type
      call lower_case(mesh_(i) % type)
      select case (mesh_(i) % type)
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

      ! Check that dimensions are all greater than zero
      if (any(mesh_(i) % dimension <= 0)) then
        message = "All entries on the <dimension> element for a tally mesh &
             &must be positive."
        call fatal_error()
      end if

      ! Read dimensions in each direction
      m % dimension = mesh_(i) % dimension

      ! Read mesh lower-left corner location
      if (m % n_dimension /= size(mesh_(i) % lower_left)) then
        message = "Number of entries on <lower_left> must be the same as &
             &the number of entries on <dimension>."
        call fatal_error()
      end if
      m % lower_left = mesh_(i) % lower_left

      ! Make sure either upper-right or width was specified
      if (associated(mesh_(i) % upper_right) .and. &
           associated(mesh_(i) % width)) then
        message = "Cannot specify both <upper_right> and <width> on a &
             &tally mesh."
        call fatal_error()
      end if

      ! Make sure either upper-right or width was specified
      if (.not. associated(mesh_(i) % upper_right) .and. &
           .not. associated(mesh_(i) % width)) then
        message = "Must specify either <upper_right> and <width> on a &
             &tally mesh."
        call fatal_error()
      end if

      if (associated(mesh_(i) % width)) then
        ! Check to ensure width has same dimensions
        if (size(mesh_(i) % width) /= size(mesh_(i) % lower_left)) then
          message = "Number of entries on <width> must be the same as the &
               &number of entries on <lower_left>."
          call fatal_error()
        end if

        ! Check for negative widths
        if (any(mesh_(i) % width < ZERO)) then
          message = "Cannot have a negative <width> on a tally mesh."
          call fatal_error()
        end if

        ! Set width and upper right coordinate
        m % width = mesh_(i) % width
        m % upper_right = m % lower_left + m % dimension * m % width

      elseif (associated(mesh_(i) % upper_right)) then
        ! Check to ensure width has same dimensions
        if (size(mesh_(i) % upper_right) /= size(mesh_(i) % lower_left)) then
          message = "Number of entries on <upper_right> must be the same as &
               &the number of entries on <lower_left>."
          call fatal_error()
        end if

        ! Check that upper-right is above lower-left
        if (any(mesh_(i) % upper_right < mesh_(i) % lower_left)) then
          message = "The <upper_right> coordinates must be greater than the &
               &<lower_left> coordinates on a tally mesh."
          call fatal_error()
        end if

        ! Set width and upper right coordinate
        m % upper_right = mesh_(i) % upper_right
        m % width = (m % upper_right - m % lower_left) / m % dimension
      end if

      ! Set volume fraction
      m % volume_frac = ONE/real(product(m % dimension),8)

      ! Add mesh to dictionary
      call dict_add_key(mesh_dict, m % id, i)
    end do

    ! ==========================================================================
    ! READ TALLY DATA

    READ_TALLIES: do i = 1, n_user_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Set tally type to volume by default
      t % type = TALLY_VOLUME

      ! It's desirable to use a track-length esimator for tallies since
      ! generally more events will score to the tally, reducing the
      ! variance. However, for tallies that require information on
      ! post-collision parameters (e.g. tally with an energyout filter) the
      ! analog esimator must be used.

      t % estimator = ESTIMATOR_TRACKLENGTH

      ! Copy material id
      t % id = tally_(i) % id

      ! Copy tally label
      t % label = tally_(i) % label

      ! =======================================================================
      ! READ DATA FOR FILTERS

      ! In older versions, tally filters were specified with a <filters>
      ! element followed by sub-elements <cell>, <mesh>, etc. This checks for
      ! the old format and if it is present, raises an error

      if (size(tally_(i) % filters) > 0) then
        message = "Tally filters should be specified with multiple <filter> &
             &elements. Did you forget to change your <filters> element?"
        call fatal_error()
      end if

      if (associated(tally_(i) % filter)) then
        ! Determine number of filters
        n_filters = size(tally_(i) % filter)

        ! Allocate filters array
        t % n_filters = n_filters
        allocate(t % filters(n_filters))

        READ_FILTERS: do j = 1, n_filters
          ! Convert filter type to lower case
          call lower_case(tally_(i) % filter(j) % type)

          ! Determine number of bins
          n_words = size(tally_(i) % filter(j) % bins)

          ! Determine type of filter
          select case (tally_(i) % filter(j) % type)
          case ('cell')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                   tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('cellborn')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELLBORN

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                   tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('material')
            ! Set type of filter
            t % filters(j) % type = FILTER_MATERIAL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                   tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('universe')
            ! Set type of filter
            t % filters(j) % type = FILTER_UNIVERSE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                   tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('surface')
            ! Set type of filter
            t % filters(j) % type = FILTER_SURFACE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            do k = 1, n_words
              t % filters(j) % int_bins(k) = int(str_to_int(&
                   tally_(i) % filter(j) % bins(k)),4)
            end do

          case ('mesh')
            ! Set type of filter
            t % filters(j) % type = FILTER_MESH

            ! Check to make sure multiple meshes weren't given
            if (n_words /= 1) then
              message = "Can only have one mesh filter specified."
              call fatal_error()
            end if

            ! Determine id of mesh
            id = int(str_to_int(tally_(i) % filter(j) % bins(1)),4)

            ! Get pointer to mesh
            if (dict_has_key(mesh_dict, id)) then
              i_mesh = dict_get_key(mesh_dict, id)
              m => meshes(i_mesh)
            else
              message = "Could not find mesh " // trim(to_str(id)) // &
                   " specified on tally " // trim(to_str(t % id))
              call fatal_error()
            end if

            ! Determine number of bins -- this is assuming that the tally is
            ! a volume tally and not a surface current tally. If it is a
            ! surface current tally, the number of bins will get reset later
            t % filters(j) % n_bins = product(m % dimension)

            ! Allocate and store index of mesh
            allocate(t % filters(j) % int_bins(1))
            t % filters(j) % int_bins(1) = i_mesh

          case ('energy')
            ! Set type of filter
            t % filters(j) % type = FILTER_ENERGYIN

            ! Set number of bins
            t % filters(j) % n_bins = n_words - 1

            ! Allocate and store bins
            allocate(t % filters(j) % real_bins(n_words))
            do k = 1, n_words
              t % filters(j) % real_bins(k) = str_to_real(&
                   tally_(i) % filter(j) % bins(k))
            end do

          case ('energyout')
            ! Set type of filter
            t % filters(j) % type = FILTER_ENERGYOUT

            ! Set number of bins
            t % filters(j) % n_bins = n_words - 1

            ! Allocate and store bins
            allocate(t % filters(j) % real_bins(n_words))
            do k = 1, n_words
              t % filters(j) % real_bins(k) = str_to_real(&
                   tally_(i) % filter(j) % bins(k))
            end do

          end select

          ! Set find_filter, e.g. if filter(3) has type FILTER_CELL, then
          ! find_filter(FILTER_CELL) would be set to 3.

          t % find_filter(t % filters(j) % type) = j

        end do READ_FILTERS

        ! Check that both cell and surface weren't specified
        if (t % find_filter(FILTER_CELL) > 0 .and. &
             t % find_filter(FILTER_SURFACE) > 0) then
          message = "Cannot specify both cell and surface filters for tally " &
               // trim(to_str(t % id))
          call fatal_error()
        end if

      else
        ! No filters were specified
        t % n_filters = 0
      end if

      ! =======================================================================
      ! READ DATA FOR NUCLIDES

      if (associated(tally_(i) % nuclides)) then
        if (tally_(i) % nuclides(1) == 'all') then
          ! Handle special case <nuclides>all</nuclides>
          allocate(t % nuclide_bins(n_nuclides_total + 1))

          ! Set bins to 1, 2, 3, ..., n_nuclides_total, -1
          t % nuclide_bins(1:n_nuclides_total) = &
               (/ (j, j=1, n_nuclides_total) /)
          t % nuclide_bins(n_nuclides_total + 1) = -1

          ! Set number of nuclide bins
          t % n_nuclide_bins = n_nuclides_total + 1

          ! Set flag so we can treat this case specially
          t % all_nuclides = .true.
        else
          ! Any other case, e.g. <nuclides>U-235 Pu-239</nuclides>
          n_words = size(tally_(i) % nuclides)
          allocate(t % nuclide_bins(n_words))
          do j = 1, n_words
            ! Check if total material was specified
            if (tally_(i) % nuclides(j) == 'total') then
              t % nuclide_bins(j) = -1
              cycle
            end if

            ! Check if xs specifier was given
            if (ends_with(tally_(i) % nuclides(j), 'c')) then
              word = tally_(i) % nuclides(j)
            else
              if (default_xs == '') then
                ! No default cross section specified, search through nuclides
                key_list => dict_keys(nuclide_dict)
                do while (associated(key_list))
                  if (starts_with(key_list % data % key, &
                       tally_(i) % nuclides(j))) then
                    word = key_list % data % key
                    exit
                  end if
                  
                  ! Advance to next
                  key_list => key_list % next
                end do

                ! Check if no nuclide was found
                if (.not. associated(key_list)) then
                  message = "Could not find the nuclide " // trim(&
                       tally_(i) % nuclides(j)) // " specified in tally " &
                       // trim(to_str(t % id)) // " in any material."
                  call fatal_error()
                end if
                deallocate(key_list)
              else
                ! Set nuclide to default xs
                word = trim(tally_(i) % nuclides(j)) // "." // default_xs
              end if
            end if

            ! Check to make sure nuclide specified is in problem
            if (.not. dict_has_key(nuclide_dict, word)) then
              message = "The nuclide " // trim(word) // " from tally " // &
                   trim(to_str(t % id)) // " is not present in any material."
              call fatal_error()
            end if

            ! Set bin to index in nuclides array
            t % nuclide_bins(j) = dict_get_key(nuclide_dict, word)
          end do

          ! Set number of nuclide bins
          t % n_nuclide_bins = n_words
        end if

      else
        ! No <nuclides> were specified -- create only one bin will be added
        ! for the total material.
        allocate(t % nuclide_bins(1))
        t % nuclide_bins(1) = -1
        t % n_nuclide_bins = 1
      end if

      ! =======================================================================
      ! READ DATA FOR SCORES

      if (associated(tally_(i) % scores)) then
        n_words = size(tally_(i) % scores)
        allocate(t % score_bins(n_words))
        do j = 1, n_words
          call lower_case(tally_(i) % scores(j))
          select case (tally_(i) % scores(j))
          case ('flux')
            ! Prohibit user from tallying flux for an individual nuclide
            if (.not. (t % n_nuclide_bins == 1 .and. &
                 t % nuclide_bins(1) == -1)) then
              message = "Cannot tally flux for an individual nuclide."
              call fatal_error()
            end if

            t % score_bins(j) = SCORE_FLUX
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              message = "Cannot tally flux with an outgoing energy filter."
              call fatal_error()
            end if
          case ('total')
            t % score_bins(j) = SCORE_TOTAL
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              message = "Cannot tally total reaction rate with an &
                   &outgoing energy filter."
              call fatal_error()
            end if
          case ('scatter')
            t % score_bins(j) = SCORE_SCATTER
          case ('nu-scatter')
            t % score_bins(j) = SCORE_NU_SCATTER

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('scatter-0') !does the same as 'scatter', for convenience
            t % score_bins(j) = SCORE_SCATTER
          case ('scatter-1')
            t % score_bins(j) = SCORE_SCATTER_1

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('scatter-2')
            t % score_bins(j) = SCORE_SCATTER_2

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('scatter-3')
            t % score_bins(j) = SCORE_SCATTER_3

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('scatter-4')
            t % score_bins(j) = SCORE_SCATTER_4

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('scatter-5')
            t % score_bins(j) = SCORE_SCATTER_5

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case('transport')
            t % score_bins(j) = SCORE_TRANSPORT

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('diffusion')
            t % score_bins(j) = SCORE_DIFFUSION

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('n1n')
            t % score_bins(j) = SCORE_N_1N

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('n2n')
            t % score_bins(j) = SCORE_N_2N

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('n3n')
            t % score_bins(j) = SCORE_N_3N

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('n4n')
            t % score_bins(j) = SCORE_N_4N

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('absorption')
            t % score_bins(j) = SCORE_ABSORPTION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              message = "Cannot tally absorption rate with an outgoing &
                   &energy filter."
              call fatal_error()
            end if
          case ('fission')
            t % score_bins(j) = SCORE_FISSION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              message = "Cannot tally fission rate with an outgoing &
                   &energy filter."
              call fatal_error()
            end if
          case ('nu-fission')
            t % score_bins(j) = SCORE_NU_FISSION
          case ('current')
            t % score_bins(j) = SCORE_CURRENT
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
            ! assuming it was a volume tally, we need to adjust the number
            ! of bins

            ! Get index of mesh filter
            k = t % find_filter(FILTER_MESH)

            ! Get pointer to mesh
            i_mesh = t % filters(k) % int_bins(1)
            m => meshes(i_mesh)

            ! We need to increase the dimension by one since we also need
            ! currents coming into and out of the boundary mesh cells.
            t % filters(k) % n_bins = product(m % dimension + 1)

            ! Copy filters to temporary array
            allocate(filters(t % n_filters + 1))
            filters(1:t % n_filters) = t % filters

            ! Move allocation back -- filters becomes deallocated during
            ! this call
            call move_alloc(FROM=filters, TO=t%filters)

            ! Add surface filter
            t % n_filters = t % n_filters + 1
            t % filters(t % n_filters) % type = FILTER_SURFACE
            t % filters(t % n_filters) % n_bins = 2 * m % n_dimension
            allocate(t % filters(t % n_filters) % int_bins(&
                 2 * m % n_dimension))
            if (m % n_dimension == 2) then
              t % filters(t % n_filters) % int_bins = (/ IN_RIGHT, &
                   OUT_RIGHT, IN_FRONT, OUT_FRONT /)
            elseif (m % n_dimension == 3) then
              t % filters(t % n_filters) % int_bins = (/ IN_RIGHT, &
                   OUT_RIGHT, IN_FRONT, OUT_FRONT, IN_TOP, OUT_TOP /)
            end if
            t % find_filter(FILTER_SURFACE) = t % n_filters

          case ('events')
            t % score_bins(j) = SCORE_EVENTS

          case default
            message = "Unknown scoring function: " // &
                 trim(tally_(i) % scores(j))
            call fatal_error()
          end select
        end do
        t % n_score_bins = n_words
      else
        message = "No <scores> specified on tally " // trim(to_str(t % id)) &
             // "."
        call fatal_error()
      end if

      ! =======================================================================
      ! SET TALLY ESTIMATOR

      ! Check if user specified estimator
      if (len_trim(tally_(i) % estimator) > 0) then
        select case(tally_(i) % estimator)
        case ('analog')
          t % estimator = ESTIMATOR_ANALOG

        case ('tracklength', 'track-length', 'pathlength', 'path-length')
          ! If the estimator was set to an analog estimator, this means the
          ! tally needs post-collision information
          if (t % estimator == ESTIMATOR_ANALOG) then
            message = "Cannot use track-length estimator for tally " &
                 // to_str(t % id)
            call fatal_error()
          end if

          ! Set estimator to track-length estimator
          t % estimator = ESTIMATOR_TRACKLENGTH

        case default
          message = "Invalid estimator '" // trim(tally_(i) % estimator) &
               // "' on tally " // to_str(t % id)
          call fatal_error()
        end select
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


    end do READ_TALLIES

    ! ==========================================================================
    ! LISTS FOR ANALOG, TRACKLENGTH, CURRENT TALLIES

    ! Determine number of types of tallies
    if (cmfd_run) then
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
          pl % colors(j) % rgb(1) = int(prn()*255)
          pl % colors(j) % rgb(2) = int(prn()*255)
          pl % colors(j) % rgb(3) = int(prn()*255)
        end do

      case ("mat", "material")

        pl % color_by = PLOT_COLOR_MATS
        allocate(pl % colors(n_materials))
        do j = 1, n_materials
          pl % colors(j) % rgb(1) = int(prn()*255)
          pl % colors(j) % rgb(2) = int(prn()*255)
          pl % colors(j) % rgb(3) = int(prn()*255)
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
