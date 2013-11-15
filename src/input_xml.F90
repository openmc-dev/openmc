module input_xml

  use cmfd_input,       only: configure_cmfd
  use constants
  use dict_header,      only: DictIntInt, ElemKeyValueCI
  use error,            only: fatal_error, warning
  use geometry_header,  only: Cell, Surface, Lattice
  use global
  use list_header,      only: ListChar, ListReal
  use mesh_header,      only: StructuredMesh
  use output,           only: write_message
  use plot_header
  use random_lcg,       only: prn
  use string,           only: lower_case, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with
  use tally_header,     only: TallyObject, TallyFilter
  use tally_initialize, only: add_tallies

  implicit none
  save

  type(DictIntInt) :: cells_in_univ_dict ! used to count how many cells each
                                         ! universe contains

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    call read_settings_xml()
    if ((run_mode /= MODE_PLOTTING)) call read_cross_sections_xml()
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
    integer :: n_tracks
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
    output_path_ = ''
    verbosity_ = 0
    threads_ = NONE
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

    if (len_trim(cross_sections_) == 0 .and. run_mode /= MODE_PLOTTING) then
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

    ! Set output directory if a path has been specified on the <output_path>
    ! element
    if (len_trim(output_path_) > 0) then
      path_output = output_path_
      if (.not. ends_with(path_output, "/")) &
           path_output = trim(path_output) // "/"
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
      allocate(k_generation(n_batches*gen_per_batch))
      allocate(entropy(n_batches*gen_per_batch))
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

    ! Number of OpenMP threads
    if (threads_ /= NONE) then
#ifdef OPENMP
      if (n_threads == NONE) then
        n_threads = threads_
        if (n_threads < 1) then
          message = "Invalid number of threads: " // to_str(n_threads)
          call fatal_error()
        end if
        call omp_set_num_threads(n_threads)
      end if
#else
      message = "Ignoring number of threads."
      call warning()
#endif
    end if

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
          message = "Not enough parameters specified for spatial &
               &distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for spatial &
               &distribution of external source."
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
          message = "Not enough parameters specified for angle &
               &distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for angle &
               &distribution of external source."
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
          message = "Not enough parameters specified for energy &
               &distribution of external source."
          call fatal_error()
        elseif (n > coeffs_reqd) then
          message = "Too many parameters specified for energy &
               &distribution of external source."
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

    ! Particle tracks
    if (associated(track_)) then
      n_tracks = size(track_)
      ! Make sure that there are three values per particle
      if (mod(n_tracks, 3) /= 0) then
        message = "Number of integers specified in 'track' is not divisible &
             &by 3.  Please provide 3 integers per particle to be tracked."
        call fatal_error()
      end if
      n_tracks = n_tracks/3
      allocate(track_identifiers(3,n_tracks))
      do i=1, n_tracks
        track_identifiers(1:3,i) = track_(3*(i-1)+1 : 3*(i-1)+3)
      end do
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
        message = "Dimension of UFS mesh must be given as three integers."
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
        do i = 1, n_state_points
          call statepoint_batch % add(state_point_(1) % batches(i))
        end do

      elseif (state_point_(1) % interval /= 0) then
        ! User gave an interval for writing state points
        n_state_points = n_batches / state_point_(1) % interval
        do i = 1, n_state_points
          call statepoint_batch % add(state_point_(1) % interval * i)
        end do
      else
        ! If neither were specified, write state point at last batch
        n_state_points = 1
        call statepoint_batch % add(n_batches)
      end if

      ! Check if the user has specified to write binary source file
      call lower_case(state_point_(1) % source_separate)
      call lower_case(state_point_(1) % source_write)
      if (state_point_(1) % source_separate == 'true' .or. &
           state_point_(1) % source_separate == '1') source_separate = .true.
      if (state_point_(1) % source_write == 'false' .or. &
           state_point_(1) % source_write == '0') source_write = .false.
    else
      ! If no <state_point> tag was present, by default write state point at
      ! last batch only
      n_state_points = 1
      call statepoint_batch % add(n_batches)
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

    ! Check for cmfd run
    call lower_case(run_cmfd_)
    if (run_cmfd_ == 'true' .or. run_cmfd_ == '1') then
      cmfd_run = .true.
#ifndef PETSC
      message = 'CMFD is not available, recompile OpenMC with PETSc'
      call fatal_error()
#endif
    end if

  end subroutine read_settings_xml

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    use xml_data_geometry_t

    integer :: i, j, k, m
    integer :: n
    integer :: n_x, n_y, n_z
    integer :: universe_num
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    integer :: mid
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    type(Cell),    pointer :: c => null()
    type(Surface), pointer :: s => null()
    type(Lattice), pointer :: lat => null()

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

    if (check_overlaps) then
      allocate(overlap_check_cnt(n_cells))
      overlap_check_cnt = 0
    end if

    n_universes = 0
    do i = 1, n_cells
      c => cells(i)

      ! Copy data into cells
      c % id       = cell_(i) % id
      c % universe = cell_(i) % universe
      c % fill     = cell_(i) % fill

      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has_key(c % id)) then
        message = "Two or more cells use the same unique ID: " // to_str(c % id)
        call fatal_error()
      end if

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
      call cell_dict % add_key(c % id, i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      universe_num = cell_(i) % universe
      if (.not. cells_in_univ_dict % has_key(universe_num)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % add_key(universe_num, n_universes)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get_key(universe_num)
      end if
      call cells_in_univ_dict % add_key(universe_num, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! This variable is used to check whether at least one boundary condition was
    ! applied to a surface
    boundary_exists = .false.

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

      ! Check to make sure 'id' hasn't been used
      if (surface_dict % has_key(s % id)) then
        message = "Two or more surfaces use the same unique ID: " // &
             to_str(s % id)
        call fatal_error()
      end if

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
        boundary_exists = .true.
      case ('reflective', 'reflect', 'reflecting')
        s % bc = BC_REFLECT
        boundary_exists = .true.
      case default
        message = "Unknown boundary condition '" // trim(word) // &
             "' specified on surface " // trim(to_str(s % id))
        call fatal_error()
      end select

      ! Add surface to dictionary
      call surface_dict % add_key(s % id, i)

    end do

    ! Check to make sure a boundary condition was applied to at least one
    ! surface
    if (.not. boundary_exists) then
      message = "No boundary conditions were applied to any surfaces!"
      call fatal_error()
    end if

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    ! Allocate lattices array
    n_lattices = size(lattice_)
    allocate(lattices(n_lattices))

    do i = 1, n_lattices
      lat => lattices(i)

      ! ID of lattice
      lat % id = lattice_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        message = "Two or more lattices use the same unique ID: " // &
             to_str(lat % id)
        call fatal_error()
      end if

      ! Read lattice type
      word = lattice_(i) % type
      call lower_case(word)
      select case (trim(word))
      case ('rect', 'rectangle', 'rectangular')
        lat % type = LATTICE_RECT
      case ('hex', 'hexagon', 'hexagonal')
        lat % type = LATTICE_HEX
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

      lat % n_dimension = n
      allocate(lat % dimension(n))
      lat % dimension = lattice_(i) % dimension

      ! Read lattice lower-left location
      if (size(lattice_(i) % dimension) /= size(lattice_(i) % lower_left)) then
        message = "Number of entries on <lower_left> must be the same as &
             &the number of entries on <dimension>."
        call fatal_error()
      end if

      allocate(lat % lower_left(n))
      lat % lower_left = lattice_(i) % lower_left

      ! Read lattice widths
      if (size(lattice_(i) % width) /= size(lattice_(i) % lower_left)) then
        message = "Number of entries on <width> must be the same as &
             &the number of entries on <lower_left>."
        call fatal_error()
      end if

      allocate(lat % width(n))
      lat % width = lattice_(i) % width

      ! Copy number of dimensions
      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if
      allocate(lat % universes(n_x, n_y, n_z))

      ! Check that number of universes matches size
      if (size(lattice_(i) % universes) /= n_x*n_y*n_z) then
        message = "Number of universes on <universes> does not match size of &
             &lattice " // trim(to_str(lat % id)) // "."
        call fatal_error()
      end if

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = lattice_(i) % &
                 universes(j + n_x*k + n_x*n_y*(m-1))
          end do
        end do
      end do

      ! Read material for area outside lattice
      mid = lattice_(i) % outside
      if (mid == 0 .or. mid == MATERIAL_VOID) then
        lat % outside = MATERIAL_VOID
      else
        lat % outside = mid
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, i)

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
    integer :: n_sab         ! number of sab tables for a material
    integer :: index_list    ! index in xs_listings array
    integer :: index_nuclide ! index in nuclides
    integer :: index_sab     ! index in sab_tables
    real(8) :: val           ! value entered for density
    logical :: file_exists   ! does materials.xml exist?
    logical :: sum_density   ! density is taken to be sum of nuclide densities
    character(12) :: name    ! name of isotope, e.g. 92235.03c
    character(12) :: alias   ! alias of nuclide, e.g. U-235.03c
    character(MAX_WORD_LEN) :: units    ! units on density
    character(MAX_LINE_LEN) :: filename ! absolute path to materials.xml
    type(ListChar) :: list_names   ! temporary list of nuclide names
    type(ListReal) :: list_density ! temporary list of nuclide densities
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

      ! Check to make sure 'id' hasn't been used
      if (material_dict % has_key(mat % id)) then
        message = "Two or more materials use the same unique ID: " // &
             to_str(mat % id)
        call fatal_error()
      end if

      if (run_mode == MODE_PLOTTING) then
        ! add to the dictionary and skip xs processing
        call material_dict % add_key(mat % id, i)
        cycle
      end if

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
      if (.not. associated(material_(i) % nuclides) .and. &
           .not. associated(material_(i) % elements)) then
        message = "No nuclides or natural elements specified on material " // &
             trim(to_str(mat % id))
        call fatal_error()
      end if

      ! Create list of nuclides based on those specified plus natural elements
      INDIVIDUAL_NUCLIDES: do j = 1, size(material_(i) % nuclides)
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

        ! store full name
        name = trim(nuc % name) // "." // trim(nuc % xs)

        ! save name and density to list
        call list_names % append(name)

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
          call list_density % append(nuc % ao)
        else
          call list_density % append(-nuc % wo)
        end if
      end do INDIVIDUAL_NUCLIDES

      ! =======================================================================
      ! READ AND PARSE <element> TAGS

      NATURAL_ELEMENTS: do j = 1, size(material_(i) % elements)
        nuc => material_(i) % elements(j)

        ! Check for empty name on natural element
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

        ! Expand element into naturally-occurring isotopes
        if (nuc % ao /= ZERO) then
          call expand_natural_element(nuc % name, nuc % xs, nuc % ao, &
               list_names, list_density)
        else
          message = "The ability to expand a natural element based on weight &
               &percentage is not yet supported."
          call fatal_error()
        end if
      end do NATURAL_ELEMENTS

      ! ========================================================================
      ! COPY NUCLIDES TO ARRAYS IN MATERIAL

      ! allocate arrays in Material object
      n = list_names % size()
      mat % n_nuclides = n
      allocate(mat % names(n))
      allocate(mat % nuclide(n))
      allocate(mat % atom_density(n))

      ALL_NUCLIDES: do j = 1, mat % n_nuclides
        ! Check that this nuclide is listed in the cross_sections.xml file
        name = list_names % get_item(j)
        if (.not. xs_listing_dict % has_key(name)) then
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
        index_list = xs_listing_dict % get_key(name)
        name       = xs_listings(index_list) % name
        alias      = xs_listings(index_list) % alias

        ! If this nuclide hasn't been encountered yet, we need to add its name
        ! and alias to the nuclide_dict
        if (.not. nuclide_dict % has_key(name)) then
          index_nuclide    = index_nuclide + 1
          mat % nuclide(j) = index_nuclide

          call nuclide_dict % add_key(name, index_nuclide)
          call nuclide_dict % add_key(alias, index_nuclide)
        else
          mat % nuclide(j) = nuclide_dict % get_key(name)
        end if

        ! Copy name and atom/weight percent
        mat % names(j) = name
        mat % atom_density(j) = list_density % get_item(j)
      end do ALL_NUCLIDES

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

      ! Clear lists
      call list_names % clear()
      call list_density % clear()

      ! =======================================================================
      ! READ AND PARSE <sab> TAG FOR S(a,b) DATA

      n_sab = size(material_(i) % sab)
      if (n_sab > 0) then
        ! Set number of S(a,b) tables
        mat % n_sab = n_sab

        ! Allocate names and indices for nuclides and tables
        allocate(mat % sab_names(n_sab))
        allocate(mat % i_sab_nuclides(n_sab))
        allocate(mat % i_sab_tables(n_sab))

        ! Initialize i_sab_nuclides
        mat % i_sab_nuclides = NONE

        do j = 1, n_sab
          ! Get pointer to S(a,b) table
          sab => material_(i) % sab(j)

          ! Determine name of S(a,b) table
          name = trim(sab % name) // "." // trim(sab % xs)
          mat % sab_names(j) = name

          ! Check that this nuclide is listed in the cross_sections.xml file
          if (.not. xs_listing_dict % has_key(name)) then
            message = "Could not find S(a,b) table " // trim(name) // &
                 " in cross_sections.xml file!"
            call fatal_error()
          end if

          ! Find index in xs_listing and set the name and alias according to the
          ! listing
          index_list = xs_listing_dict % get_key(name)
          name       = xs_listings(index_list) % name
          alias      = xs_listings(index_list) % alias

          ! If this S(a,b) table hasn't been encountered yet, we need to add its
          ! name and alias to the sab_dict
          if (.not. sab_dict % has_key(name)) then
            index_sab           = index_sab + 1
            mat % i_sab_tables(j) = index_sab

            call sab_dict % add_key(name, index_sab)
            call sab_dict % add_key(alias, index_sab)
          else
            mat % i_sab_tables(j) = sab_dict % get_key(name)
          end if
        end do
      end if

      ! Add material to dictionary
      call material_dict % add_key(mat % id, i)
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
    integer :: l             ! another loop index
    integer :: id            ! user-specified identifier
    integer :: i_mesh        ! index in meshes array
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filters     ! number of filters
    integer :: n_new         ! number of new scores to add based on Pn tally
    integer :: n_scores      ! number of tot scores after adjusting for Pn tally
    integer :: n_order       ! Scattering order requested
    integer :: n_order_pos   ! Position of Scattering order in score name string
    integer :: MT            ! user-specified MT for score
    logical :: file_exists   ! does tallies.xml file exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    type(ElemKeyValueCI), pointer :: pair_list => null()
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
    end if

    ! Allocate tally array
    if (n_user_tallies > 0) then
      call add_tallies("user", n_user_tallies)
    end if

    ! Check for <assume_separate> setting
    call lower_case(separate_)
    if (separate_ == 'true' .or. separate_ == '1') assume_separate = .true.

    ! ==========================================================================
    ! READ MESH DATA

    do i = 1, n_user_meshes
      m => meshes(i)

      ! copy mesh id
      m % id = mesh_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (mesh_dict % has_key(m % id)) then
        message = "Two or more meshes use the same unique ID: " // &
             to_str(m % id)
        call fatal_error()
      end if

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
      call mesh_dict % add_key(m % id, i)
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

      ! Check to make sure 'id' hasn't been used
      if (tally_dict % has_key(t % id)) then
        message = "Two or more tallies use the same unique ID: " // &
             to_str(t % id)
        call fatal_error()
      end if

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
            message = "Surface filter is not yet supported!"
            call fatal_error()

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
            if (mesh_dict % has_key(id)) then
              i_mesh = mesh_dict % get_key(id)
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

            ! Set to analog estimator
            t % estimator = ESTIMATOR_ANALOG

          case default
            ! Specified tally filter is invalid, raise error
            message = "Unknown filter type '" // trim(tally_(i) % &
                 filter(j) % type) // "' on tally " // &
                 trim(to_str(t % id)) // "."
            call fatal_error()

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
                pair_list => nuclide_dict % keys()
                do while (associated(pair_list))
                  if (starts_with(pair_list % key, &
                       tally_(i) % nuclides(j))) then
                    word = pair_list % key(1:150)
                    exit
                  end if
                  
                  ! Advance to next
                  pair_list => pair_list % next
                end do

                ! Check if no nuclide was found
                if (.not. associated(pair_list)) then
                  message = "Could not find the nuclide " // trim(&
                       tally_(i) % nuclides(j)) // " specified in tally " &
                       // trim(to_str(t % id)) // " in any material."
                  call fatal_error()
                end if
                deallocate(pair_list)
              else
                ! Set nuclide to default xs
                word = trim(tally_(i) % nuclides(j)) // "." // default_xs
              end if
            end if

            ! Check to make sure nuclide specified is in problem
            if (.not. nuclide_dict % has_key(word)) then
              message = "The nuclide " // trim(word) // " from tally " // &
                   trim(to_str(t % id)) // " is not present in any material."
              call fatal_error()
            end if

            ! Set bin to index in nuclides array
            t % nuclide_bins(j) = nuclide_dict % get_key(word)
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
        ! Loop through scores and determine if a scatter-p# input was used
        ! to allow for proper pre-allocating of t % score_bins
        ! This scheme allows multiple scatter-p# to be requested by the user
        ! if so desired
        n_words = size(tally_(i) % scores)
        n_new = 0
        do j = 1, n_words
          call lower_case(tally_(i) % scores(j))
          ! Find if scores(j) is of the form 'scatter-p'
          ! If so, get the number and do a select case on that.
          score_name = tally_(i) % scores(j)
          if (starts_with(score_name,'scatter-p')) then
            n_order_pos = scan(score_name,'0123456789')
            n_order = int(str_to_int( &
              score_name(n_order_pos:(len_trim(score_name)))),4)
            if (n_order > SCATT_ORDER_MAX) then
              ! Throw a warning. Set to the maximum number.
              ! The above scheme will essentially take the absolute value
              message = "Invalid scattering order of " // trim(to_str(n_order)) // &
                " requested. Setting to the maximum permissible value, " // &
                trim(to_str(SCATT_ORDER_MAX))
              call warning()
              n_order = SCATT_ORDER_MAX
              tally_(i) % scores(j) = SCATT_ORDER_MAX_PNSTR
            end if
            n_new = n_new + n_order
          end if
        end do
        n_scores = n_words + n_new
        
        ! Allocate accordingly
        allocate(t % score_bins(n_scores))
        allocate(t % scatt_order(n_scores))
        t % scatt_order = 0
        j = 0
        do l = 1, n_words
          j = j + 1
          ! Get the input string in scores(l) but if scatter-n or scatter-pn
          ! then strip off the n, and store it as an integer to be used later
          ! Peform the select case on this modified (number removed) string
          score_name = tally_(i) % scores(l)
          if (starts_with(score_name,'scatter-p')) then
            n_order_pos = scan(score_name,'0123456789')
            n_order = int(str_to_int( &
              score_name(n_order_pos:(len_trim(score_name)))),4)
            if (n_order > SCATT_ORDER_MAX) then
              ! Throw a warning. Set to the maximum number.
              ! The above scheme will essentially take the absolute value
              message = "Invalid scattering order of " // trim(to_str(n_order)) // &
                " requested. Setting to the maximum permissible value, " // &
                trim(to_str(SCATT_ORDER_MAX))
              call warning()
              n_order = SCATT_ORDER_MAX
            end if
            score_name = "scatter-pn"
          else if (starts_with(score_name,'scatter-')) then
            n_order_pos = scan(score_name,'0123456789')
            n_order = int(str_to_int( &
              score_name(n_order_pos:(len_trim(score_name)))),4)
            if (n_order > SCATT_ORDER_MAX) then
              ! Throw a warning. Set to the maximum number.
              ! The above scheme will essentially take the absolute value
              message = "Invalid scattering order of " // trim(to_str(n_order)) // &
                " requested. Setting to the maximum permissible value, " // &
                trim(to_str(SCATT_ORDER_MAX))
              call warning()
              n_order = SCATT_ORDER_MAX
            end if
            score_name = "scatter-n"
          end if
          
          select case (score_name)
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
          case ('scatter-n')
            if (n_order == 0) then
              t % score_bins(j) = SCORE_SCATTER
            else
              t % score_bins(j) = SCORE_SCATTER_N
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if
            t % scatt_order(j) = n_order
            
          case ('scatter-pn')
            t % estimator = ESTIMATOR_ANALOG
            ! Setup P0:Pn
            t % score_bins(j : j + n_order) = SCORE_SCATTER_PN
            t % scatt_order(j : j + n_order) = n_order
            j = j + n_order
            
          case('transport')
            t % score_bins(j) = SCORE_TRANSPORT

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('diffusion')
            message = "Diffusion score no longer supported for tallies, & 
                      &please remove"
            call fatal_error()
          case ('n1n')
            t % score_bins(j) = SCORE_N_1N

            ! Set tally estimator to analog
            t % estimator = ESTIMATOR_ANALOG
          case ('n2n')
            t % score_bins(j) = N_2N

          case ('n3n')
            t % score_bins(j) = N_3N

          case ('n4n')
            t % score_bins(j) = N_4N

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
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if
          case ('kappa-fission')
            t % score_bins(j) = SCORE_KAPPA_FISSION
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
            ! Assume that user has specified an MT number
            MT = int(str_to_int(score_name))

            if (MT /= ERROR_INT) then
              ! Specified score was an integer
              if (MT > 1) then
                t % score_bins(j) = MT
              else
                message = "Invalid MT on <scores>: " // &
                     trim(tally_(i) % scores(j))
                call fatal_error()
              end if

            else
              ! Specified score was not an integer
              message = "Unknown scoring function: " // &
                   trim(tally_(i) % scores(j))
              call fatal_error()
            end if

          end select
        end do
        t % n_score_bins = n_scores
        t % n_user_score_bins = n_words
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

      ! Add tally to dictionary
      call tally_dict % add_key(t % id, i)

    end do READ_TALLIES

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
    type(ObjectPlot), pointer :: pl => null()

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

    READ_PLOTS: do i = 1, n_plots
      pl => plots(i)

      ! Copy data into plots
      pl % id = plot_(i) % id

      ! Check to make sure 'id' hasn't been used
      if (plot_dict % has_key(pl % id)) then
        message = "Two or more plots use the same unique ID: " // &
             to_str(pl % id)
        call fatal_error()
      end if

      ! Copy plot type
      select case (plot_(i) % type)
      case ("slice")
        pl % type = PLOT_TYPE_SLICE
      case ("voxel")
        pl % type = PLOT_TYPE_VOXEL
      case default
        message = "Unsupported plot type '" // trim(plot_(i) % type) &
             // "' in plot " // trim(to_str(pl % id))
        call fatal_error()
      end select

      ! Set output file path
      select case (pl % type)
      case (PLOT_TYPE_SLICE)
        pl % path_plot = trim(path_input) // trim(to_str(pl % id)) // &
             "_" // trim(plot_(i) % filename) // ".ppm"
      case (PLOT_TYPE_VOXEL)
        pl % path_plot = trim(path_input) // trim(to_str(pl % id)) // &
             "_" // trim(plot_(i) % filename) // ".voxel"
      end select
      
      ! Copy plot pixel size
      if (pl % type == PLOT_TYPE_SLICE) then
        if (size(plot_(i) % pixels) == 2) then
          pl % pixels(1) = plot_(i) % pixels(1)
          pl % pixels(2) = plot_(i) % pixels(2)
        else
          message = "<pixels> must be length 2 in slice plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (size(plot_(i) % pixels) == 3) then
          pl % pixels = plot_(i) % pixels
        else
          message = "<pixels> must be length 3 in voxel plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      end if

      ! Copy plot background color
      if (associated(plot_(i) % background)) then
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Background color ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
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
      
      ! Copy plot basis
      if (pl % type == PLOT_TYPE_SLICE) then
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
      end if
      
      ! Copy plotting origin
      if (size(plot_(i) % origin) == 3) then
        pl % origin = plot_(i) % origin
      else
        message = "Origin must be length 3 " &
             // "in plot " // trim(to_str(pl % id))
        call fatal_error()
      end if

      ! Copy plotting width
      if (pl % type == PLOT_TYPE_SLICE) then
        if (size(plot_(i) % width) == 2) then
          pl % width(1) = plot_(i) % width(1)
          pl % width(2) = plot_(i) % width(2)
        else
          message = "<width> must be length 2 in slice plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (size(plot_(i) % width) == 3) then
          pl % width = plot_(i) % width
        else
          message = "<width> must be length 3 in voxel plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
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
      
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Color specifications ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
      
        n_cols = size(plot_(i) % col_spec_)
        do j = 1, n_cols
          if (size(plot_(i) % col_spec_(j) % rgb) /= 3) then
            message = "Bad RGB " &
                 // "in plot " // trim(to_str(pl % id))
            call fatal_error()          
          end if

          col_id = plot_(i) % col_spec_(j) % id

          if (pl % color_by == PLOT_COLOR_CELLS) then

            if (cell_dict % has_key(col_id)) then
              col_id = cell_dict % get_key(col_id)
              pl % colors(col_id) % rgb = plot_(i) % col_spec_(j) % rgb
            else
              message = "Could not find cell " // trim(to_str(col_id)) // &
                   " specified in plot " // trim(to_str(pl % id))
              call fatal_error()
            end if

          else if (pl % color_by == PLOT_COLOR_MATS) then

            if (material_dict % has_key(col_id)) then
              col_id = material_dict % get_key(col_id)
              pl % colors(col_id) % rgb = plot_(i) % col_spec_(j) % rgb
            else
              message = "Could not find material " // trim(to_str(col_id)) // &
                   " specified in plot " // trim(to_str(pl % id))
              call fatal_error()
            end if

          end if
        end do
      end if

      ! Deal with masks
      if (associated(plot_(i) % mask_)) then
      
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Mask ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
      
        select case(size(plot_(i) % mask_))
          case default
            message = "Mutliple masks" // &
                 " specified in plot " // trim(to_str(pl % id))
            call fatal_error()
          case (0)
          case (1)
          
            ! First we need to change the user-specified identifiers to indices
            ! in the cell and material arrays
            do j=1,size(plot_(i) % mask_(1) % components)
              col_id = plot_(i) % mask_(1) % components(j)
            
              if (pl % color_by == PLOT_COLOR_CELLS) then
              
                if (cell_dict % has_key(col_id)) then
                  plot_(i) % mask_(1) % components(j) = cell_dict % get_key(col_id)
                else
                  message = "Could not find cell " // trim(to_str(col_id)) // &
                       " specified in the mask in plot " // trim(to_str(pl % id))
                  call fatal_error()
                end if
              
              else if (pl % color_by == PLOT_COLOR_MATS) then
              
                if (material_dict % has_key(col_id)) then
                  plot_(i) % mask_(1) % components(j) = material_dict % get_key(col_id)
                else
                  message = "Could not find material " // trim(to_str(col_id)) // &
                       " specified in the mask in plot " // trim(to_str(pl % id))
                  call fatal_error()
                end if
                
              end if  
            end do
          
            ! Alter colors based on mask information
            do j=1,size(pl % colors)
              if (.not. any(j .eq. plot_(i) % mask_(1) % components)) then
                pl % colors(j) % rgb = plot_(i) % mask_(1) % background
              end if
            end do
            
        end select
        
      end if

      ! Add plot to dictionary
      call plot_dict % add_key(pl % id, i)

    end do READ_PLOTS

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
       call xs_listing_dict % add_key(listing % name, i)
       if (listing % alias /= '') then
         call xs_listing_dict % add_key(listing % alias, i)
       end if
    end do

  end subroutine read_cross_sections_xml

!===============================================================================
! EXPAND_NATURAL_ELEMENT converts natural elements specified using an <element>
! tag within a material into individual isotopes based on IUPAC Isotopic
! Compositions of the Elements 2009 (doi:10.1351/PAC-REP-10-06-02). In some
! cases, modifications have been made to work with ENDF/B-VII.1 where
! evaluations of particular isotopes don't exist.
!===============================================================================

  subroutine expand_natural_element(name, xs, density, list_names, &
       list_density)

    character(*),   intent(in)    :: name
    character(*),   intent(in)    :: xs
    real(8),        intent(in)    :: density
    type(ListChar), intent(inout) :: list_names
    type(ListReal), intent(inout) :: list_density

    character(2) :: element_name

    element_name = name(1:2)
    call lower_case(element_name)

    select case (element_name)
    case ('h')
      call list_names % append('1001.' // xs)
      call list_density % append(density * 0.999885_8)
      call list_names % append('1002.' // xs)
      call list_density % append(density * 0.000115_8)

    case ('he')
      call list_names % append('2003.' // xs)
      call list_density % append(density * 0.00000134_8)
      call list_names % append('2004.' // xs)
      call list_density % append(density * 0.99999866_8)

    case ('li')
      call list_names % append('3006.' // xs)
      call list_density % append(density * 0.0759_8)
      call list_names % append('3007.' // xs)
      call list_density % append(density * 0.9241_8)

    case ('be')
      call list_names % append('4009.' // xs)
      call list_density % append(density)

    case ('b')
      call list_names % append('5010.' // xs)
      call list_density % append(density * 0.199_8)
      call list_names % append('5011.' // xs)
      call list_density % append(density * 0.801_8)

    case ('c')
      ! The evaluation of Carbon in ENDF/B-VII.1 and JEFF 3.1.2 is a natural
      ! element, i.e. it's not possible to split into C-12 and C-13.
      call list_names % append('6000.' // xs)
      call list_density % append(density)

    case ('n')
      call list_names % append('7014.' // xs)
      call list_density % append(density * 0.99636_8)
      call list_names % append('7015.' // xs)
      call list_density % append(density * 0.00364_8)

    case ('o')
      ! O-18 does not exist in ENDF/B-VII.1 or JEFF 3.1.2 so its 0.205% has been
      ! added to O-16. The isotopic abundance for O-16 is ordinarily 99.757%.
      call list_names % append('8016.' // xs)
      call list_density % append(density * 0.99962_8)
      call list_names % append('8017.' // xs)
      call list_density % append(density * 0.00038_8)

    case ('f')
      call list_names % append('9019.' // xs)
      call list_density % append(density)

    case ('ne')
      call list_names % append('10020.' // xs)
      call list_density % append(density * 0.9048_8)
      call list_names % append('10021.' // xs)
      call list_density % append(density * 0.0027_8)
      call list_names % append('10022.' // xs)
      call list_density % append(density * 0.0925_8)

    case ('na')
      call list_names % append('11023.' // xs)
      call list_density % append(density)

    case ('mg')
      call list_names % append('12024.' // xs)
      call list_density % append(density * 0.7899_8)
      call list_names % append('12025.' // xs)
      call list_density % append(density * 0.1000_8)
      call list_names % append('12026.' // xs)
      call list_density % append(density * 0.1101_8)

    case ('al')
      call list_names % append('13027.' // xs)
      call list_density % append(density)

    case ('si')
      call list_names % append('14028.' // xs)
      call list_density % append(density * 0.92223_8)
      call list_names % append('14029.' // xs)
      call list_density % append(density * 0.04685_8)
      call list_names % append('14030.' // xs)
      call list_density % append(density * 0.03092_8)

    case ('p')
      call list_names % append('15031.' // xs)
      call list_density % append(density)

    case ('s')
      call list_names % append('16032.' // xs)
      call list_density % append(density * 0.9499_8)
      call list_names % append('16033.' // xs)
      call list_density % append(density * 0.0075_8)
      call list_names % append('16034.' // xs)
      call list_density % append(density * 0.0425_8)
      call list_names % append('16036.' // xs)
      call list_density % append(density * 0.0001_8)

    case ('cl')
      call list_names % append('17035.' // xs)
      call list_density % append(density * 0.7576_8)
      call list_names % append('17037.' // xs)
      call list_density % append(density * 0.2424_8)

    case ('ar')
      call list_names % append('18036.' // xs)
      call list_density % append(density * 0.003336_8)
      call list_names % append('18038.' // xs)
      call list_density % append(density * 0.000629_8)
      call list_names % append('18040.' // xs)
      call list_density % append(density * 0.996035_8)

    case ('k')
      call list_names % append('19039.' // xs)
      call list_density % append(density * 0.932581_8)
      call list_names % append('19040.' // xs)
      call list_density % append(density * 0.000117_8)
      call list_names % append('19041.' // xs)
      call list_density % append(density * 0.067302_8)

    case ('ca')
      call list_names % append('20040.' // xs)
      call list_density % append(density * 0.96941_8)
      call list_names % append('20042.' // xs)
      call list_density % append(density * 0.00647_8)
      call list_names % append('20043.' // xs)
      call list_density % append(density * 0.00135_8)
      call list_names % append('20044.' // xs)
      call list_density % append(density * 0.02086_8)
      call list_names % append('20046.' // xs)
      call list_density % append(density * 0.00004_8)
      call list_names % append('20048.' // xs)
      call list_density % append(density * 0.00187_8)

    case ('sc')
      call list_names % append('21045.' // xs)
      call list_density % append(density)

    case ('ti')
      call list_names % append('22046.' // xs)
      call list_density % append(density * 0.0825_8)
      call list_names % append('22047.' // xs)
      call list_density % append(density * 0.0744_8)
      call list_names % append('22048.' // xs)
      call list_density % append(density * 0.7372_8)
      call list_names % append('22049.' // xs)
      call list_density % append(density * 0.0541_8)
      call list_names % append('22050.' // xs)
      call list_density % append(density * 0.0518_8)

    case ('v')
      ! The evaluation of Vanadium in ENDF/B-VII.1 and JEFF 3.1.2 is a natural
      ! element. The IUPAC isotopic composition specifies the following
      ! breakdown which is not used:
      !   V-50 =  0.250%
      !   V-51 = 99.750%
      call list_names % append('23000.' // xs)
      call list_density % append(density)

    case ('cr')
      call list_names % append('24050.' // xs)
      call list_density % append(density * 0.04345_8)
      call list_names % append('24052.' // xs)
      call list_density % append(density * 0.83789_8)
      call list_names % append('24053.' // xs)
      call list_density % append(density * 0.09501_8)
      call list_names % append('24054.' // xs)
      call list_density % append(density * 0.02365_8)

    case ('mn')
      call list_names % append('25055.' // xs)
      call list_density % append(density)

    case ('fe')
      call list_names % append('26054.' // xs)
      call list_density % append(density * 0.05845_8)
      call list_names % append('26056.' // xs)
      call list_density % append(density * 0.91754_8)
      call list_names % append('26057.' // xs)
      call list_density % append(density * 0.02119_8)
      call list_names % append('26058.' // xs)
      call list_density % append(density * 0.00282_8)

    case ('co')
      call list_names % append('27059.' // xs)
      call list_density % append(density)

    case ('ni')
      call list_names % append('28058.' // xs)
      call list_density % append(density * 0.68077_8)
      call list_names % append('28060.' // xs)
      call list_density % append(density * 0.26223_8)
      call list_names % append('28061.' // xs)
      call list_density % append(density * 0.011399_8)
      call list_names % append('28062.' // xs)
      call list_density % append(density * 0.036346_8)
      call list_names % append('28064.' // xs)
      call list_density % append(density * 0.009255_8)

    case ('cu')
      call list_names % append('29063.' // xs)
      call list_density % append(density * 0.6915_8)
      call list_names % append('29065.' // xs)
      call list_density % append(density * 0.3085_8)

    case ('zn')
      ! The evaluation of Zinc in ENDF/B-VII.1 is a natural element. The IUPAC
      ! isotopic composition specifies the following breakdown which is not used
      ! here:
      !   Zn-64 = 48.63%
      !   Zn-66 = 27.90%
      !   Zn-67 =  4.10%
      !   Zn-68 = 18.75%
      !   Zn-70 =  0.62%
      call list_names % append('30000.' // xs)
      call list_density % append(density)

    case ('ga')
      ! JEFF 3.1.2 does not have evaluations for Ga-69 and Ga-71, only for
      ! natural Gallium, so this may cause problems.
      call list_names % append('31069.' // xs)
      call list_density % append(density * 0.60108_8)
      call list_names % append('31071.' // xs)
      call list_density % append(density * 0.39892_8)

    case ('ge')
      call list_names % append('32070.' // xs)
      call list_density % append(density * 0.2057_8)
      call list_names % append('32072.' // xs)
      call list_density % append(density * 0.2745_8)
      call list_names % append('32073.' // xs)
      call list_density % append(density * 0.0775_8)
      call list_names % append('32074.' // xs)
      call list_density % append(density * 0.3650_8)
      call list_names % append('32076.' // xs)
      call list_density % append(density * 0.0773_8)

    case ('as')
      call list_names % append('33075.' // xs)
      call list_density % append(density)

    case ('se')
      call list_names % append('34074.' // xs)
      call list_density % append(density * 0.0089_8)
      call list_names % append('34076.' // xs)
      call list_density % append(density * 0.0937_8)
      call list_names % append('34077.' // xs)
      call list_density % append(density * 0.0763_8)
      call list_names % append('34078.' // xs)
      call list_density % append(density * 0.2377_8)
      call list_names % append('34080.' // xs)
      call list_density % append(density * 0.4961_8)
      call list_names % append('34082.' // xs)
      call list_density % append(density * 0.0873_8)

    case ('br')
      call list_names % append('35079.' // xs)
      call list_density % append(density * 0.5069_8)
      call list_names % append('35081.' // xs)
      call list_density % append(density * 0.4931_8)

    case ('kr')
      call list_names % append('36078.' // xs)
      call list_density % append(density * 0.00355_8)
      call list_names % append('36080.' // xs)
      call list_density % append(density * 0.02286_8)
      call list_names % append('36082.' // xs)
      call list_density % append(density * 0.11593_8)
      call list_names % append('36083.' // xs)
      call list_density % append(density * 0.11500_8)
      call list_names % append('36084.' // xs)
      call list_density % append(density * 0.56987_8)
      call list_names % append('36086.' // xs)
      call list_density % append(density * 0.17279_8)

    case ('rb')
      call list_names % append('37085.' // xs)
      call list_density % append(density * 0.7217_8)
      call list_names % append('37087.' // xs)
      call list_density % append(density * 0.2783_8)

    case ('sr')
      call list_names % append('38084.' // xs)
      call list_density % append(density * 0.0056_8)
      call list_names % append('38086.' // xs)
      call list_density % append(density * 0.0986_8)
      call list_names % append('38087.' // xs)
      call list_density % append(density * 0.0700_8)
      call list_names % append('38088.' // xs)
      call list_density % append(density * 0.8258_8)

    case ('y')
      call list_names % append('39089.' // xs)
      call list_density % append(density)

    case ('zr')
      call list_names % append('40090.' // xs)
      call list_density % append(density * 0.5145_8)
      call list_names % append('40091.' // xs)
      call list_density % append(density * 0.1122_8)
      call list_names % append('40092.' // xs)
      call list_density % append(density * 0.1715_8)
      call list_names % append('40094.' // xs)
      call list_density % append(density * 0.1738_8)
      call list_names % append('40096.' // xs)
      call list_density % append(density * 0.0280_8)

    case ('nb')
      call list_names % append('41093.' // xs)
      call list_density % append(density)

    case ('mo')
      call list_names % append('42092.' // xs)
      call list_density % append(density * 0.1453_8)
      call list_names % append('42094.' // xs)
      call list_density % append(density * 0.0915_8)
      call list_names % append('42095.' // xs)
      call list_density % append(density * 0.1584_8)
      call list_names % append('42096.' // xs)
      call list_density % append(density * 0.1667_8)
      call list_names % append('42097.' // xs)
      call list_density % append(density * 0.0960_8)
      call list_names % append('42098.' // xs)
      call list_density % append(density * 0.2439_8)
      call list_names % append('42100.' // xs)
      call list_density % append(density * 0.0982_8)

    case ('ru')
      call list_names % append('44096.' // xs)
      call list_density % append(density * 0.0554_8)
      call list_names % append('44098.' // xs)
      call list_density % append(density * 0.0187_8)
      call list_names % append('44099.' // xs)
      call list_density % append(density * 0.1276_8)
      call list_names % append('44100.' // xs)
      call list_density % append(density * 0.1260_8)
      call list_names % append('44101.' // xs)
      call list_density % append(density * 0.1706_8)
      call list_names % append('44102.' // xs)
      call list_density % append(density * 0.3155_8)
      call list_names % append('44104.' // xs)
      call list_density % append(density * 0.1862_8)

    case ('rh')
      call list_names % append('45103.' // xs)
      call list_density % append(density)

    case ('pd')
      call list_names % append('46102.' // xs)
      call list_density % append(density * 0.0102_8)
      call list_names % append('46104.' // xs)
      call list_density % append(density * 0.1114_8)
      call list_names % append('46105.' // xs)
      call list_density % append(density * 0.2233_8)
      call list_names % append('46106.' // xs)
      call list_density % append(density * 0.2733_8)
      call list_names % append('46108.' // xs)
      call list_density % append(density * 0.2646_8)
      call list_names % append('46110.' // xs)
      call list_density % append(density * 0.1172_8)

    case ('ag')
      call list_names % append('47107.' // xs)
      call list_density % append(density * 0.51839_8)
      call list_names % append('47109.' // xs)
      call list_density % append(density * 0.48161_8)

    case ('cd')
      call list_names % append('48106.' // xs)
      call list_density % append(density * 0.0125_8)
      call list_names % append('48108.' // xs)
      call list_density % append(density * 0.0089_8)
      call list_names % append('48110.' // xs)
      call list_density % append(density * 0.1249_8)
      call list_names % append('48111.' // xs)
      call list_density % append(density * 0.1280_8)
      call list_names % append('48112.' // xs)
      call list_density % append(density * 0.2413_8)
      call list_names % append('48113.' // xs)
      call list_density % append(density * 0.1222_8)
      call list_names % append('48114.' // xs)
      call list_density % append(density * 0.2873_8)
      call list_names % append('48116.' // xs)
      call list_density % append(density * 0.0749_8)

    case ('in')
      call list_names % append('49113.' // xs)
      call list_density % append(density * 0.0429_8)
      call list_names % append('49115.' // xs)
      call list_density % append(density * 0.9571_8)

    case ('sn')
      call list_names % append('50112.' // xs)
      call list_density % append(density * 0.0097_8)
      call list_names % append('50114.' // xs)
      call list_density % append(density * 0.0066_8)
      call list_names % append('50115.' // xs)
      call list_density % append(density * 0.0034_8)
      call list_names % append('50116.' // xs)
      call list_density % append(density * 0.1454_8)
      call list_names % append('50117.' // xs)
      call list_density % append(density * 0.0768_8)
      call list_names % append('50118.' // xs)
      call list_density % append(density * 0.2422_8)
      call list_names % append('50119.' // xs)
      call list_density % append(density * 0.0859_8)
      call list_names % append('50120.' // xs)
      call list_density % append(density * 0.3258_8)
      call list_names % append('50122.' // xs)
      call list_density % append(density * 0.0463_8)
      call list_names % append('50124.' // xs)
      call list_density % append(density * 0.0579_8)

    case ('sb')
      call list_names % append('51121.' // xs)
      call list_density % append(density * 0.5721_8)
      call list_names % append('51123.' // xs)
      call list_density % append(density * 0.4279_8)

    case ('te')
      call list_names % append('52120.' // xs)
      call list_density % append(density * 0.0009_8)
      call list_names % append('52122.' // xs)
      call list_density % append(density * 0.0255_8)
      call list_names % append('52123.' // xs)
      call list_density % append(density * 0.0089_8)
      call list_names % append('52124.' // xs)
      call list_density % append(density * 0.0474_8)
      call list_names % append('52125.' // xs)
      call list_density % append(density * 0.0707_8)
      call list_names % append('52126.' // xs)
      call list_density % append(density * 0.1884_8)
      call list_names % append('52128.' // xs)
      call list_density % append(density * 0.3174_8)
      call list_names % append('52130.' // xs)
      call list_density % append(density * 0.3408_8)

    case ('i')
      call list_names % append('53127.' // xs)
      call list_density % append(density)

    case ('xe')
      call list_names % append('54124.' // xs)
      call list_density % append(density * 0.000952_8)
      call list_names % append('54126.' // xs)
      call list_density % append(density * 0.000890_8)
      call list_names % append('54128.' // xs)
      call list_density % append(density * 0.019102_8)
      call list_names % append('54129.' // xs)
      call list_density % append(density * 0.264006_8)
      call list_names % append('54130.' // xs)
      call list_density % append(density * 0.040710_8)
      call list_names % append('54131.' // xs)
      call list_density % append(density * 0.212324_8)
      call list_names % append('54132.' // xs)
      call list_density % append(density * 0.269086_8)
      call list_names % append('54134.' // xs)
      call list_density % append(density * 0.104357_8)
      call list_names % append('54136.' // xs)
      call list_density % append(density * 0.088573_8)

    case ('cs')
      call list_names % append('55133.' // xs)
      call list_density % append(density)

    case ('ba')
      call list_names % append('56130.' // xs)
      call list_density % append(density * 0.00106_8)
      call list_names % append('56132.' // xs)
      call list_density % append(density * 0.00101_8)
      call list_names % append('56134.' // xs)
      call list_density % append(density * 0.02417_8)
      call list_names % append('56135.' // xs)
      call list_density % append(density * 0.06592_8)
      call list_names % append('56136.' // xs)
      call list_density % append(density * 0.07854_8)
      call list_names % append('56137.' // xs)
      call list_density % append(density * 0.11232_8)
      call list_names % append('56138.' // xs)
      call list_density % append(density * 0.71698_8)

    case ('la')
      call list_names % append('57138.' // xs)
      call list_density % append(density * 0.0008881_8)
      call list_names % append('57139.' // xs)
      call list_density % append(density * 0.9991119_8)

    case ('ce')
      call list_names % append('58136.' // xs)
      call list_density % append(density * 0.00185_8)
      call list_names % append('58138.' // xs)
      call list_density % append(density * 0.00251_8)
      call list_names % append('58140.' // xs)
      call list_density % append(density * 0.88450_8)
      call list_names % append('58142.' // xs)
      call list_density % append(density * 0.11114_8)

    case ('pr')
      call list_names % append('59141.' // xs)
      call list_density % append(density)

    case ('nd')
      call list_names % append('60142.' // xs)
      call list_density % append(density * 0.27152_8)
      call list_names % append('60143.' // xs)
      call list_density % append(density * 0.12174_8)
      call list_names % append('60144.' // xs)
      call list_density % append(density * 0.23798_8)
      call list_names % append('60145.' // xs)
      call list_density % append(density * 0.08293_8)
      call list_names % append('60146.' // xs)
      call list_density % append(density * 0.17189_8)
      call list_names % append('60148.' // xs)
      call list_density % append(density * 0.05756_8)
      call list_names % append('60150.' // xs)
      call list_density % append(density * 0.05638_8)

    case ('sm')
      call list_names % append('62144.' // xs)
      call list_density % append(density * 0.0307_8)
      call list_names % append('62147.' // xs)
      call list_density % append(density * 0.1499_8)
      call list_names % append('62148.' // xs)
      call list_density % append(density * 0.1124_8)
      call list_names % append('62149.' // xs)
      call list_density % append(density * 0.1382_8)
      call list_names % append('62150.' // xs)
      call list_density % append(density * 0.0738_8)
      call list_names % append('62152.' // xs)
      call list_density % append(density * 0.2675_8)
      call list_names % append('62154.' // xs)
      call list_density % append(density * 0.2275_8)

    case ('eu')
      call list_names % append('63151.' // xs)
      call list_density % append(density * 0.4781_8)
      call list_names % append('63153.' // xs)
      call list_density % append(density * 0.5219_8)

    case ('gd')
      call list_names % append('64152.' // xs)
      call list_density % append(density * 0.0020_8)
      call list_names % append('64154.' // xs)
      call list_density % append(density * 0.0218_8)
      call list_names % append('64155.' // xs)
      call list_density % append(density * 0.1480_8)
      call list_names % append('64156.' // xs)
      call list_density % append(density * 0.2047_8)
      call list_names % append('64157.' // xs)
      call list_density % append(density * 0.1565_8)
      call list_names % append('64158.' // xs)
      call list_density % append(density * 0.2484_8)
      call list_names % append('64160.' // xs)
      call list_density % append(density * 0.2186_8)

    case ('tb')
      call list_names % append('65159.' // xs)
      call list_density % append(density)

    case ('dy')
      call list_names % append('66156.' // xs)
      call list_density % append(density * 0.00056_8)
      call list_names % append('66158.' // xs)
      call list_density % append(density * 0.00095_8)
      call list_names % append('66160.' // xs)
      call list_density % append(density * 0.02329_8)
      call list_names % append('66161.' // xs)
      call list_density % append(density * 0.18889_8)
      call list_names % append('66162.' // xs)
      call list_density % append(density * 0.25475_8)
      call list_names % append('66163.' // xs)
      call list_density % append(density * 0.24896_8)
      call list_names % append('66164.' // xs)
      call list_density % append(density * 0.28260_8)

    case ('ho')
      call list_names % append('67165.' // xs)
      call list_density % append(density)

    case ('er')
      call list_names % append('68162.' // xs)
      call list_density % append(density * 0.00139_8)
      call list_names % append('68164.' // xs)
      call list_density % append(density * 0.01601_8)
      call list_names % append('68166.' // xs)
      call list_density % append(density * 0.33503_8)
      call list_names % append('68167.' // xs)
      call list_density % append(density * 0.22869_8)
      call list_names % append('68168.' // xs)
      call list_density % append(density * 0.26978_8)
      call list_names % append('68170.' // xs)
      call list_density % append(density * 0.14910_8)

    case ('tm')
      call list_names % append('69169.' // xs)
      call list_density % append(density)

    case ('yb')
      call list_names % append('70168.' // xs)
      call list_density % append(density * 0.00123_8)
      call list_names % append('70170.' // xs)
      call list_density % append(density * 0.02982_8)
      call list_names % append('70171.' // xs)
      call list_density % append(density * 0.1409_8)
      call list_names % append('70172.' // xs)
      call list_density % append(density * 0.2168_8)
      call list_names % append('70173.' // xs)
      call list_density % append(density * 0.16103_8)
      call list_names % append('70174.' // xs)
      call list_density % append(density * 0.32026_8)
      call list_names % append('70176.' // xs)
      call list_density % append(density * 0.12996_8)

    case ('lu')
      call list_names % append('71175.' // xs)
      call list_density % append(density * 0.97401_8)
      call list_names % append('71176.' // xs)
      call list_density % append(density * 0.02599_8)

    case ('hf')
      call list_names % append('72174.' // xs)
      call list_density % append(density * 0.0016_8)
      call list_names % append('72176.' // xs)
      call list_density % append(density * 0.0526_8)
      call list_names % append('72177.' // xs)
      call list_density % append(density * 0.1860_8)
      call list_names % append('72178.' // xs)
      call list_density % append(density * 0.2728_8)
      call list_names % append('72179.' // xs)
      call list_density % append(density * 0.1362_8)
      call list_names % append('72180.' // xs)
      call list_density % append(density * 0.3508_8)

    case ('ta')
      call list_names % append('73180.' // xs)
      call list_density % append(density * 0.0001201_8)
      call list_names % append('73181.' // xs)
      call list_density % append(density * 0.9998799_8)

    case ('w')
      ! ENDF/B-VII.0 does not have W-180 so this may cause problems. However, it
      ! has been added as of ENDF/B-VII.1
      call list_names % append('74180.' // xs)
      call list_density % append(density * 0.0012_8)
      call list_names % append('74182.' // xs)
      call list_density % append(density * 0.2650_8)
      call list_names % append('74183.' // xs)
      call list_density % append(density * 0.1431_8)
      call list_names % append('74184.' // xs)
      call list_density % append(density * 0.3064_8)
      call list_names % append('74186.' // xs)
      call list_density % append(density * 0.2843_8)

    case ('re')
      call list_names % append('75185.' // xs)
      call list_density % append(density * 0.3740_8)
      call list_names % append('75187.' // xs)
      call list_density % append(density * 0.6260_8)

    case ('os')
      call list_names % append('76184.' // xs)
      call list_density % append(density * 0.0002_8)
      call list_names % append('76186.' // xs)
      call list_density % append(density * 0.0159_8)
      call list_names % append('76187.' // xs)
      call list_density % append(density * 0.0196_8)
      call list_names % append('76188.' // xs)
      call list_density % append(density * 0.1324_8)
      call list_names % append('76189.' // xs)
      call list_density % append(density * 0.1615_8)
      call list_names % append('76190.' // xs)
      call list_density % append(density * 0.2626_8)
      call list_names % append('76192.' // xs)
      call list_density % append(density * 0.4078_8)

    case ('ir')
      call list_names % append('77191.' // xs)
      call list_density % append(density * 0.373_8)
      call list_names % append('77193.' // xs)
      call list_density % append(density * 0.627_8)

    case ('pt')
      call list_names % append('78190.' // xs)
      call list_density % append(density * 0.00012_8)
      call list_names % append('78192.' // xs)
      call list_density % append(density * 0.00782_8)
      call list_names % append('78194.' // xs)
      call list_density % append(density * 0.3286_8)
      call list_names % append('78195.' // xs)
      call list_density % append(density * 0.3378_8)
      call list_names % append('78196.' // xs)
      call list_density % append(density * 0.2521_8)
      call list_names % append('78198.' // xs)
      call list_density % append(density * 0.07356_8)

    case ('au')
      call list_names % append('79197.' // xs)
      call list_density % append(density)

    case ('hg')
      call list_names % append('80196.' // xs)
      call list_density % append(density * 0.0015_8)
      call list_names % append('80198.' // xs)
      call list_density % append(density * 0.0997_8)
      call list_names % append('80199.' // xs)
      call list_density % append(density * 0.1687_8)
      call list_names % append('80200.' // xs)
      call list_density % append(density * 0.2310_8)
      call list_names % append('80201.' // xs)
      call list_density % append(density * 0.1318_8)
      call list_names % append('80202.' // xs)
      call list_density % append(density * 0.2986_8)
      call list_names % append('80204.' // xs)
      call list_density % append(density * 0.0687_8)

    case ('tl')
      call list_names % append('81203.' // xs)
      call list_density % append(density * 0.2952_8)
      call list_names % append('81205.' // xs)
      call list_density % append(density * 0.7048_8)

    case ('pb')
      call list_names % append('82204.' // xs)
      call list_density % append(density * 0.014_8)
      call list_names % append('82206.' // xs)
      call list_density % append(density * 0.241_8)
      call list_names % append('82207.' // xs)
      call list_density % append(density * 0.221_8)
      call list_names % append('82208.' // xs)
      call list_density % append(density * 0.524_8)

    case ('bi')
      call list_names % append('83209.' // xs)
      call list_density % append(density)

    case ('th')
      call list_names % append('90232.' // xs)
      call list_density % append(density)

    case ('pa')
      call list_names % append('91231.' // xs)
      call list_density % append(density)

    case ('u')
      call list_names % append('92234.' // xs)
      call list_density % append(density * 0.000054_8)
      call list_names % append('92235.' // xs)
      call list_density % append(density * 0.007204_8)
      call list_names % append('92238.' // xs)
      call list_density % append(density * 0.992742_8)

    case default
      message = "Cannot expand element: " // name
      call fatal_error()

    end select

  end subroutine expand_natural_element

end module input_xml
