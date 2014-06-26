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
  use xml_interface

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

    character(MAX_LINE_LEN) :: temp_str
    integer :: i
    integer :: n
    integer :: coeffs_reqd
    integer :: temp_int
    integer :: temp_int_array3(3)
    integer, allocatable :: temp_int_array(:)
    integer(8) :: temp_long
    integer :: n_tracks
    logical :: file_exists
    character(MAX_FILE_LEN) :: env_variable
    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: filename
    type(Node), pointer :: doc          => null()
    type(Node), pointer :: node_mode    => null()
    type(Node), pointer :: node_source  => null()
    type(Node), pointer :: node_dist    => null()
    type(Node), pointer :: node_cutoff  => null()
    type(Node), pointer :: node_entropy => null()
    type(Node), pointer :: node_ufs     => null()
    type(Node), pointer :: node_sp      => null()
    type(Node), pointer :: node_output  => null()
    type(Node), pointer :: node_verb    => null()

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

    ! Parse settings.xml file
    call open_xmldoc(doc, filename)

    ! Find cross_sections.xml file -- the first place to look is the
    ! settings.xml file. If no file is found there, then we check the
    ! CROSS_SECTIONS environment variable
    if (run_mode /= MODE_PLOTTING) then
      if (.not. check_for_node(doc, "cross_sections") .and. &
           run_mode /= MODE_PLOTTING) then
        ! No cross_sections.xml file specified in settings.xml, check
        ! environment variable
        call get_environment_variable("CROSS_SECTIONS", env_variable)
        if (len_trim(env_variable) == 0) then
          message = "No cross_sections.xml file was specified in settings.xml &
               &or in the CROSS_SECTIONS environment variable."
          call fatal_error()
        else
          path_cross_sections = trim(env_variable)
        end if
      else
        call get_node_value(doc, "cross_sections", path_cross_sections)
      end if
    end if

    ! Set output directory if a path has been specified on the <output_path>
    ! element
    if (check_for_node(doc, "output_path")) then
      call get_node_value(doc, "output_path", path_output)
      if (.not. ends_with(path_output, "/")) &
           path_output = trim(path_output) // "/"
    end if

    ! Make sure that either eigenvalue or fixed source was specified
    if (.not.check_for_node(doc, "eigenvalue") .and. &
         .not.check_for_node(doc, "fixed_source")) then
      message = "<eigenvalue> or <fixed_source> not specified."
      call fatal_error()
    end if

    ! Eigenvalue information
    if (check_for_node(doc, "eigenvalue")) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_EIGENVALUE

      ! Get pointer to eigenvalue XML block
      call get_node_ptr(doc, "eigenvalue", node_mode)

      ! Check number of particles
      if (.not.check_for_node(node_mode, "particles")) then
        message = "Need to specify number of particles per generation."
        call fatal_error()
      end if

      ! Get number of particles
      call get_node_value(node_mode, "particles", temp_long)

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = temp_long

      ! Copy batch and generation information
      call get_node_value(node_mode, "batches", n_batches)
      call get_node_value(node_mode, "inactive", n_inactive)
      n_active = n_batches - n_inactive
      if (check_for_node(node_mode, "generations_per_batch")) then
        call get_node_value(node_mode, "generations_per_batch", gen_per_batch)
      end if

      ! Allocate array for batch keff and entropy
      allocate(k_generation(n_batches*gen_per_batch))
      allocate(entropy(n_batches*gen_per_batch))
      entropy = ZERO
    end if

    ! Fixed source calculation information
    if (check_for_node(doc, "fixed_source")) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_FIXEDSOURCE

      ! Get pointer to fixed_source XML block
      call get_node_ptr(doc, "fixed_source", node_mode)

      ! Check number of particles
      if (.not.check_for_node(node_mode, "particles")) then
        message = "Need to specify number of particles per batch."
        call fatal_error()
      end if

      ! Get number of particles
      call get_node_value(node_mode, "particles", temp_long)

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = temp_long 

      ! Copy batch information
      call get_node_value(node_mode, "batches", n_batches)
      n_active = n_batches
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
    if (check_for_node(doc, "seed")) call get_node_value(doc, "seed", seed)

    ! Energy grid methods
    if (check_for_node(doc, "energy_grid")) then
      call get_node_value(doc, "energy_grid", temp_str)
    else
      temp_str = 'union'
    end if
    select case (trim(temp_str))
    case ('nuclide')
      grid_method = GRID_NUCLIDE
    case ('union')
      grid_method = GRID_UNION
    case ('lethargy')
      message = "Lethargy mapped energy grid not yet supported."
      call fatal_error()
    case default
      message = "Unknown energy grid method: " // trim(temp_str)
      call fatal_error()
    end select

    ! Verbosity
    if (check_for_node(doc, "verbosity")) then
      call get_node_ptr(doc, "verbosity", node_verb)
      call get_node_value(node_verb, "value", verbosity)
    end if

    ! Number of OpenMP threads
    if (check_for_node(doc, "threads")) then
#ifdef _OPENMP
      if (n_threads == NONE) then
        call get_node_value(doc, "threads", n_threads)
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

    ! Get pointer to source
    if (check_for_node(doc, "source")) then
      call get_node_ptr(doc, "source", node_source)
    else
      message = "No source specified in settings XML file."
      call fatal_error()
    end if

    ! Check for external source file
    if (check_for_node(node_source, "file")) then
      ! Copy path of source file
      call get_node_value(node_source, "file", path_source)

      ! Check if source file exists
      inquire(FILE=path_source, EXIST=file_exists)
      if (.not. file_exists) then
        message = "Binary source file '" // trim(path_source) // &
             "' does not exist!"
        call fatal_error()
      end if

    else

      ! Spatial distribution for external source
      if (check_for_node(node_source, "space")) then 

        ! Get pointer to spatial distribution
        call get_node_ptr(node_source, "space", node_dist) 

        ! Check for type of spatial distribution
        type = ''
        if (check_for_node(node_dist, "type")) &
             call get_node_value(node_dist, "type", type)
        call lower_case(type)
        select case (trim(type))
        case ('box')
          external_source % type_space = SRC_SPACE_BOX
          coeffs_reqd = 6
        case ('point')
          external_source % type_space = SRC_SPACE_POINT
          coeffs_reqd = 3
        case default
          message = "Invalid spatial distribution for external source: " &
              // trim(type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (check_for_node(node_dist, "parameters")) then
          n = get_arraysize_double(node_dist, "parameters")
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
          call get_node_array(node_dist, "parameters", &
               external_source % params_space)
        end if
      else
        message = "No spatial distribution specified for external source."
        call fatal_error()
      end if

      ! Determine external source angular distribution
      if (check_for_node(node_source, "angle")) then

        ! Get pointer to angular distribution
        call get_node_ptr(node_source, "angle", node_dist)

        ! Check for type of angular distribution
        type = ''
        if (check_for_node(node_dist, "type")) &
             call get_node_value(node_dist, "type", type)
        call lower_case(type)
        select case (trim(type))
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
               // trim(type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (check_for_node(node_dist, "parameters")) then
          n = get_arraysize_double(node_dist, "parameters")
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
          call get_node_array(node_dist, "parameters", &
               external_source % params_angle)
        end if
      else
        ! Set default angular distribution isotropic
        external_source % type_angle  = SRC_ANGLE_ISOTROPIC
      end if

      ! Determine external source energy distribution
      if (check_for_node(node_source, "energy")) then

        ! Get pointer to energy distribution
        call get_node_ptr(node_source, "energy", node_dist)

        ! Check for type of energy distribution
        type = ''
        if (check_for_node(node_dist, "type")) &
          call get_node_value(node_dist, "type", type)
        call lower_case(type)
        select case (trim(type))
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
               // trim(type)
          call fatal_error()
        end select

        ! Determine number of parameters specified
        if (check_for_node(node_dist, "parameters")) then
          n = get_arraysize_double(node_dist, "parameters")
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
          call get_node_array(node_dist, "parameters", &
               external_source % params_energy)
        end if
      else
        ! Set default energy distribution to Watt fission spectrum
        external_source % type_energy = SRC_ENERGY_WATT
        allocate(external_source % params_energy(2))
        external_source % params_energy = (/ 0.988_8, 2.249_8 /)
      end if
    end if

    ! Survival biasing
    if (check_for_node(doc, "survival_biasing")) then
      call get_node_value(doc, "survival_biasing", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
           survival_biasing = .true.
    end if

    ! Probability tables
    if (check_for_node(doc, "ptables")) then
      call get_node_value(doc, "ptables", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'false' .or. trim(temp_str) == '0') &
           urr_ptables_on = .false.
    end if

    ! Cutoffs
    if (check_for_node(doc, "cutoff")) then
      call get_node_ptr(doc, "cutoff", node_cutoff)
      call get_node_value(node_cutoff, "weight", weight_cutoff)
      call get_node_value(node_cutoff, "weight_avg", weight_survive)
    end if

    ! Particle trace
    if (check_for_node(doc, "trace")) then
      call get_node_array(doc, "trace", temp_int_array3)
      trace_batch    = temp_int_array3(1)
      trace_gen      = temp_int_array3(2)
      trace_particle = int(temp_int_array3(3), 8)
    end if

    ! Particle tracks
    if (check_for_node(doc, "track")) then
      ! Make sure that there are three values per particle
      n_tracks = get_arraysize_integer(doc, "track")
      if (mod(n_tracks, 3) /= 0) then
        message = "Number of integers specified in 'track' is not divisible &
             &by 3.  Please provide 3 integers per particle to be tracked."
        call fatal_error()
      end if

      ! Allocate space and get list of tracks
      allocate(temp_int_array(n_tracks))
      call get_node_array(doc, "track", temp_int_array)

      ! Reshape into track_identifiers
      allocate(track_identifiers(3, n_tracks/3))
      track_identifiers = reshape(temp_int_array, [3, n_tracks/3])
    end if

    ! Shannon Entropy mesh
    if (check_for_node(doc, "entropy")) then

      ! Get pointer to entropy node
      call get_node_ptr(doc, "entropy", node_entropy)

      ! Check to make sure enough values were supplied
      if (get_arraysize_double(node_entropy, "lower_left") /= 3) then
        message = "Need to specify (x,y,z) coordinates of lower-left corner &
             &of Shannon entropy mesh."
      elseif (get_arraysize_double(node_entropy, "upper_right") /= 3) then
        message = "Need to specify (x,y,z) coordinates of upper-right corner &
             &of Shannon entropy mesh."
      end if

      ! Allocate mesh object and coordinates on mesh
      allocate(entropy_mesh)
      allocate(entropy_mesh % lower_left(3))
      allocate(entropy_mesh % upper_right(3))

      ! Copy values
      call get_node_array(node_entropy, "lower_left", &
           entropy_mesh % lower_left)
      call get_node_array(node_entropy, "upper_right", &
           entropy_mesh % upper_right)

      ! Check on values provided
      if (.not. all(entropy_mesh % upper_right > entropy_mesh % lower_left)) then
        message = "Upper-right coordinate must be greater than lower-left &
             &coordinate for Shannon entropy mesh."
        call fatal_error()
      end if

      ! Check if dimensions were specified -- if not, they will be calculated
      ! automatically upon first entry into shannon_entropy
      if (check_for_node(node_entropy, "dimension")) then

        ! If so, make sure proper number of values were given
        if (get_arraysize_integer(node_entropy, "dimension") /= 3) then
          message = "Dimension of entropy mesh must be given as three &
               &integers."
          call fatal_error()
        end if

        ! Allocate dimensions
        entropy_mesh % n_dimension = 3
        allocate(entropy_mesh % dimension(3))

        ! Copy dimensions
        call get_node_array(node_entropy, "dimension", entropy_mesh % dimension)
      end if

      ! Turn on Shannon entropy calculation
      entropy_on = .true.
    end if

    ! Uniform fission source weighting mesh
    if (check_for_node(doc, "uniform_fs")) then

      ! Get pointer to ufs node
      call get_node_ptr(doc, "uniform_fs", node_ufs)

      ! Check to make sure enough values were supplied
      if (get_arraysize_double(node_ufs, "lower_left") /= 3) then
        message = "Need to specify (x,y,z) coordinates of lower-left corner &
             &of UFS mesh."
      elseif (get_arraysize_double(node_ufs, "upper_right") /= 3) then
        message = "Need to specify (x,y,z) coordinates of upper-right corner &
             &of UFS mesh."
      elseif (get_arraysize_integer(node_ufs, "dimension") /= 3) then
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
      call get_node_array(node_ufs, "dimension", ufs_mesh % dimension)

      ! Copy values
      call get_node_array(node_ufs, "lower_left", ufs_mesh % lower_left)
      call get_node_array(node_ufs, "upper_right", ufs_mesh % upper_right)

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
    if (check_for_node(doc, "state_point")) then

      ! Get pointer to state_point node
      call get_node_ptr(doc, "state_point", node_sp)

      ! Determine number of batches at which to store state points
      if (check_for_node(node_sp, "batches")) then
        n_state_points = get_arraysize_integer(node_sp, "batches")
      else
        n_state_points = 0
      end if

      if (n_state_points > 0) then
        ! User gave specific batches to write state points
        allocate(temp_int_array(n_state_points))
        call get_node_array(node_sp, "batches", temp_int_array)
        do i = 1, n_state_points
          call statepoint_batch % add(temp_int_array(i))
        end do
        deallocate(temp_int_array)
      elseif (check_for_node(node_sp, "interval")) then
        ! User gave an interval for writing state points
        call get_node_value(node_sp, "interval", temp_int)
        n_state_points = n_batches / temp_int 
        do i = 1, n_state_points
          call statepoint_batch % add(temp_int * i)
        end do
      else
        ! If neither were specified, write state point at last batch
        n_state_points = 1
        call statepoint_batch % add(n_batches)
      end if
    else
      ! If no <state_point> tag was present, by default write state point at
      ! last batch only
      n_state_points = 1
      call statepoint_batch % add(n_batches)
    end if

    ! Check if the user has specified to write source points
    if (check_for_node(doc, "source_point")) then

      ! Get pointer to source_point node
      call get_node_ptr(doc, "source_point", node_sp)

      ! Determine number of batches at which to store source points
      if (check_for_node(node_sp, "batches")) then
        n_source_points = get_arraysize_integer(node_sp, "batches") 
      else
        n_source_points = 0
      end if

      if (n_source_points > 0) then
        ! User gave specific batches to write source points
        allocate(temp_int_array(n_source_points))
        call get_node_array(node_sp, "batches", temp_int_array)
        do i = 1, n_source_points
          call sourcepoint_batch % add(temp_int_array(i))
        end do
        deallocate(temp_int_array)
      elseif (check_for_node(node_sp, "interval")) then
        ! User gave an interval for writing source points
        call get_node_value(node_sp, "interval", temp_int)
        n_source_points = n_batches / temp_int
        do i = 1, n_source_points
          call sourcepoint_batch % add(temp_int * i)
        end do
      else
        ! If neither were specified, write source points with state points
        n_source_points = n_state_points
        do i = 1, n_state_points
          call sourcepoint_batch % add(statepoint_batch % get_item(i))
        end do
      end if

      ! Check if the user has specified to write binary source file
      if (check_for_node(node_sp, "separate")) then
        call get_node_value(node_sp, "separate", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') source_separate = .true.
      end if
      if (check_for_node(node_sp, "write")) then
        call get_node_value(node_sp, "write", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'false' .or. &
             trim(temp_str) == '0') source_write = .false.
      end if
      if (check_for_node(node_sp, "overwrite_latest")) then
        call get_node_value(node_sp, "overwrite_latest", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') then
          source_latest = .true.
          source_separate = .true.
        end if
      end if
    else
      ! If no <source_point> tag was present, by default we keep source bank in
      ! statepoint file and write it out at statepoints intervals 
      source_separate = .false.
      n_source_points = n_state_points
      do i = 1, n_state_points
        call sourcepoint_batch % add(statepoint_batch % get_item(i))
      end do
    end if

    ! If source is not seperate and is to be written out in the statepoint file,
    ! make sure that the sourcepoint batch numbers are contained in the
    ! statepoint list
    if (.not. source_separate) then
      do i = 1, n_source_points
        if (.not. statepoint_batch % contains(sourcepoint_batch % &
            get_item(i))) then
          message = 'Sourcepoint batches are not a subset&
                    & of statepoint batches.'
          call fatal_error()
        end if
      end do
    end if

    ! Check if the user has specified to not reduce tallies at the end of every
    ! batch
    if (check_for_node(doc, "no_reduce")) then
      call get_node_value(doc, "no_reduce", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        reduce_tallies = .false.
    end if

    ! Check if the user has specified to use confidence intervals for
    ! uncertainties rather than standard deviations
    if (check_for_node(doc, "confidence_intervals")) then
      call get_node_value(doc, "confidence_intervals", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. &
           trim(temp_str) == '1') confidence_intervals = .true.
    end if

    ! Check for output options
    if (check_for_node(doc, "output")) then

      ! Get pointer to output node
      call get_node_ptr(doc, "output", node_output)

      ! Check for summary option
      if (check_for_node(node_output, "summary")) then
        call get_node_value(node_output, "summary", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') output_summary = .true.
      end if

      ! Check for cross sections option
      if (check_for_node(node_output, "cross_sections")) then
        call get_node_value(node_output, "cross_sections", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') output_xs = .true.
      end if

      ! Check for ASCII tallies output option
      if (check_for_node(node_output, "tallies")) then
        call get_node_value(node_output, "tallies", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'false' .or. &
             trim(temp_str) == '0') output_tallies = .false.
      end if
    end if

    ! Check for cmfd run
    if (check_for_node(doc, "run_cmfd")) then
      call get_node_value(doc, "run_cmfd", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') then
        cmfd_run = .true.
#ifndef PETSC
        if (master) then
          message = 'CMFD is not available, compile OpenMC with PETSc'
          call fatal_error()
        end if
#endif
      end if
    end if

    ! Natural element expansion option
    if (check_for_node(doc, "natural_elements")) then
      call get_node_value(doc, "natural_elements", temp_str)
      call lower_case(temp_str)
      select case (temp_str)
      case ('endf/b-vii.0')
        default_expand = ENDF_BVII0
      case ('endf/b-vii.1')
        default_expand = ENDF_BVII1
      case ('jeff-3.1.1')
        default_expand = JEFF_311
      case ('jeff-3.1.2')
        default_expand = JEFF_312
      case ('jeff-3.2')
        default_expand = JEFF_32
      case ('jendl-3.2')
        default_expand = JENDL_32
      case ('jendl-3.3')
        default_expand = JENDL_33
      case ('jendl-4.0')
        default_expand = JENDL_40
      case default
        message = "Unknown natural element expansion option: " // trim(temp_str)
        call fatal_error()
      end select
    end if

    ! Close settings XML file
    call close_xmldoc(doc)

  end subroutine read_settings_xml

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i, j, k, m
    integer :: n
    integer :: n_x, n_y, n_z
    integer :: universe_num
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    integer :: mid
    integer :: temp_int_array3(3)
    integer, allocatable :: temp_int_array(:)
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    type(Cell),    pointer :: c => null()
    type(Surface), pointer :: s => null()
    type(Lattice), pointer :: lat => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_cell => null()
    type(Node), pointer :: node_surf => null()
    type(Node), pointer :: node_lat => null()
    type(NodeList), pointer :: node_cell_list => null()
    type(NodeList), pointer :: node_surf_list => null()
    type(NodeList), pointer :: node_lat_list => null()

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
    call open_xmldoc(doc, filename)

    ! Get pointer to list of XML <cell>
    call get_node_list(doc, "cell", node_cell_list)

    ! Get number of <cell> tags
    n_cells = get_list_size(node_cell_list)

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

      ! Get pointer to i-th cell node
      call get_list_item(node_cell_list, i, node_cell)

      ! Copy data into cells
      if (check_for_node(node_cell, "id")) then
        call get_node_value(node_cell, "id", c % id)
      else
        message = "Must specify id of cell in geometry XML file."
        call fatal_error()
      end if
      if (check_for_node(node_cell, "universe")) then
        call get_node_value(node_cell, "universe", c % universe)
      else
        c % universe = NONE
      end if
      if (check_for_node(node_cell, "fill")) then
        call get_node_value(node_cell, "fill", c % fill)
      else
        c % fill = NONE
      end if

      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has_key(c % id)) then
        message = "Two or more cells use the same unique ID: " // to_str(c % id)
        call fatal_error()
      end if

      ! Read material
      word = ''
      if (check_for_node(node_cell, "material")) &
        call get_node_value(node_cell, "material", word)
      call lower_case(word)
      select case(word)
      case ('void')
        c % material = MATERIAL_VOID

      case ('')
        ! This case is called if no material was specified
        c % material = NONE

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
      if (.not. check_for_node(node_cell, "surfaces")) then
        message = "No surfaces specified for cell " // &
             trim(to_str(c % id))
        call fatal_error()
      end if

      ! Allocate array for surfaces and copy
      n = get_arraysize_integer(node_cell, "surfaces")
      c % n_surfaces = n
      allocate(c % surfaces(n))
      call get_node_array(node_cell, "surfaces", c % surfaces)

      ! Rotation matrix
      if (check_for_node(node_cell, "rotation")) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a rotation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of rotation parameters
        n = get_arraysize_double(node_cell, "rotation")
        if (n /= 3) then
          message = "Incorrect number of rotation parameters on cell " // &
               to_str(c % id)
          call fatal_error()
        end if

        ! Copy rotation angles in x,y,z directions
        call get_node_array(node_cell, "rotation", temp_int_array3)
        phi   = -temp_int_array3(1) * PI/180.0_8
        theta = -temp_int_array3(2) * PI/180.0_8
        psi   = -temp_int_array3(3) * PI/180.0_8

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
      if (check_for_node(node_cell, "translation")) then
        ! Translations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          message = "Cannot apply a translation to cell " // trim(to_str(&
               c % id)) // " because it is not filled with another universe"
          call fatal_error()
        end if

        ! Read number of translation parameters
        n = get_arraysize_double(node_cell, "translation")
        if (n /= 3) then
          message = "Incorrect number of translation parameters on cell " &
               // to_str(c % id)
          call fatal_error()
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        call get_node_array(node_cell, "translation", c % translation)
      end if

      ! Add cell to dictionary
      call cell_dict % add_key(c % id, i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      universe_num = c % universe
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

    ! get pointer to list of xml <surface>
    call get_node_list(doc, "surface", node_surf_list)

    ! Get number of <surface> tags
    n_surfaces = get_list_size(node_surf_list)

    ! Check for no surfaces
    if (n_surfaces == 0) then
      message = "No surfaces found in geometry.xml!"
      call fatal_error()
    end if

    ! Allocate cells array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      s => surfaces(i)

      ! Get pointer to i-th surface node
      call get_list_item(node_surf_list, i, node_surf)

      ! Copy data into cells
      if (check_for_node(node_surf, "id")) then
        call get_node_value(node_surf, "id", s % id)
      else
        message = "Must specify id of surface in geometry XML file."
        call fatal_error()
      end if

      ! Check to make sure 'id' hasn't been used
      if (surface_dict % has_key(s % id)) then
        message = "Two or more surfaces use the same unique ID: " // &
             to_str(s % id)
        call fatal_error()
      end if

      ! Copy and interpret surface type
      word = ''
      if (check_for_node(node_surf, "type")) &
        call get_node_value(node_surf, "type", word)
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
        message = "Invalid surface type: " // trim(word)
        call fatal_error()
      end select

      ! Check to make sure that the proper number of coefficients
      ! have been specified for the given type of surface. Then copy
      ! surface coordinates.

      n = get_arraysize_double(node_surf, "coeffs")
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
        call get_node_array(node_surf, "coeffs", s % coeffs)
      end if

      ! Boundary conditions
      word = ''
      if (check_for_node(node_surf, "boundary")) &
        call get_node_value(node_surf, "boundary", word)
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

    ! Get pointer to list of XML <lattice>
    call get_node_list(doc, "lattice", node_lat_list)

    ! Allocate lattices array
    n_lattices = get_list_size(node_lat_list)
    allocate(lattices(n_lattices))

    do i = 1, n_lattices
      lat => lattices(i)

      ! Get pointer to i-th lattice
      call get_list_item(node_lat_list, i, node_lat)

      ! ID of lattice
      if (check_for_node(node_lat, "id")) then
        call get_node_value(node_lat, "id", lat % id)
      else
        message = "Must specify id of lattice in geometry XML file."
        call fatal_error()
      end if

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        message = "Two or more lattices use the same unique ID: " // &
             to_str(lat % id)
        call fatal_error()
      end if

      ! Read lattice type
      word = ''
      if (check_for_node(node_lat, "type")) &
        call get_node_value(node_lat, "type", word)
      call lower_case(word)
      select case (trim(word))
      case ('rect', 'rectangle', 'rectangular')
        lat % type = LATTICE_RECT
      case ('hex', 'hexagon', 'hexagonal')
        lat % type = LATTICE_HEX
      case default
        message = "Invalid lattice type: " // trim(word)
        call fatal_error()
      end select

      ! Read number of lattice cells in each dimension
      n = get_arraysize_integer(node_lat, "dimension")
      if (n /= 2 .and. n /= 3) then
        message = "Lattice must be two or three dimensions."
        call fatal_error()
      end if

      lat % n_dimension = n
      allocate(lat % dimension(n))
      call get_node_array(node_lat, "dimension", lat % dimension)

      ! Read lattice lower-left location
      if (size(lat % dimension) /= &
          get_arraysize_double(node_lat, "lower_left")) then
        message = "Number of entries on <lower_left> must be the same as &
             &the number of entries on <dimension>."
        call fatal_error()
      end if

      allocate(lat % lower_left(n))
      call get_node_array(node_lat, "lower_left", lat % lower_left)

      ! Read lattice widths
      if (size(lat % dimension) /= &
          get_arraysize_double(node_lat, "width")) then
        message = "Number of entries on <width> must be the same as &
             &the number of entries on <lower_left>."
        call fatal_error()
      end if

      allocate(lat % width(n))
      call get_node_array(node_lat, "width", lat % width)

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
      n = get_arraysize_integer(node_lat, "universes")
      if (n /= n_x*n_y*n_z) then
        message = "Number of universes on <universes> does not match size of &
             &lattice " // trim(to_str(lat % id)) // "."
        call fatal_error()
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = &
               temp_int_array(j + n_x*k + n_x*n_y*(m-1))
          end do
        end do
      end do
      deallocate(temp_int_array)

      ! Read material for area outside lattice
      lat % outside = MATERIAL_VOID
      if (check_for_node(node_lat, "outside")) then
        call get_node_value(node_lat, "outside", mid)
        if (mid == 0 .or. mid == MATERIAL_VOID) then
          lat % outside = MATERIAL_VOID
        else
          lat % outside = mid
        end if
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, i)

    end do

    ! Close geometry XML file
    call close_xmldoc(doc)

  end subroutine read_geometry_xml

!===============================================================================
! READ_MATERIAL_XML reads data from a materials.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_materials_xml()

    integer :: i             ! loop index for materials
    integer :: j             ! loop index for nuclides
    integer :: n             ! number of nuclides
    integer :: n_sab         ! number of sab tables for a material
    integer :: index_list    ! index in xs_listings array
    integer :: index_nuclide ! index in nuclides
    integer :: index_sab     ! index in sab_tables
    real(8) :: val           ! value entered for density
    real(8) :: temp_dble     ! temporary double prec. real
    logical :: file_exists   ! does materials.xml exist?
    logical :: sum_density   ! density is taken to be sum of nuclide densities
    character(12) :: name    ! name of isotope, e.g. 92235.03c
    character(12) :: alias   ! alias of nuclide, e.g. U-235.03c
    character(MAX_WORD_LEN) :: units    ! units on density
    character(MAX_LINE_LEN) :: filename ! absolute path to materials.xml
    character(MAX_LINE_LEN) :: temp_str ! temporary string when reading
    type(ListChar) :: list_names   ! temporary list of nuclide names
    type(ListReal) :: list_density ! temporary list of nuclide densities
    type(Material),    pointer :: mat => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mat => null()
    type(Node), pointer :: node_dens => null()
    type(Node), pointer :: node_nuc => null()
    type(Node), pointer :: node_ele => null()
    type(Node), pointer :: node_sab => null()
    type(NodeList), pointer :: node_mat_list => null()
    type(NodeList), pointer :: node_nuc_list => null()
    type(NodeList), pointer :: node_ele_list => null()
    type(NodeList), pointer :: node_sab_list => null()

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
    default_xs = ""

    ! Parse materials.xml file
    call open_xmldoc(doc, filename)

    ! Copy default cross section if present
    if (check_for_node(doc, "default_xs")) &
      call get_node_value(doc, "default_xs", default_xs)

    ! Get pointer to list of XML <material>
    call get_node_list(doc, "material", node_mat_list)

    ! Allocate cells array
    n_materials = get_list_size(node_mat_list)
    allocate(materials(n_materials))

    ! Initialize count for number of nuclides/S(a,b) tables
    index_nuclide = 0
    index_sab = 0

    do i = 1, n_materials
      mat => materials(i)

      ! Get pointer to i-th material node
      call get_list_item(node_mat_list, i, node_mat)

      ! Copy material id
      if (check_for_node(node_mat, "id")) then
        call get_node_value(node_mat, "id", mat % id)
      else
        message = "Must specify id of material in materials XML file"
        call fatal_error()
      end if

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

      ! Get pointer to density element
      if (check_for_node(node_mat, "density")) then
        call get_node_ptr(node_mat, "density", node_dens)
      else
        message = "Must specify density element in material " // &
                  trim(to_str(mat % id))
        call fatal_error()
      end if

      ! Initialize value to zero
      val = ZERO

      ! Copy units
      call get_node_value(node_dens, "units", units)

      if (units == 'sum') then
        ! If the user gave the units as 'sum', then the total density of the
        ! material is taken to be the sum of the atom fractions listed on the
        ! nuclides

        sum_density = .true.

      else
        ! Copy value
        call get_node_value(node_dens, "value", val)

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
          message = "Unkwown units '" // trim(units) &
               // "' specified on material " // trim(to_str(mat % id))
          call fatal_error()
        end select
      end if

      ! =======================================================================
      ! READ AND PARSE <nuclide> TAGS

      ! Check to ensure material has at least one nuclide
      if (.not. check_for_node(node_mat, "nuclide") .and. &
           .not. check_for_node(node_mat, "element")) then
        message = "No nuclides or natural elements specified on material " // &
             trim(to_str(mat % id))
        call fatal_error()
      end if

      ! Get pointer list of XML <nuclide>
      call get_node_list(node_mat, "nuclide", node_nuc_list)

      ! Create list of nuclides based on those specified plus natural elements
      INDIVIDUAL_NUCLIDES: do j = 1, get_list_size(node_nuc_list)
        ! Combine nuclide identifier and cross section and copy into names
        call get_list_item(node_nuc_list, j, node_nuc)

        ! Check for empty name on nuclide
        if (.not.check_for_node(node_nuc, "name")) then
          message = "No name specified on nuclide in material " // &
               trim(to_str(mat % id))
          call fatal_error()
        end if

        ! Check for cross section
        if (.not.check_for_node(node_nuc, "xs")) then
          if (default_xs == '') then
            message = "No cross section specified for nuclide in material " &
                 // trim(to_str(mat % id))
            call fatal_error()
          else
            name = trim(default_xs)
          end if
        end if

        ! store full name
        call get_node_value(node_nuc, "name", temp_str)
        if (check_for_node(node_nuc, "xs")) &
          call get_node_value(node_nuc, "xs", name)
        name = trim(temp_str) // "." // trim(name)

        ! save name and density to list
        call list_names % append(name)

        ! Check if no atom/weight percents were specified or if both atom and
        ! weight percents were specified
        if (.not.check_for_node(node_nuc, "ao") .and. &
            .not.check_for_node(node_nuc, "wo")) then
          message = "No atom or weight percent specified for nuclide " // &
               trim(name)
          call fatal_error()
        elseif (check_for_node(node_nuc, "ao") .and. &
                check_for_node(node_nuc, "wo")) then
          message = "Cannot specify both atom and weight percents for a &
               &nuclide: " // trim(name)
          call fatal_error()
        end if

        ! Copy atom/weight percents
        if (check_for_node(node_nuc, "ao")) then
          call get_node_value(node_nuc, "ao", temp_dble)
          call list_density % append(temp_dble)
        else
          call get_node_value(node_nuc, "wo", temp_dble)
          call list_density % append(-temp_dble)
        end if
      end do INDIVIDUAL_NUCLIDES

      ! =======================================================================
      ! READ AND PARSE <element> TAGS

      ! Get pointer list of XML <element>
      call get_node_list(node_mat, "element", node_ele_list)

      NATURAL_ELEMENTS: do j = 1, get_list_size(node_ele_list)
        call get_list_item(node_ele_list, j, node_ele)

        ! Check for empty name on natural element
        if (.not.check_for_node(node_ele, "name")) then
          message = "No name specified on nuclide in material " // &
               trim(to_str(mat % id))
          call fatal_error()
        end if
        call get_node_value(node_ele, "name", name)

        ! Check for cross section
        if (check_for_node(node_ele, "xs")) then
          call get_node_value(node_ele, "xs", temp_str)
        else
          if (default_xs == '') then
            message = "No cross section specified for nuclide in material " &
                 // trim(to_str(mat % id))
            call fatal_error()
          else
            temp_str = trim(default_xs)
          end if
        end if

        ! Check if no atom/weight percents were specified or if both atom and
        ! weight percents were specified
        if (.not.check_for_node(node_ele, "ao") .and. &
            .not.check_for_node(node_ele, "wo")) then
          message = "No atom or weight percent specified for element " // &
               trim(name)
          call fatal_error()
        elseif (check_for_node(node_ele, "ao") .and. &
                check_for_node(node_ele, "wo")) then
          message = "Cannot specify both atom and weight percents for a &
               &element: " // trim(name)
          call fatal_error()
        end if

        ! Expand element into naturally-occurring isotopes
        if (check_for_node(node_ele, "ao")) then
          call get_node_value(node_ele, "ao", temp_dble)
          call expand_natural_element(name, temp_str, temp_dble, &
               list_names, list_density)
        else
          call get_node_value(node_ele, "wo", temp_dble)
          call expand_natural_element(name, temp_str, -temp_dble, &
               list_names, list_density)
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
        name = trim(list_names % get_item(j))
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

      ! Get pointer list to XML <sab>
      call get_node_list(node_mat, "sab", node_sab_list)

      n_sab = get_list_size(node_sab_list)
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
          call get_list_item(node_sab_list, j, node_sab)

          ! Determine name of S(a,b) table
          if (.not.check_for_node(node_sab, "name") .or. &
              .not.check_for_node(node_sab, "xs")) then
            message = "Need to specify <name> and <xs> for S(a,b) table."
            call fatal_error()
          end if
          call get_node_value(node_sab, "name", name)
          call get_node_value(node_sab, "xs", temp_str)
          name = trim(name) // "." // trim(temp_str)
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

          ! If this S(a,b) table hasn't been encountered yet, we need to add its
          ! name and alias to the sab_dict
          if (.not. sab_dict % has_key(name)) then
            index_sab = index_sab + 1
            mat % i_sab_tables(j) = index_sab
            call sab_dict % add_key(name, index_sab)
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

    ! Close materials XML file
    call close_xmldoc(doc)

  end subroutine read_materials_xml

!===============================================================================
! READ_TALLIES_XML reads data from a tallies.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_tallies_xml()

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
    integer :: iarray3(3)    ! temporary integer array
    logical :: file_exists   ! does tallies.xml file exist?
    real(8) :: rarray3(3)    ! temporary double prec. array
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    type(ElemKeyValueCI), pointer :: pair_list => null()
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter), allocatable :: filters(:) ! temporary filters
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mesh => null()
    type(Node), pointer :: node_tal => null()
    type(Node), pointer :: node_filt => null()
    type(NodeList), pointer :: node_mesh_list => null()
    type(NodeList), pointer :: node_tal_list => null()
    type(NodeList), pointer :: node_filt_list => null()

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
    call open_xmldoc(doc, filename)

    ! ==========================================================================
    ! DETERMINE SIZE OF ARRAYS AND ALLOCATE

    ! Get pointer list to XML <mesh>
    call get_node_list(doc, "mesh", node_mesh_list)

    ! Get pointer list to XML <tally>
    call get_node_list(doc, "tally", node_tal_list)

    ! Check for user meshes
    n_user_meshes = get_list_size(node_mesh_list)
    if (cmfd_run) then
      n_meshes = n_user_meshes + n_cmfd_meshes
    else
      n_meshes = n_user_meshes
    end if

    ! Allocate mesh array
    if (n_meshes > 0) allocate(meshes(n_meshes))

    ! Check for user tallies
    n_user_tallies = get_list_size(node_tal_list)
    if (n_user_tallies == 0) then
      message = "No tallies present in tallies.xml file!"
      call warning()
    end if

    ! Allocate tally array
    if (n_user_tallies > 0) then
      call add_tallies("user", n_user_tallies)
    end if

    ! Check for <assume_separate> setting
    if (check_for_node(doc, "assume_separate")) then
      call get_node_value(doc, "assume_separate", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        assume_separate = .true.
    end if

    ! ==========================================================================
    ! READ MESH DATA

    do i = 1, n_user_meshes
      m => meshes(i)

      ! Get pointer to mesh node
      call get_list_item(node_mesh_list, i, node_mesh)

      ! Copy mesh id
      if (check_for_node(node_mesh, "id")) then
        call get_node_value(node_mesh, "id", m % id)
      else
        message = "Must specify id for mesh in tally XML file."
        call fatal_error()
      end if

      ! Check to make sure 'id' hasn't been used
      if (mesh_dict % has_key(m % id)) then
        message = "Two or more meshes use the same unique ID: " // &
             to_str(m % id)
        call fatal_error()
      end if

      ! Read mesh type
      temp_str = ''
      if (check_for_node(node_mesh, "type")) &
        call get_node_value(node_mesh, "type", temp_str)
      call lower_case(temp_str)
      select case (trim(temp_str))
      case ('rect', 'rectangle', 'rectangular')
        m % type = LATTICE_RECT
      case ('hex', 'hexagon', 'hexagonal')
        m % type = LATTICE_HEX
      case default
        message = "Invalid mesh type: " // trim(temp_str)
        call fatal_error()
      end select

      ! Determine number of dimensions for mesh
      n = get_arraysize_integer(node_mesh, "dimension")
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
      call get_node_array(node_mesh, "dimension", iarray3(1:n))
      if (any(iarray3(1:n) <= 0)) then
        message = "All entries on the <dimension> element for a tally mesh &
             &must be positive."
        call fatal_error()
      end if

      ! Read dimensions in each direction
      m % dimension = iarray3(1:n)

      ! Read mesh lower-left corner location
      if (m % n_dimension /= get_arraysize_double(node_mesh, "lower_left")) then
        message = "Number of entries on <lower_left> must be the same as &
             &the number of entries on <dimension>."
        call fatal_error()
      end if
      call get_node_array(node_mesh, "lower_left", m % lower_left)

      ! Make sure both upper-right or width were specified
      if (check_for_node(node_mesh, "upper_right") .and. &
          check_for_node(node_mesh, "width")) then
        message = "Cannot specify both <upper_right> and <width> on a &
             &tally mesh."
        call fatal_error()
      end if

      ! Make sure either upper-right or width was specified
      if (.not.check_for_node(node_mesh, "upper_right") .and. &
          .not.check_for_node(node_mesh, "width")) then
        message = "Must specify either <upper_right> and <width> on a &
             &tally mesh."
        call fatal_error()
      end if

      if (check_for_node(node_mesh, "width")) then
        ! Check to ensure width has same dimensions
        if (get_arraysize_double(node_mesh, "width") /= &
            get_arraysize_double(node_mesh, "lower_left")) then
          message = "Number of entries on <width> must be the same as the &
               &number of entries on <lower_left>."
          call fatal_error()
        end if

        ! Check for negative widths
        call get_node_array(node_mesh, "width", rarray3(1:n))
        if (any(rarray3(1:n) < ZERO)) then
          message = "Cannot have a negative <width> on a tally mesh."
          call fatal_error()
        end if

        ! Set width and upper right coordinate
        m % width = rarray3(1:n)
        m % upper_right = m % lower_left + m % dimension * m % width

      elseif (check_for_node(node_mesh, "upper_right")) then
        ! Check to ensure width has same dimensions
        if (get_arraysize_double(node_mesh, "upper_right") /= &
            get_arraysize_double(node_mesh, "lower_left")) then
          message = "Number of entries on <upper_right> must be the same as &
               &the number of entries on <lower_left>."
          call fatal_error()
        end if

        ! Check that upper-right is above lower-left
        call get_node_array(node_mesh, "upper_right", rarray3(1:n))
        if (any(rarray3(1:n) < m % lower_left)) then
          message = "The <upper_right> coordinates must be greater than the &
               &<lower_left> coordinates on a tally mesh."
          call fatal_error()
        end if

        ! Set width and upper right coordinate
        m % upper_right = rarray3(1:n)
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

      ! Get pointer to tally xml node
      call get_list_item(node_tal_list, i, node_tal)

      ! Set tally type to volume by default
      t % type = TALLY_VOLUME

      ! It's desirable to use a track-length esimator for tallies since
      ! generally more events will score to the tally, reducing the
      ! variance. However, for tallies that require information on
      ! post-collision parameters (e.g. tally with an energyout filter) the
      ! analog esimator must be used.

      t % estimator = ESTIMATOR_TRACKLENGTH

      ! Copy material id
      if (check_for_node(node_tal, "id")) then
        call get_node_value(node_tal, "id", t % id)
      else
        message = "Must specify id for tally in tally XML file."
        call fatal_error()
      end if

      ! Check to make sure 'id' hasn't been used
      if (tally_dict % has_key(t % id)) then
        message = "Two or more tallies use the same unique ID: " // &
             to_str(t % id)
        call fatal_error()
      end if

      ! Copy tally label
      t % label = ''
      if (check_for_node(node_tal, "label")) &
        call get_node_value(node_tal, "label", t % label)

      ! =======================================================================
      ! READ DATA FOR FILTERS

      ! In older versions, tally filters were specified with a <filters>
      ! element followed by sub-elements <cell>, <mesh>, etc. This checks for
      ! the old format and if it is present, raises an error

!     if (get_number_nodes(node_tal, "filters") > 0) then
!       message = "Tally filters should be specified with multiple <filter> &
!            &elements. Did you forget to change your <filters> element?"
!       call fatal_error()
!     end if

      ! Get pointer list to XML <filter> and get number of filters
      call get_node_list(node_tal, "filter", node_filt_list)
      n_filters = get_list_size(node_filt_list)

      if (n_filters /= 0) then

        ! Allocate filters array
        t % n_filters = n_filters
        allocate(t % filters(n_filters))

        READ_FILTERS: do j = 1, n_filters
          ! Get pointer to filter xml node
          call get_list_item(node_filt_list, j, node_filt)

          ! Convert filter type to lower case
          temp_str = ''
          if (check_for_node(node_filt, "type")) &
            call get_node_value(node_filt, "type", temp_str)
          call lower_case(temp_str)

          ! Determine number of bins
          if (check_for_node(node_filt, "bins")) then
            if (trim(temp_str) == 'energy' .or. &
                trim(temp_str) == 'energyout') then
              n_words = get_arraysize_double(node_filt, "bins")
            else
              n_words = get_arraysize_integer(node_filt, "bins")
            end if
          else
            message = "Bins not set in filter on tally " // trim(to_str(t % id))
            call fatal_error()
          end if

          ! Determine type of filter
          select case (temp_str)
          case ('cell')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % int_bins)

          case ('cellborn')
            ! Set type of filter
            t % filters(j) % type = FILTER_CELLBORN

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % int_bins)

          case ('material')
            ! Set type of filter
            t % filters(j) % type = FILTER_MATERIAL

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % int_bins)

          case ('universe')
            ! Set type of filter
            t % filters(j) % type = FILTER_UNIVERSE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % int_bins)

          case ('surface')
            message = "Surface filter is not yet supported!"
            call fatal_error()

            ! Set type of filter
            t % filters(j) % type = FILTER_SURFACE

            ! Set number of bins
            t % filters(j) % n_bins = n_words

            ! Allocate and store bins
            allocate(t % filters(j) % int_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % int_bins)

          case ('mesh')
            ! Set type of filter
            t % filters(j) % type = FILTER_MESH

            ! Check to make sure multiple meshes weren't given
            if (n_words /= 1) then
              message = "Can only have one mesh filter specified."
              call fatal_error()
            end if

            ! Determine id of mesh
            call get_node_value(node_filt, "bins", id)

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
            call get_node_array(node_filt, "bins", t % filters(j) % real_bins)

          case ('energyout')
            ! Set type of filter
            t % filters(j) % type = FILTER_ENERGYOUT

            ! Set number of bins
            t % filters(j) % n_bins = n_words - 1

            ! Allocate and store bins
            allocate(t % filters(j) % real_bins(n_words))
            call get_node_array(node_filt, "bins", t % filters(j) % real_bins)

            ! Set to analog estimator
            t % estimator = ESTIMATOR_ANALOG

          case default
            ! Specified tally filter is invalid, raise error
            message = "Unknown filter type '" // & 
                 trim(temp_str) // "' on tally " // &
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

      if (check_for_node(node_tal, "nuclides")) then

        ! Allocate a temporary string array for nuclides and copy values over
        allocate(sarray(get_arraysize_string(node_tal, "nuclides")))
        call get_node_array(node_tal, "nuclides", sarray)

        if (trim(sarray(1)) == 'all') then
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
          n_words = get_arraysize_string(node_tal, "nuclides") 
          allocate(t % nuclide_bins(n_words))
          do j = 1, n_words
            ! Check if total material was specified
            if (trim(sarray(j)) == 'total') then
              t % nuclide_bins(j) = -1
              cycle
            end if

            ! Check if xs specifier was given
            if (ends_with(sarray(j), 'c')) then
              word = sarray(j)
            else
              if (default_xs == '') then
                ! No default cross section specified, search through nuclides
                pair_list => nuclide_dict % keys()
                do while (associated(pair_list))
                  if (starts_with(pair_list % key, &
                       sarray(j))) then
                    word = pair_list % key(1:150)
                    exit
                  end if
                  
                  ! Advance to next
                  pair_list => pair_list % next
                end do

                ! Check if no nuclide was found
                if (.not. associated(pair_list)) then
                  message = "Could not find the nuclide " // trim(&
                       sarray(j)) // " specified in tally " &
                       // trim(to_str(t % id)) // " in any material."
                  call fatal_error()
                end if
                deallocate(pair_list)
              else
                ! Set nuclide to default xs
                word = trim(sarray(j)) // "." // default_xs
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

        ! Deallocate temporary string array
        deallocate(sarray)

      else
        ! No <nuclides> were specified -- create only one bin will be added
        ! for the total material.
        allocate(t % nuclide_bins(1))
        t % nuclide_bins(1) = -1
        t % n_nuclide_bins = 1
      end if

      ! =======================================================================
      ! READ DATA FOR SCORES

      if (check_for_node(node_tal, "scores")) then
        ! Loop through scores and determine if a scatter-p# input was used
        ! to allow for proper pre-allocating of t % score_bins
        ! This scheme allows multiple scatter-p# to be requested by the user
        ! if so desired
        n_words = get_arraysize_string(node_tal, "scores")
        allocate(sarray(n_words))
        call get_node_array(node_tal, "scores", sarray)
        n_new = 0
        do j = 1, n_words
          call lower_case(sarray(j))
          ! Find if scores(j) is of the form 'scatter-p'
          ! If so, get the number and do a select case on that.
          score_name = trim(sarray(j))
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
              sarray(j) = SCATT_ORDER_MAX_PNSTR
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
          score_name = sarray(l)
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
          
          select case (trim(score_name))
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
                     trim(sarray(j))
                call fatal_error()
              end if

            else
              ! Specified score was not an integer
              message = "Unknown scoring function: " // &
                   trim(sarray(j))
              call fatal_error()
            end if

          end select
        end do
        t % n_score_bins = n_scores
        t % n_user_score_bins = n_words

        ! Deallocate temporary string array of scores
        deallocate(sarray)
      else
        message = "No <scores> specified on tally " // trim(to_str(t % id)) &
             // "."
        call fatal_error()
      end if

      ! =======================================================================
      ! SET TALLY ESTIMATOR

      ! Check if user specified estimator
      if (check_for_node(node_tal, "estimator")) then
        temp_str = ''
        call get_node_value(node_tal, "estimator", temp_str)
        select case(trim(temp_str))
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
          message = "Invalid estimator '" // trim(temp_str) &
               // "' on tally " // to_str(t % id)
          call fatal_error()
        end select
      end if

      ! Add tally to dictionary
      call tally_dict % add_key(t % id, i)

    end do READ_TALLIES

    ! Close XML document
    call close_xmldoc(doc)

  end subroutine read_tallies_xml

!===============================================================================
! READ_PLOTS_XML reads data from a plots.xml file
!===============================================================================

  subroutine read_plots_xml()

    integer i, j
    integer n_cols, col_id, n_comp, n_masks
    integer, allocatable :: iarray(:)
    logical :: file_exists              ! does plots.xml file exist?
    character(MAX_LINE_LEN) :: filename ! absolute path to plots.xml
    character(MAX_LINE_LEN) :: temp_str
    type(ObjectPlot), pointer :: pl => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_plot => null()
    type(Node), pointer :: node_col => null()
    type(Node), pointer :: node_mask => null()
    type(NodeList), pointer :: node_plot_list => null()
    type(NodeList), pointer :: node_col_list => null()
    type(NodeList), pointer :: node_mask_list => null()

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
    call open_xmldoc(doc, filename)

    ! Get list pointer to XML <plot>
    call get_node_list(doc, "plot", node_plot_list)

    ! Allocate plots array
    n_plots = get_list_size(node_plot_list)
    allocate(plots(n_plots))

    READ_PLOTS: do i = 1, n_plots
      pl => plots(i)

      ! Get pointer to plot XML node
      call get_list_item(node_plot_list, i, node_plot)

      ! Copy data into plots
      if (check_for_node(node_plot, "id")) then
        call get_node_value(node_plot, "id", pl % id)
      else
        message = "Must specify plot id in plots XML file."
        call fatal_error()
      end if

      ! Check to make sure 'id' hasn't been used
      if (plot_dict % has_key(pl % id)) then
        message = "Two or more plots use the same unique ID: " // &
             to_str(pl % id)
        call fatal_error()
      end if

      ! Copy plot type
      temp_str = 'slice'
      if (check_for_node(node_plot, "type")) &
        call get_node_value(node_plot, "type", temp_str)
      call lower_case(temp_str)
      select case (trim(temp_str))
      case ("slice")
        pl % type = PLOT_TYPE_SLICE
      case ("voxel")
        pl % type = PLOT_TYPE_VOXEL
      case default
        message = "Unsupported plot type '" // trim(temp_str) &
             // "' in plot " // trim(to_str(pl % id))
        call fatal_error()
      end select

      ! Set output file path
      filename = "plot"
      if (check_for_node(node_plot, "filename")) &
        call get_node_value(node_plot, "filename", filename)
      select case (pl % type)
      case (PLOT_TYPE_SLICE)
        pl % path_plot = trim(path_input) // trim(to_str(pl % id)) // &
             "_" // trim(filename) // ".ppm"
      case (PLOT_TYPE_VOXEL)
        pl % path_plot = trim(path_input) // trim(to_str(pl % id)) // &
             "_" // trim(filename) // ".voxel"
      end select
      
      ! Copy plot pixel size
      if (pl % type == PLOT_TYPE_SLICE) then
        if (get_arraysize_integer(node_plot, "pixels") == 2) then
          call get_node_array(node_plot, "pixels", pl % pixels(1:2))
        else
          message = "<pixels> must be length 2 in slice plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (get_arraysize_integer(node_plot, "pixels") == 3) then
          call get_node_array(node_plot, "pixels", pl % pixels(1:3))
        else
          message = "<pixels> must be length 3 in voxel plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      end if

      ! Copy plot background color
      if (check_for_node(node_plot, "background")) then
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Background color ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
        if (get_arraysize_integer(node_plot, "background") == 3) then
          call get_node_array(node_plot, "background", pl % not_found % rgb)
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
        temp_str = 'xy'
        if (check_for_node(node_plot, "basis")) &
          call get_node_value(node_plot, "basis", temp_str)
        call lower_case(temp_str)
        select case (trim(temp_str))
        case ("xy")
          pl % basis = PLOT_BASIS_XY
        case ("xz")
          pl % basis = PLOT_BASIS_XZ
        case ("yz")
          pl % basis = PLOT_BASIS_YZ
        case default
          message = "Unsupported plot basis '" // trim(temp_str) & 
               // "' in plot " // trim(to_str(pl % id))
          call fatal_error()
        end select
      end if
      
      ! Copy plotting origin
      if (get_arraysize_double(node_plot, "origin") == 3) then
        call get_node_array(node_plot, "origin", pl % origin)
      else
        message = "Origin must be length 3 " &
             // "in plot " // trim(to_str(pl % id))
        call fatal_error()
      end if

      ! Copy plotting width
      if (pl % type == PLOT_TYPE_SLICE) then
        if (get_arraysize_double(node_plot, "width") == 2) then
          call get_node_array(node_plot, "width", pl % width(1:2))
        else
          message = "<width> must be length 2 in slice plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (get_arraysize_double(node_plot, "width") == 3) then
          call get_node_array(node_plot, "width", pl % width(1:3))
        else
          message = "<width> must be length 3 in voxel plot " // &
                    trim(to_str(pl % id))
          call fatal_error()
        end if
      end if

      ! Copy plot color type and initialize all colors randomly
      temp_str = "cell"
      if (check_for_node(node_plot, "color")) &
        call get_node_value(node_plot, "color", temp_str)
      call lower_case(temp_str)
      select case (trim(temp_str))
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
        message = "Unsupported plot color type '" // trim(temp_str) &
             // "' in plot " // trim(to_str(pl % id))
        call fatal_error()
      end select

      ! Get the number of <col_spec> nodes and get a list of them
      call get_node_list(node_plot, "col_spec", node_col_list)
      n_cols = get_list_size(node_col_list)

      ! Copy user specified colors
      if (n_cols /= 0) then
      
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Color specifications ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
      
        do j = 1, n_cols

          ! Get pointer to color spec XML node
          call get_list_item(node_col_list, j, node_col)

          ! Check and make sure 3 values are specified for RGB
          if (get_arraysize_double(node_col, "rgb") /= 3) then
            message = "Bad RGB " &
                 // "in plot " // trim(to_str(pl % id))
            call fatal_error()          
          end if

          ! Ensure that there is an id for this color specification
          if (check_for_node(node_col, "id")) then
            call get_node_value(node_col, "id", col_id)
          else
            message = "Must specify id for color specification in plot " // &
                      trim(to_str(pl % id))
            call fatal_error()
          end if

          ! Add RGB
          if (pl % color_by == PLOT_COLOR_CELLS) then

            if (cell_dict % has_key(col_id)) then
              col_id = cell_dict % get_key(col_id)
              call get_node_array(node_col, "rgb", pl % colors(col_id) % rgb)
            else
              message = "Could not find cell " // trim(to_str(col_id)) // &
                   " specified in plot " // trim(to_str(pl % id))
              call fatal_error()
            end if

          else if (pl % color_by == PLOT_COLOR_MATS) then

            if (material_dict % has_key(col_id)) then
              col_id = material_dict % get_key(col_id)
              call get_node_array(node_col, "rgb", pl % colors(col_id) % rgb)
            else
              message = "Could not find material " // trim(to_str(col_id)) // &
                   " specified in plot " // trim(to_str(pl % id))
              call fatal_error()
            end if

          end if
        end do
      end if

      ! Deal with masks
      call get_node_list(node_plot, "mask", node_mask_list)
      n_masks = get_list_size(node_mask_list)
      if (n_masks /= 0) then
      
        if (pl % type == PLOT_TYPE_VOXEL) then
          message = "Mask ignored in voxel plot " // & 
                     trim(to_str(pl % id))
          call warning()
        end if
      
        select case(n_masks)
          case default
            message = "Mutliple masks" // &
                 " specified in plot " // trim(to_str(pl % id))
            call fatal_error()
          case (1)

            ! Get pointer to mask
            call get_list_item(node_mask_list, 1, node_mask)

            ! Determine how many components there are and allocate
            n_comp = 0
            n_comp = get_arraysize_integer(node_mask, "components")
            if (n_comp == 0) then
              message = "Missing <components> in mask of plot " // &
                        trim(to_str(pl % id))
              call fatal_error()
            end if
            allocate(iarray(n_comp))
            call get_node_array(node_mask, "components", iarray)
 
            ! First we need to change the user-specified identifiers to indices
            ! in the cell and material arrays
            do j=1, n_comp
              col_id = iarray(j)
            
              if (pl % color_by == PLOT_COLOR_CELLS) then
              
                if (cell_dict % has_key(col_id)) then
                  iarray(j) = cell_dict % get_key(col_id)
                else
                  message = "Could not find cell " // trim(to_str(col_id)) // &
                       " specified in the mask in plot " // trim(to_str(pl % id))
                  call fatal_error()
                end if
              
              else if (pl % color_by == PLOT_COLOR_MATS) then
              
                if (material_dict % has_key(col_id)) then
                  iarray(j) = material_dict % get_key(col_id)
                else
                  message = "Could not find material " // trim(to_str(col_id)) // &
                       " specified in the mask in plot " // trim(to_str(pl % id))
                  call fatal_error()
                end if
                
              end if  
            end do
          
            ! Alter colors based on mask information
            do j=1,size(pl % colors)
              if (.not. any(j .eq. iarray)) then
                if (check_for_node(node_mask, "background")) then
                  call get_node_array(node_mask, "background", pl % colors(j) % rgb)
                else
                  message = "Missing <background> in mask of plot " // &
                            trim(to_str(pl % id))
                  call fatal_error()
                end if
              end if
            end do

            deallocate(iarray)
            
        end select
        
      end if

      ! Add plot to dictionary
      call plot_dict % add_key(pl % id, i)

    end do READ_PLOTS

    ! Close plots XML file
    call close_xmldoc(doc)

  end subroutine read_plots_xml

!===============================================================================
! READ_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the ACE cross sections that may be used.
!===============================================================================

  subroutine read_cross_sections_xml()

    integer :: i           ! loop index
    integer :: filetype    ! default file type
    integer :: recl        ! default record length
    integer :: entries     ! default number of entries
    logical :: file_exists ! does cross_sections.xml exist?
    character(MAX_WORD_LEN)  :: directory ! directory with cross sections
    character(MAX_LINE_LEN)  :: temp_str
    type(XsListing), pointer :: listing => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_ace => null()
    type(NodeList), pointer :: node_ace_list => null()

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

    ! Parse cross_sections.xml file
    call open_xmldoc(doc, path_cross_sections)

    if (check_for_node(doc, "directory")) then
       ! Copy directory information if present
       call get_node_value(doc, "directory", directory)
    else
       ! If no directory is listed in cross_sections.xml, by default select the
       ! directory in which the cross_sections.xml file resides
       i = index(path_cross_sections, "/", BACK=.true.)
       directory = path_cross_sections(1:i)
    end if

    ! determine whether binary/ascii
    temp_str = ''
    if (check_for_node(doc, "filetype")) &
      call get_node_value(doc, "filetype", temp_str)
    if (trim(temp_str) == 'ascii') then
       filetype = ASCII
    elseif (trim(temp_str) == 'binary') then
       filetype = BINARY
    elseif (len_trim(temp_str) == 0) then
       filetype = ASCII
    else
       message = "Unknown filetype in cross_sections.xml: " // trim(temp_str)
       call fatal_error()
    end if

    ! copy default record length and entries for binary files
    if (filetype == BINARY) then
      call get_node_value(doc, "record_length", recl)
      call get_node_value(doc, "entries", entries)
    end if

    ! Get node list of all <ace_table>
    call get_node_list(doc, "ace_table", node_ace_list)
    n_listings = get_list_size(node_ace_list)

    ! Allocate xs_listings array
    if (n_listings == 0) then
       message = "No ACE table listings present in cross_sections.xml file!"
       call fatal_error()
    else
       allocate(xs_listings(n_listings))
    end if

    do i = 1, n_listings
       listing => xs_listings(i)

       ! Get pointer to ace table XML node
       call get_list_item(node_ace_list, i, node_ace)

       ! copy a number of attributes
       call get_node_value(node_ace, "name", listing % name)
       if (check_for_node(node_ace, "alias")) &
         call get_node_value(node_ace, "alias", listing % alias)
       call get_node_value(node_ace, "zaid", listing % zaid)
       call get_node_value(node_ace, "awr", listing % awr)
       if (check_for_node(node_ace, "temperature")) &
         call get_node_value(node_ace, "temperature", listing % kT)
       call get_node_value(node_ace, "location", listing % location)

       ! determine type of cross section
       if (ends_with(listing % name, 'c')) then
          listing % type = ACE_NEUTRON
       elseif (ends_with(listing % name, 't')) then
          listing % type = ACE_THERMAL
       end if

       ! set filetype, record length, and number of entries
       if (check_for_node(node_ace, "filetype")) then
         temp_str = ''
         call get_node_value(node_ace, "filetype", temp_str)
         if (temp_str == 'ascii') then
           listing % filetype = ASCII
         else if (temp_str == 'binary') then
           listing % filetype = BINARY
         end if
       else
         listing % filetype = filetype
       end if

       ! Set record length and entries for binary files
       if (filetype == BINARY) then
         listing % recl     = recl
         listing % entries  = entries
       end if

       ! determine metastable state
       if (.not.check_for_node(node_ace, "metastable")) then
          listing % metastable = .false.
       else
          listing % metastable = .true.
       end if

       ! determine path of cross section table
       if (check_for_node(node_ace, "path")) then
         call get_node_value(node_ace, "path", temp_str)
       else
         message = "Path missing for isotope " // listing % name
         call fatal_error()
       end if

       if (starts_with(temp_str, '/')) then
          listing % path = trim(temp_str)
       else
          if (ends_with(directory,'/')) then
             listing % path = trim(directory) // trim(temp_str)
          else
             listing % path = trim(directory) // '/' // trim(temp_str)
          end if
       end if

       ! create dictionary entry for both name and alias
       call xs_listing_dict % add_key(listing % name, i)
       if (check_for_node(node_ace, "alias")) then
         call xs_listing_dict % add_key(listing % alias, i)
       end if
    end do

    ! Close cross sections XML file
    call close_xmldoc(doc)

  end subroutine read_cross_sections_xml

!===============================================================================
! EXPAND_NATURAL_ELEMENT converts natural elements specified using an <element>
! tag within a material into individual isotopes based on data from the NIST
! Atomic Weights and Isotopic Compositions database.  It was last updated July
! 2010.  In some cases, modifications have been made to work with ENDF/B-VII.1
! where evaluations of particular isotopes don't exist.  The case statements in
! this subroutine were automatically written by the 'write_atomic_weights.py'
! script in the /src/utils directory.
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

    ! The following 'case' code is automatically generated from the script
    ! write_atomic_weights in the src/utils directory.  Running the script
    ! will automatically delete and rewrite any code between the 'Start of...'
    ! and 'End of...' comments.

    ! Start of the natural element cases.
    case ('h')
      call list_names % append('1001.' // xs)
      call list_names % append('1002.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.99988500000000_8)
        call list_density % append(density * 0.00011500000000_8)
      else
        call list_density % append(density * 0.99977095084163_8)
        call list_density % append(density * 0.00022979711535_8)
      end if

    case ('he')
      call list_names % append('2003.' // xs)
      call list_names % append('2004.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00000134000000_8)
        call list_density % append(density * 0.99999866000000_8)
      else
        call list_density % append(density * 0.00000100971300_8)
        call list_density % append(density * 0.99999897333326_8)
      end if

    case ('li')
      call list_names % append('3006.' // xs)
      call list_names % append('3007.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.07590000000000_8)
        call list_density % append(density * 0.92410000000000_8)
      else
        call list_density % append(density * 0.06577551075357_8)
        call list_density % append(density * 0.93408583844619_8)
      end if

    case ('be')
      call list_names % append('4009.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000002219218_8)
      end if

    case ('b')
      call list_names % append('5010.' // xs)
      call list_names % append('5011.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.19900000000000_8)
        call list_density % append(density * 0.80100000000000_8)
      else
        call list_density % append(density * 0.18430991240403_8)
        call list_density % append(density * 0.81569268572750_8)
      end if

    case ('c')
      ! No evaluations split up Carbon into isotopes yet
      call list_names % append('6000.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000000000_8)
      end if

    case ('n')
      call list_names % append('7014.' // xs)
      call list_names % append('7015.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.99636000000000_8)
        call list_density % append(density * 0.00364000000000_8)
      else
        call list_density % append(density * 0.99610206654119_8)
        call list_density % append(density * 0.00389816276421_8)
      end if

    case ('o')
      if (default_expand == JEFF_32) then
        call list_names % append('8016.' // xs)
        call list_names % append('8017.' // xs)
        call list_names % append('8018.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.99757000000000_8)
          call list_density % append(density * 0.00038000000000_8)
          call list_density % append(density * 0.00205000000000_8)
        else
          call list_density % append(density * 0.99729033445220_8)
          call list_density % append(density * 0.00040374451829_8)
          call list_density % append(density * 0.00230622898671_8)
        end if
      elseif (default_expand >= JENDL_32 .and. default_expand <= JENDL_40) then
        call list_names % append('8016.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('8016.' // xs)
        call list_names % append('8017.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.99962000000000_8)
          call list_density % append(density * 0.00038000000000_8)
        else
          call list_density % append(density * 0.99959656343891_8)
          call list_density % append(density * 0.00040374451829_8)
        end if
      end if

    case ('f')
      call list_names % append('9019.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000105272_8)
      end if

    case ('ne')
      call list_names % append('10020.' // xs)
      call list_names % append('10021.' // xs)
      call list_names % append('10022.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.90480000000000_8)
        call list_density % append(density * 0.00270000000000_8)
        call list_density % append(density * 0.09250000000000_8)
      else
        call list_density % append(density * 0.89640380534408_8)
        call list_density % append(density * 0.00280893105626_8)
        call list_density % append(density * 0.10080442836340_8)
      end if

    case ('na')
      call list_names % append('11023.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000003915_8)
      end if

    case ('mg')
      call list_names % append('12024.' // xs)
      call list_names % append('12025.' // xs)
      call list_names % append('12026.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.78990000000000_8)
        call list_density % append(density * 0.10000000000000_8)
        call list_density % append(density * 0.11010000000000_8)
      else
        call list_density % append(density * 0.77950151980374_8)
        call list_density % append(density * 0.10280122164164_8)
        call list_density % append(density * 0.11769938208117_8)
      end if

    case ('al')
      call list_names % append('13027.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000111187_8)
      end if

    case ('si')
      call list_names % append('14028.' // xs)
      call list_names % append('14029.' // xs)
      call list_names % append('14030.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.92223000000000_8)
        call list_density % append(density * 0.04685000000000_8)
        call list_density % append(density * 0.03092000000000_8)
      else
        call list_density % append(density * 0.91866482548174_8)
        call list_density % append(density * 0.04833628657831_8)
        call list_density % append(density * 0.03299884188127_8)
      end if

    case ('p')
      call list_names % append('15031.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999998805441_8)
      end if

    case ('s')
      call list_names % append('16032.' // xs)
      call list_names % append('16033.' // xs)
      call list_names % append('16034.' // xs)
      call list_names % append('16036.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.94990000000000_8)
        call list_density % append(density * 0.00750000000000_8)
        call list_density % append(density * 0.04250000000000_8)
        call list_density % append(density * 0.00010000000000_8)
      else
        call list_density % append(density * 0.94714705263995_8)
        call list_density % append(density * 0.00771202060502_8)
        call list_density % append(density * 0.04502212204117_8)
        call list_density % append(density * 0.00011216928352_8)
      end if

    case ('cl')
      call list_names % append('17035.' // xs)
      call list_names % append('17037.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.75760000000000_8)
        call list_density % append(density * 0.24240000000000_8)
      else
        call list_density % append(density * 0.74725418978275_8)
        call list_density % append(density * 0.25274404952517_8)
      end if

    case ('ar')
      call list_names % append('18036.' // xs)
      call list_names % append('18038.' // xs)
      call list_names % append('18040.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00336500000000_8)
        call list_density % append(density * 0.00063200000000_8)
        call list_density % append(density * 0.99600300000000_8)
      else
        call list_density % append(density * 0.00302970835290_8)
        call list_density % append(density * 0.00060059194144_8)
        call list_density % append(density * 0.99636160701811_8)
      end if

    case ('k')
      call list_names % append('19039.' // xs)
      call list_names % append('19040.' // xs)
      call list_names % append('19041.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.93258100000000_8)
        call list_density % append(density * 0.00011700000000_8)
        call list_density % append(density * 0.06730200000000_8)
      else
        call list_density % append(density * 0.92937065139254_8)
        call list_density % append(density * 0.00011959056589_8)
        call list_density % append(density * 0.07050978680146_8)
      end if

    case ('ca')
      call list_names % append('20040.' // xs)
      call list_names % append('20042.' // xs)
      call list_names % append('20043.' // xs)
      call list_names % append('20044.' // xs)
      call list_names % append('20046.' // xs)
      call list_names % append('20048.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.96941000000000_8)
        call list_density % append(density * 0.00647000000000_8)
        call list_density % append(density * 0.00135000000000_8)
        call list_density % append(density * 0.02086000000000_8)
        call list_density % append(density * 0.00004000000000_8)
        call list_density % append(density * 0.00187000000000_8)
      else
        call list_density % append(density * 0.96661847701786_8)
        call list_density % append(density * 0.00677359794712_8)
        call list_density % append(density * 0.00144703665128_8)
        call list_density % append(density * 0.02287817132462_8)
        call list_density % append(density * 0.00004586425730_8)
        call list_density % append(density * 0.00223741799940_8)
      end if

    case ('sc')
      call list_names % append('21045.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999999777560_8)
      end if

    case ('ti')
      call list_names % append('22046.' // xs)
      call list_names % append('22047.' // xs)
      call list_names % append('22048.' // xs)
      call list_names % append('22049.' // xs)
      call list_names % append('22050.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.08250000000000_8)
        call list_density % append(density * 0.07440000000000_8)
        call list_density % append(density * 0.73720000000000_8)
        call list_density % append(density * 0.05410000000000_8)
        call list_density % append(density * 0.05180000000000_8)
      else
        call list_density % append(density * 0.07920053705058_8)
        call list_density % append(density * 0.07297744113147_8)
        call list_density % append(density * 0.73844665452943_8)
        call list_density % append(density * 0.05532161545532_8)
        call list_density % append(density * 0.05404851325882_8)
      end if

    case ('v')
      if (default_expand == ENDF_BVII0 .or. default_expand == JEFF_311 &
           .or. default_expand == JEFF_32 .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_33)) then
        call list_names % append('23000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('23050.' // xs)
        call list_names % append('23051.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00250000000000_8)
          call list_density % append(density * 0.99750000000000_8)
        else
          call list_density % append(density * 0.00245120179520_8)
          call list_density % append(density * 0.99754816016902_8)
        end if
      end if

    case ('cr')
      call list_names % append('24050.' // xs)
      call list_names % append('24052.' // xs)
      call list_names % append('24053.' // xs)
      call list_names % append('24054.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.04345000000000_8)
        call list_density % append(density * 0.83789000000000_8)
        call list_density % append(density * 0.09501000000000_8)
        call list_density % append(density * 0.02365000000000_8)
      else
        call list_density % append(density * 0.04173689219941_8)
        call list_density % append(density * 0.83699415589198_8)
        call list_density % append(density * 0.09673593018503_8)
        call list_density % append(density * 0.02453365774472_8)
      end if

    case ('mn')
      call list_names % append('25055.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000182023_8)
      end if

    case ('fe')
      call list_names % append('26054.' // xs)
      call list_names % append('26056.' // xs)
      call list_names % append('26057.' // xs)
      call list_names % append('26058.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.05845000000000_8)
        call list_density % append(density * 0.91754000000000_8)
        call list_density % append(density * 0.02119000000000_8)
        call list_density % append(density * 0.00282000000000_8)
      else
        call list_density % append(density * 0.05645572985451_8)
        call list_density % append(density * 0.91901768383472_8)
        call list_density % append(density * 0.02160374248115_8)
        call list_density % append(density * 0.00292545146731_8)
      end if

    case ('co')
      call list_names % append('27059.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000000000_8)
      end if

    case ('ni')
      call list_names % append('28058.' // xs)
      call list_names % append('28060.' // xs)
      call list_names % append('28061.' // xs)
      call list_names % append('28062.' // xs)
      call list_names % append('28064.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.68076900000000_8)
        call list_density % append(density * 0.26223100000000_8)
        call list_density % append(density * 0.01139900000000_8)
        call list_density % append(density * 0.03634500000000_8)
        call list_density % append(density * 0.00925600000000_8)
      else
        call list_density % append(density * 0.67197649907298_8)
        call list_density % append(density * 0.26775940818658_8)
        call list_density % append(density * 0.01183358107290_8)
        call list_density % append(density * 0.03834819081293_8)
        call list_density % append(density * 0.01008149559058_8)
      end if

    case ('cu')
      call list_names % append('29063.' // xs)
      call list_names % append('29065.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.69150000000000_8)
        call list_density % append(density * 0.30850000000000_8)
      else
        call list_density % append(density * 0.68479238144415_8)
        call list_density % append(density * 0.31520824380370_8)
      end if

    case ('zn')
      if (default_expand == ENDF_BVII0 .or. default_expand == &
           JEFF_311 .or. default_expand == JEFF_312) then
        call list_names % append('30000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('30064.' // xs)
        call list_names % append('30066.' // xs)
        call list_names % append('30067.' // xs)
        call list_names % append('30068.' // xs)
        call list_names % append('30070.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.48268000000000_8)
          call list_density % append(density * 0.27975000000000_8)
          call list_density % append(density * 0.04102000000000_8)
          call list_density % append(density * 0.19024000000000_8)
          call list_density % append(density * 0.00631000000000_8)
        else
          call list_density % append(density * 0.47196877266895_8)
          call list_density % append(density * 0.28208638488299_8)
          call list_density % append(density * 0.04199068158223_8)
          call list_density % append(density * 0.19764488162447_8)
          call list_density % append(density * 0.00674868101534_8)
        end if
      end if

    case ('ga')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call list_names % append('31000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('31069.' // xs)
        call list_names % append('31071.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.60108000000000_8)
          call list_density % append(density * 0.39892000000000_8)
        else
          call list_density % append(density * 0.59420540968530_8)
          call list_density % append(density * 0.40579553149744_8)
        end if
      end if

    case ('ge')
      call list_names % append('32070.' // xs)
      call list_names % append('32072.' // xs)
      call list_names % append('32073.' // xs)
      call list_names % append('32074.' // xs)
      call list_names % append('32076.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.20380000000000_8)
        call list_density % append(density * 0.27310000000000_8)
        call list_density % append(density * 0.07760000000000_8)
        call list_density % append(density * 0.36720000000000_8)
        call list_density % append(density * 0.07830000000000_8)
      else
        call list_density % append(density * 0.19618063904350_8)
        call list_density % append(density * 0.27040086592759_8)
        call list_density % append(density * 0.07790281402313_8)
        call list_density % append(density * 0.37367643843833_8)
        call list_density % append(density * 0.08183708457572_8)
      end if

    case ('as')
      call list_names % append('33075.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999995328450_8)
      end if

    case ('se')
      call list_names % append('34074.' // xs)
      call list_names % append('34076.' // xs)
      call list_names % append('34077.' // xs)
      call list_names % append('34078.' // xs)
      call list_names % append('34080.' // xs)
      call list_names % append('34082.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00890000000000_8)
        call list_density % append(density * 0.09370000000000_8)
        call list_density % append(density * 0.07630000000000_8)
        call list_density % append(density * 0.23770000000000_8)
        call list_density % append(density * 0.49610000000000_8)
        call list_density % append(density * 0.08730000000000_8)
      else
        call list_density % append(density * 0.00833219402178_8)
        call list_density % append(density * 0.09009156933029_8)
        call list_density % append(density * 0.07432864030142_8)
        call list_density % append(density * 0.23456109894972_8)
        call list_density % append(density * 0.50210975452039_8)
        call list_density % append(density * 0.09056899515729_8)
      end if

    case ('br')
      call list_names % append('35079.' // xs)
      call list_names % append('35081.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.50690000000000_8)
        call list_density % append(density * 0.49310000000000_8)
      else
        call list_density % append(density * 0.50064708995782_8)
        call list_density % append(density * 0.49934700258886_8)
      end if

    case ('kr')
      call list_names % append('36078.' // xs)
      call list_names % append('36080.' // xs)
      call list_names % append('36082.' // xs)
      call list_names % append('36083.' // xs)
      call list_names % append('36084.' // xs)
      call list_names % append('36086.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00355000000000_8)
        call list_density % append(density * 0.02286000000000_8)
        call list_density % append(density * 0.11593000000000_8)
        call list_density % append(density * 0.11500000000000_8)
        call list_density % append(density * 0.56987000000000_8)
        call list_density % append(density * 0.17279000000000_8)
      else
        call list_density % append(density * 0.00330100115802_8)
        call list_density % append(density * 0.02180109816392_8)
        call list_density % append(density * 0.11332287350233_8)
        call list_density % append(density * 0.11378703119406_8)
        call list_density % append(density * 0.57064190665756_8)
        call list_density % append(density * 0.17714616611419_8)
      end if

    case ('rb')
      call list_names % append('37085.' // xs)
      call list_names % append('37087.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.72170000000000_8)
        call list_density % append(density * 0.27830000000000_8)
      else
        call list_density % append(density * 0.71700498496410_8)
        call list_density % append(density * 0.28299341904980_8)
      end if

    case ('sr')
      call list_names % append('38084.' // xs)
      call list_names % append('38086.' // xs)
      call list_names % append('38087.' // xs)
      call list_names % append('38088.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00560000000000_8)
        call list_density % append(density * 0.09860000000000_8)
        call list_density % append(density * 0.07000000000000_8)
        call list_density % append(density * 0.82580000000000_8)
      else
        call list_density % append(density * 0.00536310408583_8)
        call list_density % append(density * 0.09667488080027_8)
        call list_density % append(density * 0.06943188081488_8)
        call list_density % append(density * 0.82849183373864_8)
      end if

    case ('y')
      call list_names % append('39089.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999998087865_8)
      end if

    case ('zr')
      call list_names % append('40090.' // xs)
      call list_names % append('40091.' // xs)
      call list_names % append('40092.' // xs)
      call list_names % append('40094.' // xs)
      call list_names % append('40096.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.51450000000000_8)
        call list_density % append(density * 0.11220000000000_8)
        call list_density % append(density * 0.17150000000000_8)
        call list_density % append(density * 0.17380000000000_8)
        call list_density % append(density * 0.02800000000000_8)
      else
        call list_density % append(density * 0.50705922140884_8)
        call list_density % append(density * 0.11180844359774_8)
        call list_density % append(density * 0.17278034834254_8)
        call list_density % append(density * 0.17891034795405_8)
        call list_density % append(density * 0.02943777575200_8)
      end if

    case ('nb')
      call list_names % append('41093.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999997954931_8)
      end if

    case ('mo')
      call list_names % append('42092.' // xs)
      call list_names % append('42094.' // xs)
      call list_names % append('42095.' // xs)
      call list_names % append('42096.' // xs)
      call list_names % append('42097.' // xs)
      call list_names % append('42098.' // xs)
      call list_names % append('42100.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.14770000000000_8)
        call list_density % append(density * 0.09230000000000_8)
        call list_density % append(density * 0.15900000000000_8)
        call list_density % append(density * 0.16680000000000_8)
        call list_density % append(density * 0.09560000000000_8)
        call list_density % append(density * 0.24190000000000_8)
        call list_density % append(density * 0.09670000000000_8)
      else
        call list_density % append(density * 0.14146140042414_8)
        call list_density % append(density * 0.09032346446530_8)
        call list_density % append(density * 0.15725332319612_8)
        call list_density % append(density * 0.16670384056482_8)
        call list_density % append(density * 0.09654247244060_8)
        call list_density % append(density * 0.24680406673176_8)
        call list_density % append(density * 0.10067791815236_8)
      end if

    case ('ru')
      call list_names % append('44096.' // xs)
      call list_names % append('44098.' // xs)
      call list_names % append('44099.' // xs)
      call list_names % append('44100.' // xs)
      call list_names % append('44101.' // xs)
      call list_names % append('44102.' // xs)
      call list_names % append('44104.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.05540000000000_8)
        call list_density % append(density * 0.01870000000000_8)
        call list_density % append(density * 0.12760000000000_8)
        call list_density % append(density * 0.12600000000000_8)
        call list_density % append(density * 0.17060000000000_8)
        call list_density % append(density * 0.31550000000000_8)
        call list_density % append(density * 0.18620000000000_8)
      else
        call list_density % append(density * 0.05257030700702_8)
        call list_density % append(density * 0.01811446390521_8)
        call list_density % append(density * 0.12486789210132_8)
        call list_density % append(density * 0.12454666723063_8)
        call list_density % append(density * 0.17032247260572_8)
        call list_density % append(density * 0.31810450385030_8)
        call list_density % append(density * 0.19142368283962_8)
      end if

    case ('rh')
      call list_names % append('45103.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000003887061_8)
      end if

    case ('pd')
      call list_names % append('46102.' // xs)
      call list_names % append('46104.' // xs)
      call list_names % append('46105.' // xs)
      call list_names % append('46106.' // xs)
      call list_names % append('46108.' // xs)
      call list_names % append('46110.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.01020000000000_8)
        call list_density % append(density * 0.11140000000000_8)
        call list_density % append(density * 0.22330000000000_8)
        call list_density % append(density * 0.27330000000000_8)
        call list_density % append(density * 0.26460000000000_8)
        call list_density % append(density * 0.11720000000000_8)
      else
        call list_density % append(density * 0.00976731076677_8)
        call list_density % append(density * 0.10876629966548_8)
        call list_density % append(density * 0.22012126931498_8)
        call list_density % append(density * 0.27197352681639_8)
        call list_density % append(density * 0.26828951158805_8)
        call list_density % append(density * 0.12103818766773_8)
      end if

    case ('ag')
      call list_names % append('47107.' // xs)
      call list_names % append('47109.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.51839000000000_8)
        call list_density % append(density * 0.48161000000000_8)
      else
        call list_density % append(density * 0.51376154634851_8)
        call list_density % append(density * 0.48623799795232_8)
      end if

    case ('cd')
      call list_names % append('48106.' // xs)
      call list_names % append('48108.' // xs)
      call list_names % append('48110.' // xs)
      call list_names % append('48111.' // xs)
      call list_names % append('48112.' // xs)
      call list_names % append('48113.' // xs)
      call list_names % append('48114.' // xs)
      call list_names % append('48116.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.01250000000000_8)
        call list_density % append(density * 0.00890000000000_8)
        call list_density % append(density * 0.12490000000000_8)
        call list_density % append(density * 0.12800000000000_8)
        call list_density % append(density * 0.24130000000000_8)
        call list_density % append(density * 0.12220000000000_8)
        call list_density % append(density * 0.28730000000000_8)
        call list_density % append(density * 0.07490000000000_8)
      else
        call list_density % append(density * 0.01177670101236_8)
        call list_density % append(density * 0.00854317849321_8)
        call list_density % append(density * 0.12211336045663_8)
        call list_density % append(density * 0.12628421414986_8)
        call list_density % append(density * 0.24020901386110_8)
        call list_density % append(density * 0.12273636821788_8)
        call list_density % append(density * 0.29111416940557_8)
        call list_density % append(density * 0.07722790673866_8)
      end if

    case ('in')
      call list_names % append('49113.' // xs)
      call list_names % append('49115.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.04290000000000_8)
        call list_density % append(density * 0.95710000000000_8)
      else
        call list_density % append(density * 0.04218488467139_8)
        call list_density % append(density * 0.95781586191886_8)
      end if

    case ('sn')
      call list_names % append('50112.' // xs)
      call list_names % append('50114.' // xs)
      call list_names % append('50115.' // xs)
      call list_names % append('50116.' // xs)
      call list_names % append('50117.' // xs)
      call list_names % append('50118.' // xs)
      call list_names % append('50119.' // xs)
      call list_names % append('50120.' // xs)
      call list_names % append('50122.' // xs)
      call list_names % append('50124.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00970000000000_8)
        call list_density % append(density * 0.00660000000000_8)
        call list_density % append(density * 0.00340000000000_8)
        call list_density % append(density * 0.14540000000000_8)
        call list_density % append(density * 0.07680000000000_8)
        call list_density % append(density * 0.24220000000000_8)
        call list_density % append(density * 0.08590000000000_8)
        call list_density % append(density * 0.32580000000000_8)
        call list_density % append(density * 0.04630000000000_8)
        call list_density % append(density * 0.05790000000000_8)
      else
        call list_density % append(density * 0.00914393677533_8)
        call list_density % append(density * 0.00633272968916_8)
        call list_density % append(density * 0.00329097264594_8)
        call list_density % append(density * 0.14196034994019_8)
        call list_density % append(density * 0.07563092168815_8)
        call list_density % append(density * 0.24055065492882_8)
        call list_density % append(density * 0.08603988002022_8)
        call list_density % append(density * 0.32907198242153_8)
        call list_density % append(density * 0.04754552460366_8)
        call list_density % append(density * 0.06043395972378_8)
      end if

    case ('sb')
      call list_names % append('51121.' // xs)
      call list_names % append('51123.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.57210000000000_8)
        call list_density % append(density * 0.42790000000000_8)
      else
        call list_density % append(density * 0.56807714324877_8)
        call list_density % append(density * 0.43192110028417_8)
      end if

    case ('te')
      call list_names % append('52120.' // xs)
      call list_names % append('52122.' // xs)
      call list_names % append('52123.' // xs)
      call list_names % append('52124.' // xs)
      call list_names % append('52125.' // xs)
      call list_names % append('52126.' // xs)
      call list_names % append('52128.' // xs)
      call list_names % append('52130.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00090000000000_8)
        call list_density % append(density * 0.02550000000000_8)
        call list_density % append(density * 0.00890000000000_8)
        call list_density % append(density * 0.04740000000000_8)
        call list_density % append(density * 0.07070000000000_8)
        call list_density % append(density * 0.18840000000000_8)
        call list_density % append(density * 0.31740000000000_8)
        call list_density % append(density * 0.34080000000000_8)
      else
        call list_density % append(density * 0.00084571800940_8)
        call list_density % append(density * 0.02436150171983_8)
        call list_density % append(density * 0.00857247651254_8)
        call list_density % append(density * 0.04602659536411_8)
        call list_density % append(density * 0.06920645180635_8)
        call list_density % append(density * 0.18589485834075_8)
        call list_density % append(density * 0.31815734003088_8)
        call list_density % append(density * 0.34695957112476_8)
      end if

    case ('i')
      call list_names % append('53127.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000002363983_8)
      end if

    case ('xe')
      call list_names % append('54124.' // xs)
      call list_names % append('54126.' // xs)
      call list_names % append('54128.' // xs)
      call list_names % append('54129.' // xs)
      call list_names % append('54130.' // xs)
      call list_names % append('54131.' // xs)
      call list_names % append('54132.' // xs)
      call list_names % append('54134.' // xs)
      call list_names % append('54136.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00095200000000_8)
        call list_density % append(density * 0.00089000000000_8)
        call list_density % append(density * 0.01910200000000_8)
        call list_density % append(density * 0.26400600000000_8)
        call list_density % append(density * 0.04071000000000_8)
        call list_density % append(density * 0.21232400000000_8)
        call list_density % append(density * 0.26908600000000_8)
        call list_density % append(density * 0.10435700000000_8)
        call list_density % append(density * 0.08857300000000_8)
      else
        call list_density % append(density * 0.00089843639902_8)
        call list_density % append(density * 0.00085347127311_8)
        call list_density % append(density * 0.01860886151503_8)
        call list_density % append(density * 0.25920372898994_8)
        call list_density % append(density * 0.04027916043262_8)
        call list_density % append(density * 0.21169666863807_8)
        call list_density % append(density * 0.27033856373684_8)
        call list_density % append(density * 0.10643343707461_8)
        call list_density % append(density * 0.09168584851048_8)
      end if

    case ('cs')
      call list_names % append('55133.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000000024830_8)
      end if

    case ('ba')
      call list_names % append('56130.' // xs)
      call list_names % append('56132.' // xs)
      call list_names % append('56134.' // xs)
      call list_names % append('56135.' // xs)
      call list_names % append('56136.' // xs)
      call list_names % append('56137.' // xs)
      call list_names % append('56138.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00106000000000_8)
        call list_density % append(density * 0.00101000000000_8)
        call list_density % append(density * 0.02417000000000_8)
        call list_density % append(density * 0.06592000000000_8)
        call list_density % append(density * 0.07854000000000_8)
        call list_density % append(density * 0.11232000000000_8)
        call list_density % append(density * 0.71698000000000_8)
      else
        call list_density % append(density * 0.00100272124235_8)
        call list_density % append(density * 0.00097012322350_8)
        call list_density % append(density * 0.02356763031325_8)
        call list_density % append(density * 0.06475771692757_8)
        call list_density % append(density * 0.07772648780783_8)
        call list_density % append(density * 0.11197552217385_8)
        call list_density % append(density * 0.71999901066401_8)
      end if

    case ('la')
      call list_names % append('57138.' // xs)
      call list_names % append('57139.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00090000000000_8)
        call list_density % append(density * 0.99910000000000_8)
      else
        call list_density % append(density * 0.00089353141241_8)
        call list_density % append(density * 0.99910635327774_8)
      end if

    case ('ce')
      call list_names % append('58136.' // xs)
      call list_names % append('58138.' // xs)
      call list_names % append('58140.' // xs)
      call list_names % append('58142.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00185000000000_8)
        call list_density % append(density * 0.00251000000000_8)
        call list_density % append(density * 0.88450000000000_8)
        call list_density % append(density * 0.11114000000000_8)
      else
        call list_density % append(density * 0.00179442938851_8)
        call list_density % append(density * 0.00247041049852_8)
        call list_density % append(density * 0.88317080511969_8)
        call list_density % append(density * 0.11256240099746_8)
      end if

    case ('pr')
      call list_names % append('59141.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000001987117_8)
      end if

    case ('nd')
      call list_names % append('60142.' // xs)
      call list_names % append('60143.' // xs)
      call list_names % append('60144.' // xs)
      call list_names % append('60145.' // xs)
      call list_names % append('60146.' // xs)
      call list_names % append('60148.' // xs)
      call list_names % append('60150.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.27200000000000_8)
        call list_density % append(density * 0.12200000000000_8)
        call list_density % append(density * 0.23800000000000_8)
        call list_density % append(density * 0.08300000000000_8)
        call list_density % append(density * 0.17200000000000_8)
        call list_density % append(density * 0.05700000000000_8)
        call list_density % append(density * 0.05600000000000_8)
      else
        call list_density % append(density * 0.26759820813355_8)
        call list_density % append(density * 0.12087323625990_8)
        call list_density % append(density * 0.23745234243424_8)
        call list_density % append(density * 0.08338586270850_8)
        call list_density % append(density * 0.17399270744166_8)
        call list_density % append(density * 0.05845220463527_8)
        call list_density % append(density * 0.05820475240221_8)
      end if

    case ('sm')
      call list_names % append('62144.' // xs)
      call list_names % append('62147.' // xs)
      call list_names % append('62148.' // xs)
      call list_names % append('62149.' // xs)
      call list_names % append('62150.' // xs)
      call list_names % append('62152.' // xs)
      call list_names % append('62154.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.03070000000000_8)
        call list_density % append(density * 0.14990000000000_8)
        call list_density % append(density * 0.11240000000000_8)
        call list_density % append(density * 0.13820000000000_8)
        call list_density % append(density * 0.07380000000000_8)
        call list_density % append(density * 0.26750000000000_8)
        call list_density % append(density * 0.22750000000000_8)
      else
        call list_density % append(density * 0.02938346880354_8)
        call list_density % append(density * 0.14646543758453_8)
        call list_density % append(density * 0.11057213402155_8)
        call list_density % append(density * 0.13687386888494_8)
        call list_density % append(density * 0.07358270106345_8)
        call list_density % append(density * 0.27027486310854_8)
        call list_density % append(density * 0.23288974870810_8)
      end if

    case ('eu')
      call list_names % append('63151.' // xs)
      call list_names % append('63153.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.47810000000000_8)
        call list_density % append(density * 0.52190000000000_8)
      else
        call list_density % append(density * 0.47481495867850_8)
        call list_density % append(density * 0.52518747922909_8)
      end if

    case ('gd')
      call list_names % append('64152.' // xs)
      call list_names % append('64154.' // xs)
      call list_names % append('64155.' // xs)
      call list_names % append('64156.' // xs)
      call list_names % append('64157.' // xs)
      call list_names % append('64158.' // xs)
      call list_names % append('64160.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00200000000000_8)
        call list_density % append(density * 0.02180000000000_8)
        call list_density % append(density * 0.14800000000000_8)
        call list_density % append(density * 0.20470000000000_8)
        call list_density % append(density * 0.15650000000000_8)
        call list_density % append(density * 0.24840000000000_8)
        call list_density % append(density * 0.21860000000000_8)
      else
        call list_density % append(density * 0.00193220719873_8)
        call list_density % append(density * 0.02133847294169_8)
        call list_density % append(density * 0.14580952658824_8)
        call list_density % append(density * 0.20297143730804_8)
        call list_density % append(density * 0.15617551513927_8)
        call list_density % append(density * 0.24946484838639_8)
        call list_density % append(density * 0.22232148824331_8)
      end if

    case ('tb')
      call list_names % append('65159.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999997986476_8)
      end if

    case ('dy')
      call list_names % append('66156.' // xs)
      call list_names % append('66158.' // xs)
      call list_names % append('66160.' // xs)
      call list_names % append('66161.' // xs)
      call list_names % append('66162.' // xs)
      call list_names % append('66163.' // xs)
      call list_names % append('66164.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00056000000000_8)
        call list_density % append(density * 0.00095000000000_8)
        call list_density % append(density * 0.02329000000000_8)
        call list_density % append(density * 0.18889000000000_8)
        call list_density % append(density * 0.25475000000000_8)
        call list_density % append(density * 0.24896000000000_8)
        call list_density % append(density * 0.28260000000000_8)
      else
        call list_density % append(density * 0.00053733906757_8)
        call list_density % append(density * 0.00092325039108_8)
        call list_density % append(density * 0.02292097138323_8)
        call list_density % append(density * 0.18706146738416_8)
        call list_density % append(density * 0.25385139626092_8)
        call list_density % append(density * 0.24961684258186_8)
        call list_density % append(density * 0.28508544491372_8)
      end if

    case ('ho')
      call list_names % append('67165.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000001273265_8)
      end if

    case ('er')
      call list_names % append('68162.' // xs)
      call list_names % append('68164.' // xs)
      call list_names % append('68166.' // xs)
      call list_names % append('68167.' // xs)
      call list_names % append('68168.' // xs)
      call list_names % append('68170.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00139000000000_8)
        call list_density % append(density * 0.01601000000000_8)
        call list_density % append(density * 0.33503000000000_8)
        call list_density % append(density * 0.22869000000000_8)
        call list_density % append(density * 0.26978000000000_8)
        call list_density % append(density * 0.14910000000000_8)
      else
        call list_density % append(density * 0.00134570337871_8)
        call list_density % append(density * 0.01569127217071_8)
        call list_density % append(density * 0.33236851886770_8)
        call list_density % append(density * 0.22824296511911_8)
        call list_density % append(density * 0.27086611083742_8)
        call list_density % append(density * 0.15148588552562_8)
      end if

    case ('tm')
      call list_names % append('69169.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000001953423_8)
      end if

    case ('yb')
      call list_names % append('70168.' // xs)
      call list_names % append('70170.' // xs)
      call list_names % append('70171.' // xs)
      call list_names % append('70172.' // xs)
      call list_names % append('70173.' // xs)
      call list_names % append('70174.' // xs)
      call list_names % append('70176.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00130000000000_8)
        call list_density % append(density * 0.03040000000000_8)
        call list_density % append(density * 0.14280000000000_8)
        call list_density % append(density * 0.21830000000000_8)
        call list_density % append(density * 0.16130000000000_8)
        call list_density % append(density * 0.31830000000000_8)
        call list_density % append(density * 0.12760000000000_8)
      else
        call list_density % append(density * 0.00126153724329_8)
        call list_density % append(density * 0.02985205056641_8)
        call list_density % append(density * 0.14105254616617_8)
        call list_density % append(density * 0.21689017348024_8)
        call list_density % append(density * 0.16119207531776_8)
        call list_density % append(density * 0.31992753595080_8)
        call list_density % append(density * 0.12972986552706_8)
      end if

    case ('lu')
      call list_names % append('71175.' // xs)
      call list_names % append('71176.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.97410000000000_8)
        call list_density % append(density * 0.02590000000000_8)
      else
        call list_density % append(density * 0.97395509211107_8)
        call list_density % append(density * 0.02604445857826_8)
      end if

    case ('hf')
      call list_names % append('72174.' // xs)
      call list_names % append('72176.' // xs)
      call list_names % append('72177.' // xs)
      call list_names % append('72178.' // xs)
      call list_names % append('72179.' // xs)
      call list_names % append('72180.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00160000000000_8)
        call list_density % append(density * 0.05260000000000_8)
        call list_density % append(density * 0.18600000000000_8)
        call list_density % append(density * 0.27280000000000_8)
        call list_density % append(density * 0.13620000000000_8)
        call list_density % append(density * 0.35080000000000_8)
      else
        call list_density % append(density * 0.00155921381366_8)
        call list_density % append(density * 0.05184894443588_8)
        call list_density % append(density * 0.18438813967281_8)
        call list_density % append(density * 0.27196504584369_8)
        call list_density % append(density * 0.13654781866110_8)
        call list_density % append(density * 0.35366266872094_8)
      end if

    case ('ta')
      if (default_expand == ENDF_BVII0 .or. &
           (default_expand >= JEFF_311 .and. default_expand <= JEFF_312) .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_40)) then
        call list_names % append('73181.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('73180.' // xs)
        call list_names % append('73181.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00012000000000_8)
          call list_density % append(density * 0.99988000000000_8)
        else
          call list_density % append(density * 0.00011933655026_8)
          call list_density % append(density * 0.99988063988649_8)
        end if
      end if

    case ('w')
      if (default_expand == ENDF_BVII0 .or. default_expand == JEFF_311 &
           .or. default_expand == JEFF_312 .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_33)) then
        ! Combine W-180 with W-182
        call list_names % append('74180.' // xs)
        call list_names % append('74182.' // xs)
        call list_names % append('74183.' // xs)
        call list_names % append('74184.' // xs)
        call list_names % append('74186.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00120000000000_8)
          call list_density % append(density * 0.26500000000000_8)
          call list_density % append(density * 0.14310000000000_8)
          call list_density % append(density * 0.30640000000000_8)
          call list_density % append(density * 0.28430000000000_8)
        else
          call list_density % append(density * 0.00117458684073_8)
          call list_density % append(density * 0.26227303151110_8)
          call list_density % append(density * 0.14240740269419_8)
          call list_density % append(density * 0.30658488533333_8)
          call list_density % append(density * 0.28756976563115_8)
        end if
      else
        call list_names % append('74180.' // xs)
        call list_names % append('74182.' // xs)
        call list_names % append('74183.' // xs)
        call list_names % append('74184.' // xs)
        call list_names % append('74186.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00120000000000_8)
          call list_density % append(density * 0.26500000000000_8)
          call list_density % append(density * 0.14310000000000_8)
          call list_density % append(density * 0.30640000000000_8)
          call list_density % append(density * 0.28430000000000_8)
        else
          call list_density % append(density * 0.00117458684073_8)
          call list_density % append(density * 0.26227303151110_8)
          call list_density % append(density * 0.14240740269419_8)
          call list_density % append(density * 0.30658488533333_8)
          call list_density % append(density * 0.28756976563115_8)
        end if
      end if

    case ('re')
      call list_names % append('75185.' // xs)
      call list_names % append('75187.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.37400000000000_8)
        call list_density % append(density * 0.62600000000000_8)
      else
        call list_density % append(density * 0.37148122879376_8)
        call list_density % append(density * 0.62851719559737_8)
      end if

    case ('os')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call list_names % append('76000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('76184.' // xs)
        call list_names % append('76186.' // xs)
        call list_names % append('76187.' // xs)
        call list_names % append('76188.' // xs)
        call list_names % append('76189.' // xs)
        call list_names % append('76190.' // xs)
        call list_names % append('76192.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00020000000000_8)
          call list_density % append(density * 0.01590000000000_8)
          call list_density % append(density * 0.01960000000000_8)
          call list_density % append(density * 0.13240000000000_8)
          call list_density % append(density * 0.16150000000000_8)
          call list_density % append(density * 0.26260000000000_8)
          call list_density % append(density * 0.40780000000000_8)
        else
          call list_density % append(density * 0.00019340008316_8)
          call list_density % append(density * 0.01554258543542_8)
          call list_density % append(density * 0.01926264369342_8)
          call list_density % append(density * 0.13081718434358_8)
          call list_density % append(density * 0.16042023246202_8)
          call list_density % append(density * 0.26222513894864_8)
          call list_density % append(density * 0.41151181112054_8)
        end if
      end if

    case ('ir')
      call list_names % append('77191.' // xs)
      call list_names % append('77193.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.37300000000000_8)
        call list_density % append(density * 0.62700000000000_8)
      else
        call list_density % append(density * 0.37056192512629_8)
        call list_density % append(density * 0.62943316591561_8)
      end if

    case ('pt')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call list_names % append('78000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('78190.' // xs)
        call list_names % append('78192.' // xs)
        call list_names % append('78194.' // xs)
        call list_names % append('78195.' // xs)
        call list_names % append('78196.' // xs)
        call list_names % append('78198.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.00014000000000_8)
          call list_density % append(density * 0.00782000000000_8)
          call list_density % append(density * 0.32967000000000_8)
          call list_density % append(density * 0.33832000000000_8)
          call list_density % append(density * 0.25242000000000_8)
          call list_density % append(density * 0.07163000000000_8)
        else
          call list_density % append(density * 0.00013632276599_8)
          call list_density % append(density * 0.00769481514199_8)
          call list_density % append(density * 0.32777509593048_8)
          call list_density % append(density * 0.33811326467036_8)
          call list_density % append(density * 0.25355986681445_8)
          call list_density % append(density * 0.07268889388976_8)
        end if
      end if

    case ('au')
      call list_names % append('79197.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999999847690_8)
      end if

    case ('hg')
      call list_names % append('80196.' // xs)
      call list_names % append('80198.' // xs)
      call list_names % append('80199.' // xs)
      call list_names % append('80200.' // xs)
      call list_names % append('80201.' // xs)
      call list_names % append('80202.' // xs)
      call list_names % append('80204.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00150000000000_8)
        call list_density % append(density * 0.09970000000000_8)
        call list_density % append(density * 0.16870000000000_8)
        call list_density % append(density * 0.23100000000000_8)
        call list_density % append(density * 0.13180000000000_8)
        call list_density % append(density * 0.29860000000000_8)
        call list_density % append(density * 0.06870000000000_8)
      else
        call list_density % append(density * 0.00146542075627_8)
        call list_density % append(density * 0.09839616565781_8)
        call list_density % append(density * 0.16733610259300_8)
        call list_density % append(density * 0.23028407849843_8)
        call list_density % append(density * 0.13204988206361_8)
        call list_density % append(density * 0.30065523704970_8)
        call list_density % append(density * 0.06985881166025_8)
      end if

    case ('tl')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call list_names % append('81000.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 1.00000000000000_8)
        else
          call list_density % append(density * 1.00000000000000_8)
        end if
      else
        call list_names % append('81203.' // xs)
        call list_names % append('81205.' // xs)
        if (density > 0.0) then
          call list_density % append(density * 0.29520000000000_8)
          call list_density % append(density * 0.70480000000000_8)
        else
          call list_density % append(density * 0.29316209302737_8)
          call list_density % append(density * 0.70683845745714_8)
        end if
      end if

    case ('pb')
      call list_names % append('82204.' // xs)
      call list_names % append('82206.' // xs)
      call list_names % append('82207.' // xs)
      call list_names % append('82208.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.01400000000000_8)
        call list_density % append(density * 0.24100000000000_8)
        call list_density % append(density * 0.22100000000000_8)
        call list_density % append(density * 0.52400000000000_8)
      else
        call list_density % append(density * 0.01378196240541_8)
        call list_density % append(density * 0.23957454699469_8)
        call list_density % append(density * 0.22076097111438_8)
        call list_density % append(density * 0.52596412017568_8)
      end if

    case ('bi')
      call list_names % append('83209.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999999377932_8)
      end if

    case ('th')
      call list_names % append('90232.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 0.99999997974470_8)
      end if

    case ('pa')
      call list_names % append('91231.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 1.00000000000000_8)
      else
        call list_density % append(density * 1.00000001731333_8)
      end if

    case ('u')
      call list_names % append('92234.' // xs)
      call list_names % append('92235.' // xs)
      call list_names % append('92238.' // xs)
      if (density > 0.0) then
        call list_density % append(density * 0.00005400000000_8)
        call list_density % append(density * 0.00720400000000_8)
        call list_density % append(density * 0.99274200000000_8)
      else
        call list_density % append(density * 0.00005309527911_8)
        call list_density % append(density * 0.00711365888706_8)
        call list_density % append(density * 0.99283324693309_8)
      end if

    ! End of the natural element cases.
    case default
      message = "Cannot expand element: " // name
      call fatal_error()

    end select

  end subroutine expand_natural_element

end module input_xml
