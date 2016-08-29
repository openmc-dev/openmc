module input_xml

  use hdf5

  use cmfd_input,       only: configure_cmfd
  use constants
  use dict_header,      only: DictIntInt, ElemKeyValueCI
  use distribution_multivariate
  use distribution_univariate
  use endf,             only: reaction_name
  use energy_grid,      only: grid_method, n_log_bins
  use error,            only: fatal_error, warning
  use geometry_header,  only: Cell, Lattice, RectLattice, HexLattice
  use global
  use hdf5_interface
  use list_header,      only: ListChar, ListInt, ListReal
  use mesh_header,      only: RegularMesh
  use multipole,        only: multipole_read
  use output,           only: write_message
  use plot_header
  use random_lcg,       only: prn, seed
  use surface_header
  use set_header,       only: SetChar
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with, tokenize, split_string
  use tally_header,     only: TallyObject
  use tally_filter
  use tally_initialize, only: add_tallies
  use xml_interface

  implicit none
  save

  type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                         ! universe contains

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    call read_settings_xml()
    call read_geometry_xml()
    call read_materials()
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
    integer :: temp_int
    integer :: temp_int_array3(3)
    integer, allocatable :: temp_int_array(:)
    integer(8) :: temp_long
    real(8), allocatable :: temp_real(:)
    integer :: n_tracks
    logical :: file_exists
    character(MAX_FILE_LEN) :: env_variable
    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: filename
    type(Node), pointer :: doc            => null()
    type(Node), pointer :: node_mode      => null()
    type(Node), pointer :: node_source    => null()
    type(Node), pointer :: node_space     => null()
    type(Node), pointer :: node_angle     => null()
    type(Node), pointer :: node_dist      => null()
    type(Node), pointer :: node_cutoff    => null()
    type(Node), pointer :: node_entropy   => null()
    type(Node), pointer :: node_ufs       => null()
    type(Node), pointer :: node_sp        => null()
    type(Node), pointer :: node_output    => null()
    type(Node), pointer :: node_verb      => null()
    type(Node), pointer :: node_res_scat  => null()
    type(Node), pointer :: node_scatterer => null()
    type(Node), pointer :: node_trigger   => null()
    type(Node), pointer :: node_keff_trigger => null()
    type(Node), pointer :: node_vol => null()
    type(NodeList), pointer :: node_scat_list => null()
    type(NodeList), pointer :: node_source_list => null()
    type(NodeList), pointer :: node_vol_list => null()

    ! Check if settings.xml exists
    filename = trim(path_input) // "settings.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      if (run_mode /= MODE_PLOTTING) then
        call fatal_error("Settings XML file '" // trim(filename) // "' does &
             &not exist! In order to run OpenMC, you first need a set of input &
             &files; at a minimum, this includes settings.xml, geometry.xml, &
             &and materials.xml. Please consult the user's guide at &
             &http://mit-crpg.github.io/openmc for further information.")
      else
        ! The settings.xml file is optional if we just want to make a plot.
        return
      end if
    else
      call write_message("Reading settings XML file...", 5)
    end if

    ! Parse settings.xml file
    call open_xmldoc(doc, filename)

    ! Find if a multi-group or continuous-energy simulation is desired
    if (check_for_node(doc, "energy_mode")) then
      call get_node_value(doc, "energy_mode", temp_str)
      temp_str = trim(to_lower(temp_str))
      if (temp_str == "mg" .or. temp_str == "multi-group") then
        run_CE = .false.
      else if (temp_str == "ce" .or. temp_str == "continuous-energy") then
        run_CE = .true.
      end if
    end if

    ! Find cross_sections.xml file -- the first place to look is the
    ! settings.xml file. If no file is found there, then we check the
    ! CROSS_SECTIONS environment variable
    if (.not. check_for_node(doc, "cross_sections")) then
      ! No cross_sections.xml file specified in settings.xml, check
      ! environment variable
      if (run_CE) then
        call get_environment_variable("OPENMC_CROSS_SECTIONS", env_variable)
        if (len_trim(env_variable) == 0) then
          call get_environment_variable("CROSS_SECTIONS", env_variable)
          if (len_trim(env_variable) == 0) then
            call fatal_error("No cross_sections.xml file was specified in &
                 &settings.xml or in the OPENMC_CROSS_SECTIONS environment &
                 &variable. OpenMC needs such a file to identify where to &
                 &find ACE cross section libraries. Please consult the &
                 &user's guide at http://mit-crpg.github.io/openmc for &
                 &information on how to set up ACE cross section libraries.")
          else
            call warning("The CROSS_SECTIONS environment variable is &
                 &deprecated. Please update your environment to use &
                 &OPENMC_CROSS_SECTIONS instead.")
          end if
        end if
        path_cross_sections = trim(env_variable)
      else
        call get_environment_variable("OPENMC_MG_CROSS_SECTIONS", env_variable)
        if (len_trim(env_variable) == 0) then
          call fatal_error("No mgxs.xml file was specified in &
               &settings.xml or in the OPENMC_MG_CROSS_SECTIONS environment &
               &variable. OpenMC needs such a file to identify where to &
               &find ACE cross section libraries. Please consult the user's &
               &guide at http://mit-crpg.github.io/openmc for information on &
               &how to set up ACE cross section libraries.")
        else
          path_cross_sections = trim(env_variable)
        end if
      end if
    else
      call get_node_value(doc, "cross_sections", path_cross_sections)
    end if

    ! Find the windowed multipole library
    if (run_mode /= MODE_PLOTTING) then
      if (.not. check_for_node(doc, "multipole_library")) then
        ! No library location specified in settings.xml, check
        ! environment variable
        call get_environment_variable("OPENMC_MULTIPOLE_LIBRARY", env_variable)
        path_multipole = trim(env_variable)
      else
        call get_node_value(doc, "multipole_library", path_multipole)
      end if
      if (.not. ends_with(path_multipole, "/")) &
           path_multipole = trim(path_multipole) // "/"
    end if

    if (.not. run_CE) then
      ! Scattering Treatments
      if (check_for_node(doc, "max_order")) then
        call get_node_value(doc, "max_order", max_order)
      else
        ! Set to default of largest int - 1, which means to use whatever is
        ! contained in library.
        ! This is largest int - 1 because for legendre scattering, a value of
        ! 1 is added to the order; adding 1 to huge(0) gets you the largest
        ! negative integer, which is not what we want.
        max_order = huge(0) - 1
      end if
    else
      max_order = 0
    end if

    ! Set output directory if a path has been specified on the <output_path>
    ! element
    if (check_for_node(doc, "output_path")) then
      call get_node_value(doc, "output_path", path_output)
      if (.not. ends_with(path_output, "/")) &
           path_output = trim(path_output) // "/"
    end if

    ! Check for a trigger node and get trigger information
    if (check_for_node(doc, "trigger")) then
      call get_node_ptr(doc, "trigger", node_trigger)

      ! Check if trigger(s) are to be turned on
      call get_node_value(node_trigger, "active", temp_str)
      temp_str = trim(to_lower(temp_str))

      if (temp_str == 'true' .or. temp_str == '1') then
        trigger_on = .true.
      elseif (temp_str == 'false' .or. temp_str == '0') then
        trigger_on = .false.
      else
        call fatal_error("Unrecognized trigger active: " // temp_str)
      end if

      if (trigger_on) then

        if (check_for_node(node_trigger, "max_batches") )then
          call get_node_value(node_trigger, "max_batches", n_max_batches)
        else
          call fatal_error("The max_batches must be specified with triggers")
        end if

        ! Get the batch interval to check triggers
        if (.not. check_for_node(node_trigger, "batch_interval"))then
          pred_batches = .true.
        else
          call get_node_value(node_trigger, "batch_interval", temp_int)
          n_batch_interval = temp_int
          if (n_batch_interval <= 0) then
            call fatal_error("The batch interval must be greater than zero")
          end if
        end if
      end if
    end if

    ! Make sure that either eigenvalue or fixed source was specified
    if (.not. check_for_node(doc, "eigenvalue") .and. &
         .not. check_for_node(doc, "fixed_source")) then
      call fatal_error("<eigenvalue> or <fixed_source> not specified.")
    end if

    ! Eigenvalue information
    if (check_for_node(doc, "eigenvalue")) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_EIGENVALUE

      ! Get pointer to eigenvalue XML block
      call get_node_ptr(doc, "eigenvalue", node_mode)

      ! Check number of particles
      if (.not. check_for_node(node_mode, "particles")) then
        call fatal_error("Need to specify number of particles per generation.")
      end if

      ! Get number of particles
      call get_node_value(node_mode, "particles", temp_long)

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = temp_long

      ! Get number of basic batches
      call get_node_value(node_mode, "batches", n_batches)
      if (.not. trigger_on) then
        n_max_batches = n_batches
      end if

      ! Get number of inactive batches
      call get_node_value(node_mode, "inactive", n_inactive)
      n_active = n_batches - n_inactive
      if (check_for_node(node_mode, "generations_per_batch")) then
        call get_node_value(node_mode, "generations_per_batch", gen_per_batch)
      end if

      ! Allocate array for batch keff and entropy
      allocate(k_generation(n_max_batches*gen_per_batch))
      allocate(entropy(n_max_batches*gen_per_batch))
      entropy = ZERO

      ! Get the trigger information for keff
      if (check_for_node(node_mode, "keff_trigger")) then
        call get_node_ptr(node_mode, "keff_trigger", node_keff_trigger)

        if (check_for_node(node_keff_trigger, "type")) then
          call get_node_value(node_keff_trigger, "type", temp_str)
          temp_str = trim(to_lower(temp_str))

          select case (temp_str)
          case ('std_dev')
            keff_trigger % trigger_type = STANDARD_DEVIATION
          case ('variance')
            keff_trigger % trigger_type = VARIANCE
          case ('rel_err')
            keff_trigger % trigger_type = RELATIVE_ERROR
          case default
            call fatal_error("Unrecognized keff trigger type " // temp_str)
          end select

        else
          call fatal_error("Specify keff trigger type in settings XML")
        end if

        if (check_for_node(node_keff_trigger, "threshold")) then
          call get_node_value(node_keff_trigger, "threshold", &
               keff_trigger % threshold)
        else
          call fatal_error("Specify keff trigger threshold in settings XML")
        end if
      end if
    end if

    ! Fixed source calculation information
    if (check_for_node(doc, "fixed_source")) then
      ! Set run mode
      if (run_mode == NONE) run_mode = MODE_FIXEDSOURCE

      ! Get pointer to fixed_source XML block
      call get_node_ptr(doc, "fixed_source", node_mode)

      ! Check number of particles
      if (.not. check_for_node(node_mode, "particles")) then
        call fatal_error("Need to specify number of particles per batch.")
      end if

      ! Get number of particles
      call get_node_value(node_mode, "particles", temp_long)

      ! If the number of particles was specified as a command-line argument, we
      ! don't set it here
      if (n_particles == 0) n_particles = temp_long

      ! Copy batch information
      call get_node_value(node_mode, "batches", n_batches)
      if (.not. trigger_on) then
        n_max_batches = n_batches
      end if
      n_active = n_batches
      n_inactive    = 0
      gen_per_batch = 1
    end if

    ! Check number of active batches, inactive batches, and particles
    if (n_active <= 0) then
      call fatal_error("Number of active batches must be greater than zero.")
    elseif (n_inactive < 0) then
      call fatal_error("Number of inactive batches must be non-negative.")
    elseif (n_particles <= 0) then
      call fatal_error("Number of particles must be greater than zero.")
    end if

    ! Copy random number seed if specified
    if (check_for_node(doc, "seed")) call get_node_value(doc, "seed", seed)

    ! Energy grid methods
    if (check_for_node(doc, "energy_grid")) then
      call get_node_value(doc, "energy_grid", temp_str)
    else
      temp_str = 'logarithm'
    end if
    select case (trim(temp_str))
    case ('nuclide')
      grid_method = GRID_NUCLIDE
    case ('material-union', 'union')
      grid_method = GRID_MAT_UNION
      if (trim(temp_str) == 'union') &
           call warning('Energy grids will be unionized by material.  Global&
           & energy grid unionization is no longer an allowed option.')
    case ('logarithm', 'logarithmic', 'log')
      grid_method = GRID_LOGARITHM
    case default
      call fatal_error("Unknown energy grid method: " // trim(temp_str))
    end select

    ! Number of bins for logarithmic grid
    if (check_for_node(doc, "log_grid_bins")) then
      call get_node_value(doc, "log_grid_bins", n_log_bins)
      if (n_log_bins < 1) then
        call fatal_error("Number of bins for logarithmic grid must be &
             &greater than zero.")
      end if
    else
      n_log_bins = 8000
    end if

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
          call fatal_error("Invalid number of threads: " // to_str(n_threads))
        end if
        call omp_set_num_threads(n_threads)
      end if
#else
      if (master) call warning("Ignoring number of threads.")
#endif
    end if

    ! ==========================================================================
    ! EXTERNAL SOURCE

    ! Get point to list of <source> elements and make sure there is at least one
    call get_node_list(doc, "source", node_source_list)
    n = get_list_size(node_source_list)
    if (n == 0) call fatal_error("No source specified in settings XML file.")

    ! Allocate array for sources
    allocate(external_source(n))

    ! Read each source
    do i = 1, n
      ! Get pointer to source
      call get_list_item(node_source_list, i, node_source)

      ! Check if we want to write out source
      if (check_for_node(node_source, "write_initial")) then
        call get_node_value(node_source, "write_initial", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
             write_initial_source = .true.
      end if

      ! Check for source strength
      if (check_for_node(node_source, "strength")) then
        call get_node_value(node_source, "strength", external_source(i)%strength)
      else
        external_source(i)%strength = ONE
      end if

      ! Check for external source file
      if (check_for_node(node_source, "file")) then
        ! Copy path of source file
        call get_node_value(node_source, "file", path_source)

        ! Check if source file exists
        inquire(FILE=path_source, EXIST=file_exists)
        if (.not. file_exists) then
          call fatal_error("Binary source file '" // trim(path_source) &
               &// "' does not exist!")
        end if

      else

        ! Spatial distribution for external source
        if (check_for_node(node_source, "space")) then

          ! Get pointer to spatial distribution
          call get_node_ptr(node_source, "space", node_space)

          ! Check for type of spatial distribution
          type = ''
          if (check_for_node(node_space, "type")) &
               call get_node_value(node_space, "type", type)
          select case (to_lower(type))
          case ('cartesian')
            allocate(CartesianIndependent :: external_source(i)%space)

          case ('box')
            allocate(SpatialBox :: external_source(i)%space)

          case ('fission')
            allocate(SpatialBox :: external_source(i)%space)
            select type(space => external_source(i)%space)
            type is (SpatialBox)
              space%only_fissionable = .true.
            end select

          case ('point')
            allocate(SpatialPoint :: external_source(i)%space)

          case default
            call fatal_error("Invalid spatial distribution for external source: "&
                 // trim(type))
          end select

          select type (space => external_source(i)%space)
          type is (CartesianIndependent)
            ! Read distribution for x coordinate
            if (check_for_node(node_space, "x")) then
              call get_node_ptr(node_space, "x", node_dist)
              call distribution_from_xml(space%x, node_dist)
            else
              allocate(Discrete :: space%x)
              select type (dist => space%x)
              type is (Discrete)
                allocate(dist%x(1), dist%p(1))
                dist%x(1) = ZERO
                dist%p(1) = ONE
              end select
            end if

            ! Read distribution for y coordinate
            if (check_for_node(node_space, "y")) then
              call get_node_ptr(node_space, "y", node_dist)
              call distribution_from_xml(space%y, node_dist)
            else
              allocate(Discrete :: space%y)
              select type (dist => space%y)
              type is (Discrete)
                allocate(dist%x(1), dist%p(1))
                dist%x(1) = ZERO
                dist%p(1) = ONE
              end select
            end if

            if (check_for_node(node_space, "z")) then
              call get_node_ptr(node_space, "z", node_dist)
              call distribution_from_xml(space%z, node_dist)
            else
              allocate(Discrete :: space%z)
              select type (dist => space%z)
              type is (Discrete)
                allocate(dist%x(1), dist%p(1))
                dist%x(1) = ZERO
                dist%p(1) = ONE
              end select
            end if

          type is (SpatialBox)
            ! Make sure correct number of parameters are given
            if (get_arraysize_double(node_space, "parameters") /= 6) then
              call fatal_error('Box/fission spatial source must have &
                   &six parameters specified.')
            end if

            ! Read lower-right/upper-left coordinates
            allocate(temp_real(6))
            call get_node_array(node_space, "parameters", temp_real)
            space%lower_left(:) = temp_real(1:3)
            space%upper_right(:) = temp_real(4:6)
            deallocate(temp_real)

          type is (SpatialPoint)
            ! Make sure correct number of parameters are given
            if (get_arraysize_double(node_space, "parameters") /= 3) then
              call fatal_error('Point spatial source must have &
                   &three parameters specified.')
            end if

            ! Read location of point source
            allocate(temp_real(3))
            call get_node_array(node_space, "parameters", temp_real)
            space%xyz(:) = temp_real
            deallocate(temp_real)

          end select

        else
          call fatal_error("No spatial distribution specified for external &
               &source.")
        end if

        ! Determine external source angular distribution
        if (check_for_node(node_source, "angle")) then

          ! Get pointer to angular distribution
          call get_node_ptr(node_source, "angle", node_angle)

          ! Check for type of angular distribution
          type = ''
          if (check_for_node(node_angle, "type")) &
               call get_node_value(node_angle, "type", type)
          select case (to_lower(type))
          case ('isotropic')
            allocate(Isotropic :: external_source(i)%angle)

          case ('monodirectional')
            allocate(Monodirectional :: external_source(i)%angle)

          case ('mu-phi')
            allocate(PolarAzimuthal :: external_source(i)%angle)

          case default
            call fatal_error("Invalid angular distribution for external source: "&
                 // trim(type))
          end select

          ! Read reference directional unit vector
          if (check_for_node(node_angle, "reference_uvw")) then
            n = get_arraysize_double(node_angle, "reference_uvw")
            if (n /= 3) then
              call fatal_error('Angular distribution reference direction must have &
                   &three parameters specified.')
            end if
            call get_node_array(node_angle, "reference_uvw", &
                 external_source(i)%angle%reference_uvw)
          else
            ! By default, set reference unit vector to be positive z-direction
            external_source(i)%angle%reference_uvw(:) = [ZERO, ZERO, ONE]
          end if

          ! Read parameters for angle distribution
          select type (angle => external_source(i)%angle)
          type is (Monodirectional)
            call get_node_array(node_angle, "reference_uvw", &
                 external_source(i)%angle%reference_uvw)

          type is (PolarAzimuthal)
            if (check_for_node(node_angle, "mu")) then
              call get_node_ptr(node_angle, "mu", node_dist)
              call distribution_from_xml(angle%mu, node_dist)
            else
              allocate(Uniform :: angle%mu)
              select type (mu => angle%mu)
              type is (Uniform)
                mu%a = -ONE
                mu%b = ONE
              end select
            end if

            if (check_for_node(node_angle, "phi")) then
              call get_node_ptr(node_angle, "phi", node_dist)
              call distribution_from_xml(angle%phi, node_dist)
            else
              allocate(Uniform :: angle%phi)
              select type (phi => angle%phi)
              type is (Uniform)
                phi%a = ZERO
                phi%b = TWO*PI
              end select
            end if
          end select

        else
          ! Set default angular distribution isotropic
          allocate(Isotropic :: external_source(i)%angle)
          external_source(i)%angle%reference_uvw(:) = [ZERO, ZERO, ONE]
        end if

        ! Determine external source energy distribution
        if (check_for_node(node_source, "energy")) then
          call get_node_ptr(node_source, "energy", node_dist)
          call distribution_from_xml(external_source(i)%energy, node_dist)
        else
          ! Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
          allocate(Watt :: external_source(i)%energy)
          select type(energy => external_source(i)%energy)
          type is (Watt)
            energy%a = 0.988_8
            energy%b = 2.249_8
          end select
        end if
      end if
    end do

    ! Survival biasing
    if (check_for_node(doc, "survival_biasing")) then
      call get_node_value(doc, "survival_biasing", temp_str)
      temp_str = to_lower(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
           survival_biasing = .true.
    end if

    ! Probability tables
    if (check_for_node(doc, "ptables")) then
      call get_node_value(doc, "ptables", temp_str)
      temp_str = to_lower(temp_str)
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
        call fatal_error("Number of integers specified in 'track' is not &
             &divisible by 3.  Please provide 3 integers per particle to be &
             &tracked.")
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
        call fatal_error("Need to specify (x,y,z) coordinates of lower-left &
             &corner of Shannon entropy mesh.")
      elseif (get_arraysize_double(node_entropy, "upper_right") /= 3) then
        call fatal_error("Need to specify (x,y,z) coordinates of upper-right &
             &corner of Shannon entropy mesh.")
      end if

      ! Allocate mesh object and coordinates on mesh
      allocate(entropy_mesh)
      allocate(entropy_mesh % lower_left(3))
      allocate(entropy_mesh % upper_right(3))
      allocate(entropy_mesh % width(3))

      ! Copy values
      call get_node_array(node_entropy, "lower_left", &
           entropy_mesh % lower_left)
      call get_node_array(node_entropy, "upper_right", &
           entropy_mesh % upper_right)

      ! Check on values provided
      if (.not. all(entropy_mesh % upper_right > entropy_mesh % lower_left)) &
           &then
        call fatal_error("Upper-right coordinate must be greater than &
             &lower-left coordinate for Shannon entropy mesh.")
      end if

      ! Check if dimensions were specified -- if not, they will be calculated
      ! automatically upon first entry into shannon_entropy
      if (check_for_node(node_entropy, "dimension")) then

        ! If so, make sure proper number of values were given
        if (get_arraysize_integer(node_entropy, "dimension") /= 3) then
          call fatal_error("Dimension of entropy mesh must be given as three &
               &integers.")
        end if

        ! Allocate dimensions
        entropy_mesh % n_dimension = 3
        allocate(entropy_mesh % dimension(3))

        ! Copy dimensions
        call get_node_array(node_entropy, "dimension", entropy_mesh % dimension)

        ! Calculate width
        entropy_mesh % width = (entropy_mesh % upper_right - &
             entropy_mesh % lower_left) / entropy_mesh % dimension

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
        call fatal_error("Need to specify (x,y,z) coordinates of lower-left &
             &corner of UFS mesh.")
      elseif (get_arraysize_double(node_ufs, "upper_right") /= 3) then
        call fatal_error("Need to specify (x,y,z) coordinates of upper-right &
             &corner of UFS mesh.")
      elseif (get_arraysize_integer(node_ufs, "dimension") /= 3) then
        call fatal_error("Dimension of UFS mesh must be given as three &
             &integers.")
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
        call fatal_error("Upper-right coordinate must be greater than &
             &lower-left coordinate for UFS mesh.")
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
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') source_separate = .true.
      end if
      if (check_for_node(node_sp, "write")) then
        call get_node_value(node_sp, "write", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'false' .or. &
             trim(temp_str) == '0') source_write = .false.
      end if
      if (check_for_node(node_sp, "overwrite_latest")) then
        call get_node_value(node_sp, "overwrite_latest", temp_str)
        temp_str = to_lower(temp_str)
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
          call fatal_error('Sourcepoint batches are not a subset&
               & of statepoint batches.')
        end if
      end do
    end if

    ! Check if the user has specified to not reduce tallies at the end of every
    ! batch
    if (check_for_node(doc, "no_reduce")) then
      call get_node_value(doc, "no_reduce", temp_str)
      temp_str = to_lower(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
           reduce_tallies = .false.
    end if

    ! Check if the user has specified to use confidence intervals for
    ! uncertainties rather than standard deviations
    if (check_for_node(doc, "confidence_intervals")) then
      call get_node_value(doc, "confidence_intervals", temp_str)
      temp_str = to_lower(temp_str)
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
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'false' .or. &
             trim(temp_str) == '0') output_summary = .false.
      end if

      ! Check for cross sections option
      if (check_for_node(node_output, "cross_sections")) then
        call get_node_value(node_output, "cross_sections", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'true' .or. &
             trim(temp_str) == '1') output_xs = .true.
      end if

      ! Check for ASCII tallies output option
      if (check_for_node(node_output, "tallies")) then
        call get_node_value(node_output, "tallies", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'false' .or. &
             trim(temp_str) == '0') output_tallies = .false.
      end if
    end if

    ! Check for cmfd run
    if (check_for_node(doc, "run_cmfd")) then
      call get_node_value(doc, "run_cmfd", temp_str)
      temp_str = to_lower(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') then
        cmfd_run = .true.
      end if
    end if

    ! Resonance scattering parameters
    if (check_for_node(doc, "resonance_scattering")) then
      call get_node_ptr(doc, "resonance_scattering", node_res_scat)
      call get_node_list(node_res_scat, "scatterer", node_scat_list)

      ! check that a nuclide is specified
      if (get_list_size(node_scat_list) >= 1) then
        treat_res_scat = .true.
        n_res_scatterers_total = get_list_size(node_scat_list)

        ! store 0K info for resonant scatterers
        allocate(nuclides_0K(n_res_scatterers_total))
        do i = 1, n_res_scatterers_total
          call get_list_item(node_scat_list, i, node_scatterer)

          ! check to make sure a nuclide is specified
          if (.not. check_for_node(node_scatterer, "nuclide")) then
            call fatal_error("No nuclide specified for scatterer " &
                 // trim(to_str(i)) // " in settings.xml file!")
          end if
          call get_node_value(node_scatterer, "nuclide", &
               nuclides_0K(i) % nuclide)

          if (check_for_node(node_scatterer, "method")) then
            call get_node_value(node_scatterer, "method", &
                 nuclides_0K(i) % scheme)
          end if

          ! check to make sure xs name for which method is applied is given
          if (.not. check_for_node(node_scatterer, "xs_label")) then
            call fatal_error("Must specify the temperature dependent name of &
                 &scatterer " // trim(to_str(i)) &
                 // " given in cross_sections.xml")
          end if
          call get_node_value(node_scatterer, "xs_label", &
               nuclides_0K(i) % name)

          ! check to make sure 0K xs name for which method is applied is given
          if (.not. check_for_node(node_scatterer, "xs_label_0K")) then
            call fatal_error("Must specify the 0K name of scatterer " &
                 // trim(to_str(i)) // " given in cross_sections.xml")
          end if
          call get_node_value(node_scatterer, "xs_label_0K", &
               nuclides_0K(i) % name_0K)

          if (check_for_node(node_scatterer, "E_min")) then
            call get_node_value(node_scatterer, "E_min", &
                 nuclides_0K(i) % E_min)
          end if

          ! check that E_min is non-negative
          if (nuclides_0K(i) % E_min < ZERO) then
            call fatal_error("Lower resonance scattering energy bound is &
                 &negative")
          end if

          if (check_for_node(node_scatterer, "E_max")) then
            call get_node_value(node_scatterer, "E_max", &
                 nuclides_0K(i) % E_max)
          end if

          ! check that E_max is not less than E_min
          if (nuclides_0K(i) % E_max < nuclides_0K(i) % E_min) then
            call fatal_error("Lower resonance scattering energy bound exceeds &
                 &upper")
          end if

          nuclides_0K(i) % nuclide = trim(nuclides_0K(i) % nuclide)
          nuclides_0K(i) % scheme  = to_lower(trim(nuclides_0K(i) % scheme))
          nuclides_0K(i) % name    = trim(nuclides_0K(i) % name)
          nuclides_0K(i) % name_0K = trim(nuclides_0K(i) % name_0K)
        end do
      else
        call fatal_error("No resonant scatterers are specified within the &
             &resonance_scattering element in settings.xml")
      end if
    end if

    ! Natural element expansion option
    if (check_for_node(doc, "natural_elements")) then
      call get_node_value(doc, "natural_elements", temp_str)
      select case (to_lower(temp_str))
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
        call fatal_error("Unknown natural element expansion option: " &
             // trim(temp_str))
      end select
    end if

    ! Check to see if windowed multipole functionality is requested
    if (check_for_node(doc, "use_windowed_multipole")) then
      call get_node_value(doc, "use_windowed_multipole", temp_str)
      select case (to_lower(temp_str))
      case ('true', '1')
        multipole_active = .true.
      case ('false', '0')
        multipole_active = .false.
      case default
        call fatal_error("Unrecognized value for <use_windowed_multipole> in &
             &settings.xml")
      end select
    end if

    call get_node_list(doc, "volume_calc", node_vol_list)
    n = get_list_size(node_vol_list)
    allocate(volume_calcs(n))
    do i = 1, n
      call get_list_item(node_vol_list, i, node_vol)
      call volume_calcs(i) % from_xml(node_vol)
    end do

    ! Close settings XML file
    call close_xmldoc(doc)

  end subroutine read_settings_xml

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i, j, k, m, i_x, i_a, input_index
    integer :: n, n_mats, n_x, n_y, n_z, n_rings, n_rlats, n_hlats
    integer :: universe_num
    integer :: n_cells_in_univ
    integer :: coeffs_reqd
    integer :: i_xmin, i_xmax, i_ymin, i_ymax, i_zmin, i_zmax
    real(8) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer, allocatable :: temp_int_array(:)
    real(8) :: phi, theta, psi
    real(8), allocatable :: coeffs(:)
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    character(REGION_SPEC_LEN) :: region_spec
    type(Cell),     pointer :: c
    class(Surface), pointer :: s
    class(Lattice), pointer :: lat
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_cell => null()
    type(Node), pointer :: node_surf => null()
    type(Node), pointer :: node_lat => null()
    type(NodeList), pointer :: node_cell_list => null()
    type(NodeList), pointer :: node_surf_list => null()
    type(NodeList), pointer :: node_rlat_list => null()
    type(NodeList), pointer :: node_hlat_list => null()
    type(VectorInt) :: tokens
    type(VectorInt) :: rpn

    ! Display output message
    call write_message("Reading geometry XML file...", 5)

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Geometry XML file '" // trim(filename) // "' does not &
           &exist!")
    end if

    ! Parse geometry.xml file
    call open_xmldoc(doc, filename)

    ! Get pointer to list of XML <cell>
    call get_node_list(doc, "cell", node_cell_list)

    ! Get number of <cell> tags
    n_cells = get_list_size(node_cell_list)

    ! Check for no cells
    if (n_cells == 0) then
      call fatal_error("No cells found in geometry.xml!")
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

      ! Initialize distribcell instances and distribcell index
      c % instances = 0
      c % distribcell_index = NONE

      ! Get pointer to i-th cell node
      call get_list_item(node_cell_list, i, node_cell)

      ! Copy data into cells
      if (check_for_node(node_cell, "id")) then
        call get_node_value(node_cell, "id", c % id)
      else
        call fatal_error("Must specify id of cell in geometry XML file.")
      end if

      ! Copy cell name
      if (check_for_node(node_cell, "name")) then
        call get_node_value(node_cell, "name", c % name)
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
        call fatal_error("Two or more cells use the same unique ID: " &
             // to_str(c % id))
      end if

      ! Read material
      if (check_for_node(node_cell, "material")) then
        n_mats = get_arraysize_string(node_cell, "material")

        if (n_mats > 0) then
          allocate(sarray(n_mats))
          call get_node_array(node_cell, "material", sarray)

          allocate(c % material(n_mats))
          do j = 1, n_mats
            select case(trim(to_lower(sarray(j))))
            case ('void')
              c % material(j) = MATERIAL_VOID
            case default
              c % material(j) = int(str_to_int(sarray(j)), 4)

              ! Check for error
              if (c % material(j) == ERROR_INT) then
                call fatal_error("Invalid material specified on cell " &
                     // to_str(c % id))
              end if
            end select
          end do

          deallocate(sarray)

        else
          allocate(c % material(1))
          c % material(1) = NONE
        end if

      else
        allocate(c % material(1))
        c % material(1) = NONE
      end if

      ! Check to make sure that either material or fill was specified
      if (c % material(1) == NONE .and. c % fill == NONE) then
        call fatal_error("Neither material nor fill was specified for cell " &
             // trim(to_str(c % id)))
      end if

      ! Check to make sure that both material and fill haven't been
      ! specified simultaneously
      if (c % material(1) /= NONE .and. c % fill /= NONE) then
        call fatal_error("Cannot specify material and fill simultaneously")
      end if

      ! Check for region specification (also under deprecated name surfaces)
      region_spec = ''
      if (check_for_node(node_cell, "surfaces")) then
        call warning("The use of 'surfaces' is deprecated and will be &
             &disallowed in a future release.  Use 'region' instead. The &
             &openmc-update-inputs utility can be used to automatically &
             &update geometry.xml files.")
        call get_node_value(node_cell, "surfaces", region_spec)
      elseif (check_for_node(node_cell, "region")) then
        call get_node_value(node_cell, "region", region_spec)
      end if

      if (len_trim(region_spec) > 0) then
        ! Create surfaces array from string
        call tokenize(region_spec, tokens)

        ! Use shunting-yard algorithm to determine RPN for surface algorithm
        call generate_rpn(c%id, tokens, rpn)

        ! Copy region spec and RPN form to cell arrays
        allocate(c % region(tokens%size()))
        allocate(c % rpn(rpn%size()))
        c % region(:) = tokens%data(1:tokens%size())
        c % rpn(:) = rpn%data(1:rpn%size())

        call tokens%clear()
        call rpn%clear()
      end if
      if (.not. allocated(c%region)) allocate(c%region(0))
      if (.not. allocated(c%rpn)) allocate(c%rpn(0))

      ! Check if this is a simple cell
      if (any(c%rpn == OP_COMPLEMENT) .or. any(c%rpn == OP_UNION)) then
        c%simple = .false.
      else
        c%simple = .true.
      end if

      ! Rotation matrix
      if (check_for_node(node_cell, "rotation")) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill == NONE) then
          call fatal_error("Cannot apply a rotation to cell " // trim(to_str(&
               &c % id)) // " because it is not filled with another universe")
        end if

        ! Read number of rotation parameters
        n = get_arraysize_double(node_cell, "rotation")
        if (n /= 3) then
          call fatal_error("Incorrect number of rotation parameters on cell " &
               // to_str(c % id))
        end if

        ! Copy rotation angles in x,y,z directions
        allocate(c % rotation(3))
        call get_node_array(node_cell, "rotation", c % rotation)
        phi   = -c % rotation(1) * PI/180.0_8
        theta = -c % rotation(2) * PI/180.0_8
        psi   = -c % rotation(3) * PI/180.0_8

        ! Calculate rotation matrix based on angles given
        allocate(c % rotation_matrix(3,3))
        c % rotation_matrix = reshape((/ &
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
          call fatal_error("Cannot apply a translation to cell " &
               // trim(to_str(c % id)) // " because it is not filled with &
               &another universe")
        end if

        ! Read number of translation parameters
        n = get_arraysize_double(node_cell, "translation")
        if (n /= 3) then
          call fatal_error("Incorrect number of translation parameters on &
               &cell " // to_str(c % id))
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        call get_node_array(node_cell, "translation", c % translation)
      end if

      ! Read cell temperatures.  If the temperature is not specified, set it to
      ! ERROR_REAL for now.  During initialization we'll replace ERROR_REAL with
      ! the temperature from the material data.
      if (.not. run_CE) then
        ! Cell temperatures are not used for MG mode.
        allocate(c % sqrtkT(1))
        c % sqrtkT(1) = ZERO
      else if (check_for_node(node_cell, "temperature")) then
        n = get_arraysize_double(node_cell, "temperature")
        if (n > 0) then
          ! Make sure this is a "normal" cell.
          if (c % material(1) == NONE) call fatal_error("Cell " &
               // trim(to_str(c % id)) // " was specified with a temperature &
               &but no material. Temperature specification is only valid for &
               &cells filled with a material.")

          ! Copy in temperatures
          allocate(c % sqrtkT(n))
          call get_node_array(node_cell, "temperature", c % sqrtkT)

          ! Make sure all temperatues are positive
          do j = 1, size(c % sqrtkT)
            if (c % sqrtkT(j) < ZERO) call fatal_error("Cell " &
                 // trim(to_str(c % id)) // " was specified with a negative &
                 &temperature. All cell temperatures must be non-negative.")
          end do

          ! Convert to sqrt(kT)
          c % sqrtkT(:) = sqrt(K_BOLTZMANN * c % sqrtkT(:))
        else
          allocate(c % sqrtkT(1))
          c % sqrtkT(1) = ERROR_REAL
        end if
      else
        allocate(c % sqrtkT(1))
        c % sqrtkT = ERROR_REAL
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
      call fatal_error("No surfaces found in geometry.xml!")
    end if

    xmin = INFINITY
    xmax = -INFINITY
    ymin = INFINITY
    ymax = -INFINITY
    zmin = INFINITY
    zmax = -INFINITY

    ! Allocate cells array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      ! Get pointer to i-th surface node
      call get_list_item(node_surf_list, i, node_surf)

      ! Copy and interpret surface type
      word = ''
      if (check_for_node(node_surf, "type")) &
           call get_node_value(node_surf, "type", word)
      select case(to_lower(word))
      case ('x-plane')
        coeffs_reqd  = 1
        allocate(SurfaceXPlane :: surfaces(i)%obj)
      case ('y-plane')
        coeffs_reqd  = 1
        allocate(SurfaceYPlane :: surfaces(i)%obj)
      case ('z-plane')
        coeffs_reqd  = 1
        allocate(SurfaceZPlane :: surfaces(i)%obj)
      case ('plane')
        coeffs_reqd  = 4
        allocate(SurfacePlane :: surfaces(i)%obj)
      case ('x-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceXCylinder :: surfaces(i)%obj)
      case ('y-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceYCylinder :: surfaces(i)%obj)
      case ('z-cylinder')
        coeffs_reqd  = 3
        allocate(SurfaceZCylinder :: surfaces(i)%obj)
      case ('sphere')
        coeffs_reqd  = 4
        allocate(SurfaceSphere :: surfaces(i)%obj)
      case ('x-cone')
        coeffs_reqd  = 4
        allocate(SurfaceXCone :: surfaces(i)%obj)
      case ('y-cone')
        coeffs_reqd  = 4
        allocate(SurfaceYCone :: surfaces(i)%obj)
      case ('z-cone')
        coeffs_reqd  = 4
        allocate(SurfaceZCone :: surfaces(i)%obj)
      case ('quadric')
        coeffs_reqd  = 10
        allocate(SurfaceQuadric :: surfaces(i)%obj)
      case default
        call fatal_error("Invalid surface type: " // trim(word))
      end select

      s => surfaces(i)%obj

      ! Copy data into cells
      if (check_for_node(node_surf, "id")) then
        call get_node_value(node_surf, "id", s%id)
      else
        call fatal_error("Must specify id of surface in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (surface_dict % has_key(s%id)) then
        call fatal_error("Two or more surfaces use the same unique ID: " &
             // to_str(s%id))
      end if

      ! Copy surface name
      if (check_for_node(node_surf, "name")) then
        call get_node_value(node_surf, "name", s%name)
      end if

      ! Check to make sure that the proper number of coefficients
      ! have been specified for the given type of surface. Then copy
      ! surface coordinates.

      n = get_arraysize_double(node_surf, "coeffs")
      if (n < coeffs_reqd) then
        call fatal_error("Not enough coefficients specified for surface: " &
             // trim(to_str(s%id)))
      elseif (n > coeffs_reqd) then
        call fatal_error("Too many coefficients specified for surface: " &
             // trim(to_str(s%id)))
      end if

      allocate(coeffs(n))
      call get_node_array(node_surf, "coeffs", coeffs)

      select type(s)
      type is (SurfaceXPlane)
        s%x0 = coeffs(1)

        ! Determine outer surfaces
        xmin = min(xmin, s % x0)
        xmax = max(xmax, s % x0)
        if (xmin == s % x0) i_xmin = i
        if (xmax == s % x0) i_xmax = i
      type is (SurfaceYPlane)
        s%y0 = coeffs(1)

        ! Determine outer surfaces
        ymin = min(ymin, s % y0)
        ymax = max(ymax, s % y0)
        if (ymin == s % y0) i_ymin = i
        if (ymax == s % y0) i_ymax = i
      type is (SurfaceZPlane)
        s%z0 = coeffs(1)

        ! Determine outer surfaces
        zmin = min(zmin, s % z0)
        zmax = max(zmax, s % z0)
        if (zmin == s % z0) i_zmin = i
        if (zmax == s % z0) i_zmax = i
      type is (SurfacePlane)
        s%A = coeffs(1)
        s%B = coeffs(2)
        s%C = coeffs(3)
        s%D = coeffs(4)
      type is (SurfaceXCylinder)
        s%y0 = coeffs(1)
        s%z0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceYCylinder)
        s%x0 = coeffs(1)
        s%z0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceZCylinder)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%r = coeffs(3)
      type is (SurfaceSphere)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r = coeffs(4)
      type is (SurfaceXCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceYCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceZCone)
        s%x0 = coeffs(1)
        s%y0 = coeffs(2)
        s%z0 = coeffs(3)
        s%r2 = coeffs(4)
      type is (SurfaceQuadric)
        s%A = coeffs(1)
        s%B = coeffs(2)
        s%C = coeffs(3)
        s%D = coeffs(4)
        s%E = coeffs(5)
        s%F = coeffs(6)
        s%G = coeffs(7)
        s%H = coeffs(8)
        s%J = coeffs(9)
        s%K = coeffs(10)
      end select

      ! No longer need coefficients
      deallocate(coeffs)

      ! Boundary conditions
      word = ''
      if (check_for_node(node_surf, "boundary")) &
           call get_node_value(node_surf, "boundary", word)
      select case (to_lower(word))
      case ('transmission', 'transmit', '')
        s%bc = BC_TRANSMIT
      case ('vacuum')
        s%bc = BC_VACUUM
        boundary_exists = .true.
      case ('reflective', 'reflect', 'reflecting')
        s%bc = BC_REFLECT
        boundary_exists = .true.
      case ('periodic')
        s%bc = BC_PERIODIC
        boundary_exists = .true.

        ! Check for specification of periodic surface
        if (check_for_node(node_surf, "periodic_surface_id")) then
          call get_node_value(node_surf, "periodic_surface_id", &
               s % i_periodic)
        end if
      case default
        call fatal_error("Unknown boundary condition '" // trim(word) // &
             &"' specified on surface " // trim(to_str(s%id)))
      end select
      ! Add surface to dictionary
      call surface_dict % add_key(s%id, i)
    end do

    ! Check to make sure a boundary condition was applied to at least one
    ! surface
    if (.not. boundary_exists) then
      call fatal_error("No boundary conditions were applied to any surfaces!")
    end if

    ! Determine opposite side for periodic boundaries
    do i = 1, size(surfaces)
      if (surfaces(i) % obj % bc == BC_PERIODIC) then
        select type (surf => surfaces(i) % obj)
        type is (SurfaceXPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_xmin) then
              surf % i_periodic = i_xmax
            elseif (i == i_xmax) then
              surf % i_periodic = i_xmin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        type is (SurfaceYPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_ymin) then
              surf % i_periodic = i_ymax
            elseif (i == i_ymax) then
              surf % i_periodic = i_ymin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        type is (SurfaceZPlane)
          if (surf % i_periodic == NONE) then
            if (i == i_zmin) then
              surf % i_periodic = i_zmax
            elseif (i == i_zmax) then
              surf % i_periodic = i_zmin
            else
              call fatal_error("Periodic boundary condition applied to &
                   &interior surface.")
            end if
          else
            surf % i_periodic = surface_dict % get_key(surf % i_periodic)
          end if

        class default
          call fatal_error("Periodic boundary condition applied to &
               &non-planar surface.")
        end select

        ! Make sure opposite surface is also periodic
        associate (surf => surfaces(i) % obj)
          if (surfaces(surf % i_periodic) % obj % bc /= BC_PERIODIC) then
            call fatal_error("Could not find matching surface for periodic &
                 &boundary on surface " // trim(to_str(surf % id)) // ".")
          end if
        end associate
      end if
    end do

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    ! Get pointer to list of XML <lattice>
    call get_node_list(doc, "lattice", node_rlat_list)
    call get_node_list(doc, "hex_lattice", node_hlat_list)

    ! Allocate lattices array
    n_rlats = get_list_size(node_rlat_list)
    n_hlats = get_list_size(node_hlat_list)
    n_lattices = n_rlats + n_hlats
    allocate(lattices(n_lattices))

    RECT_LATTICES: do i = 1, n_rlats
      allocate(RectLattice::lattices(i) % obj)
      lat => lattices(i) % obj
      select type(lat)
      type is (RectLattice)

      ! Get pointer to i-th lattice
      call get_list_item(node_rlat_list, i, node_lat)

      ! ID of lattice
      if (check_for_node(node_lat, "id")) then
        call get_node_value(node_lat, "id", lat % id)
      else
        call fatal_error("Must specify id of lattice in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        call fatal_error("Two or more lattices use the same unique ID: " &
             // to_str(lat % id))
      end if

      ! Copy lattice name
      if (check_for_node(node_lat, "name")) then
        call get_node_value(node_lat, "name", lat % name)
      end if

      ! Read number of lattice cells in each dimension
      n = get_arraysize_integer(node_lat, "dimension")
      if (n == 2) then
        call get_node_array(node_lat, "dimension", lat % n_cells(1:2))
        lat % n_cells(3) = 1
        lat % is_3d = .false.
      else if (n == 3) then
        call get_node_array(node_lat, "dimension", lat % n_cells)
        lat % is_3d = .true.
      else
        call fatal_error("Rectangular lattice must be two or three dimensions.")
      end if

      ! Read lattice lower-left location
      if (get_arraysize_double(node_lat, "lower_left") /= n) then
        call fatal_error("Number of entries on <lower_left> must be the same &
             &as the number of entries on <dimension>.")
      end if

      allocate(lat % lower_left(n))
      call get_node_array(node_lat, "lower_left", lat % lower_left)

      ! Read lattice pitches.
      ! TODO: Remove this deprecation warning in a future release.
      if (check_for_node(node_lat, "width")) then
        call warning("The use of 'width' is deprecated and will be disallowed &
             &in a future release.  Use 'pitch' instead.  The utility openmc/&
             &src/utils/update_inputs.py can be used to automatically update &
             &geometry.xml files.")
        if (get_arraysize_double(node_lat, "width") /= n) then
          call fatal_error("Number of entries on <pitch> must be the same as &
               &the number of entries on <dimension>.")
        end if

      else if (get_arraysize_double(node_lat, "pitch") /= n) then
        call fatal_error("Number of entries on <pitch> must be the same as &
             &the number of entries on <dimension>.")
      end if

      allocate(lat % pitch(n))
      ! TODO: Remove the 'width' code in a future release.
      if (check_for_node(node_lat, "width")) then
        call get_node_array(node_lat, "width", lat % pitch)
      else
        call get_node_array(node_lat, "pitch", lat % pitch)
      end if

      ! TODO: Remove deprecation warning in a future release.
      if (check_for_node(node_lat, "type")) then
        call warning("The use of 'type' is no longer needed.  The utility &
             &openmc/src/utils/update_inputs.py can be used to automatically &
             &update geometry.xml files.")
      end if

      ! Copy number of dimensions
      n_x = lat % n_cells(1)
      n_y = lat % n_cells(2)
      n_z = lat % n_cells(3)
      allocate(lat % universes(n_x, n_y, n_z))

      ! Check that number of universes matches size
      n = get_arraysize_integer(node_lat, "universes")
      if (n /= n_x*n_y*n_z) then
        call fatal_error("Number of universes on <universes> does not match &
             &size of lattice " // trim(to_str(lat % id)) // ".")
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = &
                 &temp_int_array(j + n_x*k + n_x*n_y*(m-1))
          end do
        end do
      end do
      deallocate(temp_int_array)

      ! Read outer universe for area outside lattice.
      lat % outer = NO_OUTER_UNIVERSE
      if (check_for_node(node_lat, "outer")) then
        call get_node_value(node_lat, "outer", lat % outer)
      end if

      ! Check for 'outside' nodes which are no longer supported.
      if (check_for_node(node_lat, "outside")) then
        call fatal_error("The use of 'outside' in lattices is no longer &
             &supported.  Instead, use 'outer' which defines a universe rather &
             &than a material.  The utility openmc/src/utils/update_inputs.py &
             &can be used automatically replace 'outside' with 'outer'.")
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, i)

      end select
    end do RECT_LATTICES

    HEX_LATTICES: do i = 1, n_hlats
      allocate(HexLattice::lattices(n_rlats + i) % obj)
      lat => lattices(n_rlats + i) % obj
      select type (lat)
      type is (HexLattice)

      ! Get pointer to i-th lattice
      call get_list_item(node_hlat_list, i, node_lat)

      ! ID of lattice
      if (check_for_node(node_lat, "id")) then
        call get_node_value(node_lat, "id", lat % id)
      else
        call fatal_error("Must specify id of lattice in geometry XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (lattice_dict % has_key(lat % id)) then
        call fatal_error("Two or more lattices use the same unique ID: " &
             // to_str(lat % id))
      end if

      ! Copy lattice name
      if (check_for_node(node_lat, "name")) then
        call get_node_value(node_lat, "name", lat % name)
      end if

      ! Read number of lattice cells in each dimension
      call get_node_value(node_lat, "n_rings", lat % n_rings)
      if (check_for_node(node_lat, "n_axial")) then
        call get_node_value(node_lat, "n_axial", lat % n_axial)
        lat % is_3d = .true.
      else
        lat % n_axial = 1
        lat % is_3d = .false.
      end if

      ! Read lattice lower-left location
      n = get_arraysize_double(node_lat, "center")
      if (lat % is_3d .and. n /= 3) then
        call fatal_error("A hexagonal lattice with <n_axial> must have &
             &<center> specified by 3 numbers.")
      else if ((.not. lat % is_3d) .and. n /= 2) then
        call fatal_error("A hexagonal lattice without <n_axial> must have &
             &<center> specified by 2 numbers.")
      end if

      allocate(lat % center(n))
      call get_node_array(node_lat, "center", lat % center)

      ! Read lattice pitches
      n = get_arraysize_double(node_lat, "pitch")
      if (lat % is_3d .and. n /= 2) then
        call fatal_error("A hexagonal lattice with <n_axial> must have <pitch> &
              &specified by 2 numbers.")
      else if ((.not. lat % is_3d) .and. n /= 1) then
        call fatal_error("A hexagonal lattice without <n_axial> must have &
             &<pitch> specified by 1 number.")
      end if

      allocate(lat % pitch(n))
      call get_node_array(node_lat, "pitch", lat % pitch)

      ! Copy number of dimensions
      n_rings = lat % n_rings
      n_z = lat % n_axial
      allocate(lat % universes(2*n_rings - 1, 2*n_rings - 1, n_z))

      ! Check that number of universes matches size
      n = get_arraysize_integer(node_lat, "universes")
      if (n /= (3*n_rings**2 - 3*n_rings + 1)*n_z) then
        call fatal_error("Number of universes on <universes> does not match &
             &size of lattice " // trim(to_str(lat % id)) // ".")
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      ! Universes in hexagonal lattices are stored in a manner that represents
      ! a skewed coordinate system: (x, alpha) rather than (x, y).  There is
      ! no obvious, direct relationship between the order of universes in the
      ! input and the order that they will be stored in the skewed array so
      ! the following code walks a set of index values across the skewed array
      ! in a manner that matches the input order.  Note that i_x = 0, i_a = 0
      ! corresponds to the center of the hexagonal lattice.

      input_index = 1
      do m = 1, n_z
        ! Initialize lattice indecies.
        i_x = 1
        i_a = n_rings - 1

        ! Map upper triangular region of hexagonal lattice.
        do k = 1, n_rings-1
          ! Walk index to lower-left neighbor of last row start.
          i_x = i_x - 1
          do j = 1, k
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 &temp_int_array(input_index)
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*k
          i_a = i_a + k
        end do

        ! Map middle square region of hexagonal lattice.
        do k = 1, 2*n_rings - 1
          if (mod(k, 2) == 1) then
            ! Walk index to lower-left neighbor of last row start.
            i_x = i_x - 1
          else
            ! Walk index to lower-right neighbor of last row start
            i_x = i_x + 1
            i_a = i_a - 1
          end if
          do j = 1, n_rings - mod(k-1, 2)
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 &temp_int_array(input_index)
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*(n_rings - mod(k-1, 2))
          i_a = i_a + n_rings - mod(k-1, 2)
        end do

        ! Map lower triangular region of hexagonal lattice.
        do k = 1, n_rings-1
          ! Walk index to lower-right neighbor of last row start.
          i_x = i_x + 1
          i_a = i_a - 1
          do j = 1, n_rings - k
            ! Place universe in array.
            lat % universes(i_x + n_rings, i_a + n_rings, m) = &
                 &temp_int_array(input_index)
            ! Walk index to closest non-adjacent right neighbor.
            i_x = i_x + 2
            i_a = i_a - 1
            ! Increment XML array index.
            input_index = input_index + 1
          end do
          ! Return lattice index to start of current row.
          i_x = i_x - 2*(n_rings - k)
          i_a = i_a + n_rings - k
        end do
      end do
      deallocate(temp_int_array)

      ! Read outer universe for area outside lattice.
      lat % outer = NO_OUTER_UNIVERSE
      if (check_for_node(node_lat, "outer")) then
        call get_node_value(node_lat, "outer", lat % outer)
      end if

      ! Check for 'outside' nodes which are no longer supported.
      if (check_for_node(node_lat, "outside")) then
        call fatal_error("The use of 'outside' in lattices is no longer &
             &supported.  Instead, use 'outer' which defines a universe rather &
             &than a material.  The utility openmc/src/utils/update_inputs.py &
             &can be used automatically replace 'outside' with 'outer'.")
      end if

      ! Add lattice to dictionary
      call lattice_dict % add_key(lat % id, n_rlats + i)

      end select
    end do HEX_LATTICES

    ! Close geometry XML file
    call close_xmldoc(doc)

  end subroutine read_geometry_xml

!===============================================================================
! READ_MATERIAL_XML reads data from a materials.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_materials()
    integer :: i, j
    type(DictCharInt) :: library_dict
    type(Library), allocatable :: libraries(:)

    if (run_CE) then
      call read_ce_cross_sections_xml(libraries)
    else
      call read_mg_cross_sections_xml(libraries)
    end if

    ! Creating dictionary that maps the name of the material to the entry
    do i = 1, size(libraries)
      do j = 1, size(libraries(i) % materials)
        call library_dict % add_key(to_lower(libraries(i) % materials(j)), i)
      end do
    end do

    ! Check that 0K nuclides are listed in the cross_sections.xml file
    if (allocated(nuclides_0K)) then
      do i = 1, size(nuclides_0K)
        if (.not. library_dict % has_key(to_lower(nuclides_0K(i) % name_0K))) then
          call fatal_error("Could not find resonant scatterer " &
               // trim(nuclides_0K(i) % name_0K) &
               // " in cross_sections.xml file!")
        end if
      end do
    end if

    ! Parse data from materials.xml
    call read_materials_xml(libraries, library_dict)

    ! Read continuous-energy cross sections
    if (run_CE .and. run_mode /= MODE_PLOTTING) then
      call time_read_xs%start()
      call read_ce_cross_sections(libraries, library_dict)
      call time_read_xs%stop()
    end if

    ! Normalize atom/weight percents
    if (run_mode /= MODE_PLOTTING) call normalize_ao()

    ! Clear dictionary
    call library_dict % clear()
  end subroutine read_materials

  subroutine read_materials_xml(libraries, library_dict)
    type(Library), intent(in) :: libraries(:)
    type(DictCharInt), intent(inout) :: library_dict

    integer :: i              ! loop index for materials
    integer :: j              ! loop index for nuclides
    integer :: k              ! loop index for elements
    integer :: n              ! number of nuclides
    integer :: n_sab          ! number of sab tables for a material
    integer :: n_nuc_ele      ! number of nuclides in an element
    integer :: i_library      ! index in libraries array
    integer :: index_nuclide  ! index in nuclides
    integer :: index_sab      ! index in sab_tables
    real(8) :: val            ! value entered for density
    real(8) :: temp_dble      ! temporary double prec. real
    logical :: file_exists    ! does materials.xml exist?
    logical :: sum_density    ! density is taken to be sum of nuclide densities
    character(20) :: name     ! name of isotope, e.g. 92235.03c
    character(MAX_WORD_LEN) :: units    ! units on density
    character(MAX_LINE_LEN) :: filename ! absolute path to materials.xml
    character(MAX_LINE_LEN) :: temp_str ! temporary string when reading
    type(VectorChar) :: names     ! temporary list of nuclide names
    type(VectorReal) :: densities ! temporary list of nuclide densities
    type(VectorInt)  :: list_iso_lab ! temporary list of isotropic lab scatterers
    type(Material),    pointer :: mat => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mat => null()
    type(Node), pointer :: node_dens => null()
    type(Node), pointer :: node_nuc => null()
    type(Node), pointer :: node_ele => null()
    type(Node), pointer :: node_sab => null()
    type(NodeList), pointer :: node_mat_list => null()
    type(NodeList), pointer :: node_nuc_list => null()
    type(NodeList), pointer :: node_macro_list => null()
    type(NodeList), pointer :: node_ele_list => null()
    type(NodeList), pointer :: node_sab_list => null()

    ! Display output message
    call write_message("Reading materials XML file...", 5)

    ! Check is materials.xml exists
    filename = trim(path_input) // "materials.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Material XML file '" // trim(filename) // "' does not &
           &exist!")
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
        call fatal_error("Must specify id of material in materials XML file")
      end if

      ! Check to make sure 'id' hasn't been used
      if (material_dict % has_key(mat % id)) then
        call fatal_error("Two or more materials use the same unique ID: " &
             // to_str(mat % id))
      end if

      ! Copy material name
      if (check_for_node(node_mat, "name")) then
        call get_node_value(node_mat, "name", mat % name)
      end if

      ! =======================================================================
      ! READ AND PARSE <density> TAG

      ! Get pointer to density element
      if (check_for_node(node_mat, "density")) then
        call get_node_ptr(node_mat, "density", node_dens)
      else
        call fatal_error("Must specify density element in material " &
             // trim(to_str(mat % id)))
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

      else if (units == 'macro') then
        if (check_for_node(node_dens, "value")) then
          ! Copy value
          call get_node_value(node_dens, "value", val)
        else
          val = ONE
        end if

        ! Set density
        mat % density = val

        sum_density = .false.

      else
        ! Copy value
        call get_node_value(node_dens, "value", val)

        ! Check for erroneous density
        sum_density = .false.
        if (val <= ZERO) then
          call fatal_error("Need to specify a positive density on material " &
               // trim(to_str(mat % id)) // ".")
        end if

        ! Adjust material density based on specified units
        select case(to_lower(units))
        case ('g/cc', 'g/cm3')
          mat % density = -val
        case ('kg/m3')
          mat % density = -0.001_8 * val
        case ('atom/b-cm')
          mat % density = val
        case ('atom/cm3', 'atom/cc')
          mat % density = 1.0e-24_8 * val
        case default
          call fatal_error("Unkwown units '" // trim(units) &
               // "' specified on material " // trim(to_str(mat % id)))
        end select
      end if

      ! =======================================================================
      ! READ AND PARSE <nuclide> TAGS

      ! Check to ensure material has at least one nuclide
      if (.not. check_for_node(node_mat, "nuclide") .and. &
           .not. check_for_node(node_mat, "element") .and. &
           .not. check_for_node(node_mat, "macroscopic")) then
        call fatal_error("No macroscopic data, nuclides or natural elements &
                         &specified on material " // trim(to_str(mat % id)))
      end if

      ! Create list of macroscopic x/s based on those specified, just treat
      ! them as nuclides. This is all really a facade so the user thinks they
      ! are entering in macroscopic data but the code treats them the same
      ! as nuclides internally.
      ! Get pointer list of XML <macroscopic>
      call get_node_list(node_mat, "macroscopic", node_macro_list)
      if (run_CE .and. (get_list_size(node_macro_list) > 0)) then
        call fatal_error("Macroscopic can not be used in continuous-energy&
                         & mode!")
      else if (get_list_size(node_macro_list) > 1) then
        call fatal_error("Only one macroscopic object permitted per material, " &
             // trim(to_str(mat % id)))
      else if (get_list_size(node_macro_list) == 1) then

        call get_list_item(node_macro_list, 1, node_nuc)

        ! Check for empty name on nuclide
        if (.not. check_for_node(node_nuc, "name")) then
          call fatal_error("No name specified on macroscopic data in material " &
               // trim(to_str(mat % id)))
        end if

        ! Check for cross section
        if (.not. check_for_node(node_nuc, "xs")) then
          if (default_xs == '') then
            call fatal_error("No cross section specified for macroscopic data &
                 & in material " // trim(to_str(mat % id)))
          else
            name = to_lower(trim(default_xs))
          end if
        end if

        ! store full name
        call get_node_value(node_nuc, "name", temp_str)
        if (check_for_node(node_nuc, "xs")) &
             call get_node_value(node_nuc, "xs", name)
        name = trim(temp_str) // "." // trim(name)
        name = to_lower(name)

        ! save name and density to list
        call names % push_back(name)

        ! Check if no atom/weight percents were specified or if both atom and
        ! weight percents were specified
        if (units == 'macro') then
          call densities % push_back(ONE)
        else
          call fatal_error("Units can only be macro for macroscopic data " &
               // trim(name))
        end if
      else

        ! Get pointer list of XML <nuclide>
        call get_node_list(node_mat, "nuclide", node_nuc_list)

        ! Create list of nuclides based on those specified plus natural elements
        INDIVIDUAL_NUCLIDES: do j = 1, get_list_size(node_nuc_list)
          ! Combine nuclide identifier and cross section and copy into names
          call get_list_item(node_nuc_list, j, node_nuc)

          ! Check for empty name on nuclide
          if (.not. check_for_node(node_nuc, "name")) then
            call fatal_error("No name specified on nuclide in material " &
                 // trim(to_str(mat % id)))
          end if

          ! Check for cross section
          if (.not. check_for_node(node_nuc, "xs")) then
            if (default_xs == '') then
              call fatal_error("No cross section specified for nuclide in &
                   &material " // trim(to_str(mat % id)))
            else
              name = to_lower(trim(default_xs))
            end if
          end if

          ! Check enforced isotropic lab scattering
          if (run_CE) then
            if (check_for_node(node_nuc, "scattering")) then
              call get_node_value(node_nuc, "scattering", temp_str)
              if (adjustl(to_lower(temp_str)) == "iso-in-lab") then
                call list_iso_lab % push_back(1)
              else if (adjustl(to_lower(temp_str)) == "data") then
                call list_iso_lab % push_back(0)
              else
                call fatal_error("Scattering must be isotropic in lab or follow&
                     & the ACE file data")
              end if
            else
              call list_iso_lab % push_back(0)
            end if
          end if

          ! store full name
          call get_node_value(node_nuc, "name", temp_str)
          if (check_for_node(node_nuc, "xs")) &
               call get_node_value(node_nuc, "xs", name)
          name = trim(temp_str) // "." // trim(name)

          ! save name and density to list
          call names % push_back(name)

          ! Check if no atom/weight percents were specified or if both atom and
          ! weight percents were specified
          if (units == 'macro') then
            call densities % push_back(ONE)
          else
            if (.not. check_for_node(node_nuc, "ao") .and. &
                 .not. check_for_node(node_nuc, "wo")) then
              call fatal_error("No atom or weight percent specified for nuclide " &
                   // trim(name))
            elseif (check_for_node(node_nuc, "ao") .and. &
                    check_for_node(node_nuc, "wo")) then
              call fatal_error("Cannot specify both atom and weight percents for a &
                   &nuclide: " // trim(name))
            end if

            ! Copy atom/weight percents
            if (check_for_node(node_nuc, "ao")) then
              call get_node_value(node_nuc, "ao", temp_dble)
              call densities % push_back(temp_dble)
            else
              call get_node_value(node_nuc, "wo", temp_dble)
              call densities % push_back(-temp_dble)
            end if
          end if
        end do INDIVIDUAL_NUCLIDES
      end if

      ! =======================================================================
      ! READ AND PARSE <element> TAGS

      ! Get pointer list of XML <element>
      call get_node_list(node_mat, "element", node_ele_list)

      NATURAL_ELEMENTS: do j = 1, get_list_size(node_ele_list)
        call get_list_item(node_ele_list, j, node_ele)

        ! Check for empty name on natural element
        if (.not. check_for_node(node_ele, "name")) then
          call fatal_error("No name specified on nuclide in material " &
               // trim(to_str(mat % id)))
        end if
        call get_node_value(node_ele, "name", name)

        ! Check for cross section
        if (check_for_node(node_ele, "xs")) then
          call get_node_value(node_ele, "xs", temp_str)
        else
          if (default_xs == '') then
            call fatal_error("No cross section specified for nuclide in &
                 &material " // trim(to_str(mat % id)))
          else
            temp_str = to_lower(trim(default_xs))
          end if
        end if

        ! Check if no atom/weight percents were specified or if both atom and
        ! weight percents were specified
        if (.not. check_for_node(node_ele, "ao") .and. &
             .not. check_for_node(node_ele, "wo")) then
          call fatal_error("No atom or weight percent specified for element " &
               // trim(name))
        elseif (check_for_node(node_ele, "ao") .and. &
                check_for_node(node_ele, "wo")) then
          call fatal_error("Cannot specify both atom and weight percents for &
               &element: " // trim(name))
        end if

        ! Get current number of nuclides
        n_nuc_ele = names % size()

        ! Expand element into naturally-occurring isotopes
        if (check_for_node(node_ele, "ao")) then
          call get_node_value(node_ele, "ao", temp_dble)
          call expand_natural_element(name, temp_str, temp_dble, names, &
               densities)
        else
          call fatal_error("The ability to expand a natural element based on &
               &weight percentage is not yet supported.")
        end if

        ! Compute number of new nuclides from the natural element expansion
        n_nuc_ele = names % size() - n_nuc_ele

        ! Check enforced isotropic lab scattering
        if (run_CE) then
          if (check_for_node(node_ele, "scattering")) then
            call get_node_value(node_ele, "scattering", temp_str)
          else
            temp_str = "data"
          end if

          ! Set ace or iso-in-lab scattering for each nuclide in element
          do k = 1, n_nuc_ele
            if (adjustl(to_lower(temp_str)) == "iso-in-lab") then
              call list_iso_lab % push_back(1)
            else if (adjustl(to_lower(temp_str)) == "data") then
              call list_iso_lab % push_back(0)
            else
              call fatal_error("Scattering must be isotropic in lab or follow&
                   & the ACE file data")
            end if
          end do
        end if

      end do NATURAL_ELEMENTS

      ! ========================================================================
      ! COPY NUCLIDES TO ARRAYS IN MATERIAL

      ! allocate arrays in Material object
      n = names % size()
      mat % n_nuclides = n
      allocate(mat % names(n))
      allocate(mat % nuclide(n))
      allocate(mat % atom_density(n))
      allocate(mat % p0(n))

      ALL_NUCLIDES: do j = 1, mat % n_nuclides
        ! Check that this nuclide is listed in the cross_sections.xml file
        name = trim(names % data(j))
        if (.not. library_dict % has_key(to_lower(name))) then
          call fatal_error("Could not find nuclide " // trim(name) &
               // " in cross_sections data file!")
        end if
        i_library = library_dict % get_key(to_lower(name))

        if (run_CE) then
          ! Check to make sure cross-section is continuous energy neutron table
          if (libraries(i_library) % type /= LIBRARY_NEUTRON) then
            call fatal_error("Cross-section table " // trim(name) &
                 // " is not a continuous-energy neutron table.")
          end if
        end if

        ! If this nuclide hasn't been encountered yet, we need to add its name
        ! and alias to the nuclide_dict
        if (.not. nuclide_dict % has_key(to_lower(name))) then
          index_nuclide    = index_nuclide + 1
          mat % nuclide(j) = index_nuclide

          call nuclide_dict % add_key(to_lower(name), index_nuclide)
        else
          mat % nuclide(j) = nuclide_dict % get_key(to_lower(name))
        end if

        ! Copy name and atom/weight percent
        mat % names(j) = name
        mat % atom_density(j) = densities % data(j)

        ! Cast integer isotropic lab scattering flag to boolean
        if (run_CE) then
          if (list_iso_lab % data(j) == 1) then
            mat % p0(j) = .true.
          else
            mat % p0(j) = .false.
          end if
        end if

      end do ALL_NUCLIDES

      ! Check to make sure either all atom percents or all weight percents are
      ! given
      if (.not. (all(mat % atom_density >= ZERO) .or. &
           all(mat % atom_density <= ZERO))) then
        call fatal_error("Cannot mix atom and weight percents in material " &
             // to_str(mat % id))
      end if

      ! Determine density if it is a sum value
      if (sum_density) mat % density = sum(mat % atom_density)

      ! Clear lists
      call names % clear()
      call densities % clear()
      call list_iso_lab % clear()

      ! =======================================================================
      ! READ AND PARSE <sab> TAG FOR S(a,b) DATA
      if (run_CE) then
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
            if (.not. check_for_node(node_sab, "name") .or. &
                 .not. check_for_node(node_sab, "xs")) then
              call fatal_error("Need to specify <name> and <xs> for S(a,b) &
                   &table.")
            end if
            call get_node_value(node_sab, "name", name)
            call get_node_value(node_sab, "xs", temp_str)
            name = trim(name) // "." // trim(temp_str)
            mat % sab_names(j) = name

            ! Check that this nuclide is listed in the cross_sections.xml file
            if (.not. library_dict % has_key(to_lower(name))) then
              call fatal_error("Could not find S(a,b) table " // trim(name) &
                   // " in cross_sections.xml file!")
            end if

            ! Find index in xs_listing and set the name and alias according to the
            ! listing
            i_library = library_dict % get_key(to_lower(name))

            if (run_CE) then
              ! Check to make sure cross-section is continuous energy neutron table
              if (libraries(i_library) % type /= LIBRARY_THERMAL) then
                call fatal_error("Cross-section table " // trim(name) &
                     // " is not a S(a,b) table.")
              end if
            end if

            ! If this S(a,b) table hasn't been encountered yet, we need to add its
            ! name and alias to the sab_dict
            if (.not. sab_dict % has_key(to_lower(name))) then
              index_sab = index_sab + 1
              mat % i_sab_tables(j) = index_sab
              call sab_dict % add_key(to_lower(name), index_sab)
            else
              mat % i_sab_tables(j) = sab_dict % get_key(to_lower(name))
            end if
          end do
        end if
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

    integer :: d             ! delayed group index
    integer :: i             ! loop over user-specified tallies
    integer :: j             ! loop over words
    integer :: k             ! another loop index
    integer :: l             ! another loop index
    integer :: id            ! user-specified identifier
    integer :: i_mesh        ! index in meshes array
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filters     ! number of filters
    integer :: n_new         ! number of new scores to add based on Yn/Pn tally
    integer :: n_scores      ! number of tot scores after adjusting for Yn/Pn tally
    integer :: n_bins        ! total new bins for this score
    integer :: n_user_trig   ! number of user-specified tally triggers
    integer :: trig_ind      ! index of triggers array for each tally
    integer :: user_trig_ind ! index of user-specified triggers for each tally
    real(8) :: threshold     ! trigger convergence threshold
    integer :: n_order       ! moment order requested
    integer :: n_order_pos   ! oosition of Scattering order in score name string
    integer :: MT            ! user-specified MT for score
    integer :: iarray3(3)    ! temporary integer array
    integer :: imomstr       ! Index of MOMENT_STRS & MOMENT_N_STRS
    logical :: file_exists   ! does tallies.xml file exist?
    real(8) :: rarray3(3)    ! temporary double prec. array
    integer :: Nangle        ! Number of angular bins
    real(8) :: dangle        ! Mu spacing if using automatic allocation
    integer :: iangle        ! Loop counter for building mu filter bins
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    type(DictCharInt) :: trigger_scores
    type(ElemKeyValueCI), pointer :: pair_list
    type(TallyObject),    pointer :: t
    type(RegularMesh), pointer :: m
    type(TallyFilterContainer), allocatable :: filters(:) ! temporary filters
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mesh => null()
    type(Node), pointer :: node_tal => null()
    type(Node), pointer :: node_filt => null()
    type(Node), pointer :: node_trigger=>null()
    type(NodeList), pointer :: node_mesh_list => null()
    type(NodeList), pointer :: node_tal_list => null()
    type(NodeList), pointer :: node_filt_list => null()
    type(NodeList), pointer :: node_trigger_list => null()
    type(ElemKeyValueCI), pointer :: scores
    type(ElemKeyValueCI), pointer :: next

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! Since a tallies.xml file is optional, no error is issued here
      return
    end if

    ! Display output message
    call write_message("Reading tallies XML file...", 5)

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
      if (master) call warning("No tallies present in tallies.xml file!")
    end if

    ! Allocate tally array
    if (n_user_tallies > 0 .and. run_mode /= MODE_PLOTTING) then
      call add_tallies("user", n_user_tallies)
    end if

    ! Check for <assume_separate> setting
    if (check_for_node(doc, "assume_separate")) then
      call get_node_value(doc, "assume_separate", temp_str)
      temp_str = to_lower(temp_str)
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
        call fatal_error("Must specify id for mesh in tally XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (mesh_dict % has_key(m % id)) then
        call fatal_error("Two or more meshes use the same unique ID: " &
             // to_str(m % id))
      end if

      ! Read mesh type
      temp_str = ''
      if (check_for_node(node_mesh, "type")) &
           call get_node_value(node_mesh, "type", temp_str)
      select case (to_lower(temp_str))
      case ('rect', 'rectangle', 'rectangular')
        call warning("Mesh type '" // trim(temp_str) // "' is deprecated. &
             &Please use 'regular' instead.")
        m % type = MESH_REGULAR
      case ('regular')
        m % type = MESH_REGULAR
      case default
        call fatal_error("Invalid mesh type: " // trim(temp_str))
      end select

      ! Determine number of dimensions for mesh
      n = get_arraysize_integer(node_mesh, "dimension")
      if (n /= 2 .and. n /= 3) then
        call fatal_error("Mesh must be two or three dimensions.")
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
        call fatal_error("All entries on the <dimension> element for a tally &
             &mesh must be positive.")
      end if

      ! Read dimensions in each direction
      m % dimension = iarray3(1:n)

      ! Read mesh lower-left corner location
      if (m % n_dimension /= get_arraysize_double(node_mesh, "lower_left")) then
        call fatal_error("Number of entries on <lower_left> must be the same &
             &as the number of entries on <dimension>.")
      end if
      call get_node_array(node_mesh, "lower_left", m % lower_left)

      ! Make sure both upper-right or width were specified
      if (check_for_node(node_mesh, "upper_right") .and. &
           check_for_node(node_mesh, "width")) then
        call fatal_error("Cannot specify both <upper_right> and <width> on a &
             &tally mesh.")
      end if

      ! Make sure either upper-right or width was specified
      if (.not. check_for_node(node_mesh, "upper_right") .and. &
           .not. check_for_node(node_mesh, "width")) then
        call fatal_error("Must specify either <upper_right> and <width> on a &
             &tally mesh.")
      end if

      if (check_for_node(node_mesh, "width")) then
        ! Check to ensure width has same dimensions
        if (get_arraysize_double(node_mesh, "width") /= &
             get_arraysize_double(node_mesh, "lower_left")) then
          call fatal_error("Number of entries on <width> must be the same as &
               &the number of entries on <lower_left>.")
        end if

        ! Check for negative widths
        call get_node_array(node_mesh, "width", rarray3(1:n))
        if (any(rarray3(1:n) < ZERO)) then
          call fatal_error("Cannot have a negative <width> on a tally mesh.")
        end if

        ! Set width and upper right coordinate
        m % width = rarray3(1:n)
        m % upper_right = m % lower_left + m % dimension * m % width

      elseif (check_for_node(node_mesh, "upper_right")) then
        ! Check to ensure width has same dimensions
        if (get_arraysize_double(node_mesh, "upper_right") /= &
             get_arraysize_double(node_mesh, "lower_left")) then
          call fatal_error("Number of entries on <upper_right> must be the &
               &same as the number of entries on <lower_left>.")
        end if

        ! Check that upper-right is above lower-left
        call get_node_array(node_mesh, "upper_right", rarray3(1:n))
        if (any(rarray3(1:n) < m % lower_left)) then
          call fatal_error("The <upper_right> coordinates must be greater than &
               &the <lower_left> coordinates on a tally mesh.")
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

    ! We only need the mesh info for plotting
    if (run_mode == MODE_PLOTTING) return

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
        call fatal_error("Must specify id for tally in tally XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (tally_dict % has_key(t % id)) then
        call fatal_error("Two or more tallies use the same unique ID: " &
             // to_str(t % id))
      end if

      ! Copy tally name
      if (check_for_node(node_tal, "name")) &
           call get_node_value(node_tal, "name", t % name)

      ! =======================================================================
      ! READ DATA FOR FILTERS

      ! Get pointer list to XML <filter> and get number of filters
      call get_node_list(node_tal, "filter", node_filt_list)
      n_filters = get_list_size(node_filt_list)

      ! Allocate filters array
      allocate(t % filters(n_filters))

      READ_FILTERS: do j = 1, n_filters
        ! Get pointer to filter xml node
        call get_list_item(node_filt_list, j, node_filt)

        ! Convert filter type to lower case
        temp_str = ''
        if (check_for_node(node_filt, "type")) &
             call get_node_value(node_filt, "type", temp_str)
        temp_str = to_lower(temp_str)

        ! Determine number of bins
        if (check_for_node(node_filt, "bins")) then
          if (temp_str == 'energy' .or. temp_str == 'energyout' .or. &
               temp_str == 'mu' .or. temp_str == 'polar' .or. &
               temp_str == 'azimuthal') then
            n_words = get_arraysize_double(node_filt, "bins")
          else
            n_words = get_arraysize_integer(node_filt, "bins")
          end if
        else
          call fatal_error("Bins not set in filter on tally " &
               // trim(to_str(t % id)))
        end if

        ! Determine type of filter
        select case (temp_str)

        case ('distribcell')
          ! Allocate and declare the filter type
          allocate(DistribcellFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (DistribcellFilter)
            if (n_words /= 1) call fatal_error("Only one cell can be &
                 &specified per distribcell filter.")
            ! Store bins
            call get_node_value(node_filt, "bins", filt % cell)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_DISTRIBCELL) = j

        case ('cell')
          ! Allocate and declare the filter type
          allocate(CellFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (CellFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % cells(n_words))
            call get_node_array(node_filt, "bins", filt % cells)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_CELL) = j

        case ('cellborn')
          ! Allocate and declare the filter type
          allocate(CellbornFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (CellbornFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % cells(n_words))
            call get_node_array(node_filt, "bins", filt % cells)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_CELLBORN) = j

        case ('material')
          ! Allocate and declare the filter type
          allocate(MaterialFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (MaterialFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % materials(n_words))
            call get_node_array(node_filt, "bins", filt % materials)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_MATERIAL) = j

        case ('universe')
          ! Allocate and declare the filter type
          allocate(UniverseFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (UniverseFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % universes(n_words))
            call get_node_array(node_filt, "bins", filt % universes)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_UNIVERSE) = j

        case ('surface')
          call fatal_error("Surface filter is not yet supported!")
          ! Allocate and declare the filter type
          allocate(SurfaceFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (SurfaceFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % surfaces(n_words))
            call get_node_array(node_filt, "bins", filt % surfaces)
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_SURFACE) = j

        case ('mesh')
          ! Allocate and declare the filter type
          allocate(MeshFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (MeshFilter)
            if (n_words /= 1) call fatal_error("Only one mesh can be &
                 &specified per mesh filter.")

            ! Determine id of mesh
            call get_node_value(node_filt, "bins", id)

            ! Get pointer to mesh
            if (mesh_dict % has_key(id)) then
              i_mesh = mesh_dict % get_key(id)
              m => meshes(i_mesh)
            else
              call fatal_error("Could not find mesh " // trim(to_str(id)) &
                   // " specified on tally " // trim(to_str(t % id)))
            end if

            ! Determine number of bins
            filt % n_bins = product(m % dimension)

            ! Store the index of the mesh
            filt % mesh = i_mesh
          end select

          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_MESH) = j

        case ('energy')

          ! Allocate and declare the filter type
          allocate(EnergyFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (EnergyFilter)
            ! Allocate and store bins
            filt % n_bins = n_words - 1
            allocate(filt % bins(n_words))
            call get_node_array(node_filt, "bins", filt % bins)

            ! We can save tallying time if we know that the tally bins match
            ! the energy group structure.  In that case, the matching bin
            ! index is simply the group (after flipping for the different
            ! ordering of the library and tallying systems).
            if (.not. run_CE) then
              if (n_words == energy_groups + 1) then
                if (all(filt % bins == energy_bins(energy_groups + 1:1:-1))) &
                     then
                  filt % matches_transport_groups = .true.
                end if
              end if
            end if
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_ENERGYIN) = j

        case ('energyout')
          ! Allocate and declare the filter type
          allocate(EnergyoutFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (EnergyoutFilter)
            ! Allocate and store bins
            filt % n_bins = n_words - 1
            allocate(filt % bins(n_words))
            call get_node_array(node_filt, "bins", filt % bins)

            ! We can save tallying time if we know that the tally bins match
            ! the energy group structure.  In that case, the matching bin
            ! index is simply the group (after flipping for the different
            ! ordering of the library and tallying systems).
            if (.not. run_CE) then
              if (n_words == energy_groups + 1) then
                if (all(filt % bins == energy_bins(energy_groups + 1:1:-1))) &
                     then
                  filt % matches_transport_groups = .true.
                end if
              end if
            end if
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_ENERGYOUT) = j

          ! Set to analog estimator
          t % estimator = ESTIMATOR_ANALOG

        case ('delayedgroup')
          ! Check to see if running in MG mode, because if so, the current
          ! system isnt set up yet to support delayed group data and thus
          ! these tallies
          if (.not. run_CE) then
            call fatal_error("delayedgroup filter on tally " &
                             // trim(to_str(t % id)) // " not yet supported&
                             & for multi-group mode.")
          end if

          ! Allocate and declare the filter type
          allocate(DelayedGroupFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (DelayedGroupFilter)
            ! Allocate and store bins
            filt % n_bins = n_words
            allocate(filt % groups(n_words))
            call get_node_array(node_filt, "bins", filt % groups)

            ! Check that bins are all are between 1 and MAX_DELAYED_GROUPS
            do d = 1, n_words
              if (filt % groups(d) < 1 .or. &
                   filt % groups(d) > MAX_DELAYED_GROUPS) then
                call fatal_error("Encountered delayedgroup bin with index " &
                     // trim(to_str(filt % groups(d))) // " that is outside &
                     &the range of 1 to MAX_DELAYED_GROUPS ( " &
                     // trim(to_str(MAX_DELAYED_GROUPS)) // ")")
              end if
            end do
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_DELAYEDGROUP) = j

        case ('mu')
          ! Allocate and declare the filter type
          allocate(MuFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (MuFilter)
            ! Allocate and store bins
            filt % n_bins = n_words - 1
            allocate(filt % bins(n_words))
            call get_node_array(node_filt, "bins", filt % bins)

            ! Allow a user to input a lone number which will mean that you
            ! subdivide [-1,1] evenly with the input being the number of bins
            if (n_words == 1) then
              Nangle = int(filt % bins(1))
              if (Nangle > 1) then
                filt % n_bins = Nangle
                dangle = TWO / real(Nangle,8)
                deallocate(filt % bins)
                allocate(filt % bins(Nangle + 1))
                do iangle = 1, Nangle
                  filt % bins(iangle) = -ONE + (iangle - 1) * dangle
                end do
                filt % bins(Nangle + 1) = ONE
              else
                call fatal_error("Number of bins for mu filter must be&
                     & greater than 1 on tally " &
                     // trim(to_str(t % id)) // ".")
              end if
            end if
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_MU) = j

          ! Set to analog estimator
          t % estimator = ESTIMATOR_ANALOG

        case ('polar')
          ! Allocate and declare the filter type
          allocate(PolarFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (PolarFilter)
            ! Allocate and store bins
            filt % n_bins = n_words - 1
            allocate(filt % bins(n_words))
            call get_node_array(node_filt, "bins", filt % bins)

            ! Allow a user to input a lone number which will mean that you
            ! subdivide [0,pi] evenly with the input being the number of bins
            if (n_words == 1) then
              Nangle = int(filt % bins(1))
              if (Nangle > 1) then
                filt % n_bins = Nangle
                dangle = PI / real(Nangle,8)
                deallocate(filt % bins)
                allocate(filt % bins(Nangle + 1))
                do iangle = 1, Nangle
                  filt % bins(iangle) = (iangle - 1) * dangle
                end do
                filt % bins(Nangle + 1) = PI
              else
                call fatal_error("Number of bins for mu filter must be&
                     & greater than 1 on tally " &
                     // trim(to_str(t % id)) // ".")
              end if
            end if
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_POLAR) = j

        case ('azimuthal')
          ! Allocate and declare the filter type
          allocate(AzimuthalFilter::t % filters(j) % obj)
          select type (filt => t % filters(j) % obj)
          type is (AzimuthalFilter)
            ! Allocate and store bins
            filt % n_bins = n_words - 1
            allocate(filt % bins(n_words))
            call get_node_array(node_filt, "bins", filt % bins)

            ! Allow a user to input a lone number which will mean that you
            ! subdivide [-pi,pi) evenly with the input being the number of
            ! bins
            if (n_words == 1) then
              Nangle = int(filt % bins(1))
              if (Nangle > 1) then
                filt % n_bins = Nangle
                dangle = TWO * PI / real(Nangle,8)
                deallocate(filt % bins)
                allocate(filt % bins(Nangle + 1))
                do iangle = 1, Nangle
                  filt % bins(iangle) = -PI + (iangle - 1) * dangle
                end do
                filt % bins(Nangle + 1) = PI
              else
                call fatal_error("Number of bins for mu filter must be&
                     & greater than 1 on tally " &
                     // trim(to_str(t % id)) // ".")
              end if
            end if
          end select
          ! Set the filter index in the tally find_filter array
          t % find_filter(FILTER_AZIMUTHAL) = j

        case default
          ! Specified tally filter is invalid, raise error
          call fatal_error("Unknown filter type '" &
               // trim(temp_str) // "' on tally " &
               // trim(to_str(t % id)) // ".")

        end select

      end do READ_FILTERS

      ! Check that both cell and surface weren't specified
      if (t % find_filter(FILTER_CELL) > 0 .and. &
           t % find_filter(FILTER_SURFACE) > 0) then
        call fatal_error("Cannot specify both cell and surface filters for &
             &tally " // trim(to_str(t % id)))
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

              ! Check if a delayedgroup filter is present for this tally
              do l = 1, size(t % filters)
                select type(filt => t % filters(l) % obj)
                type is (DelayedGroupFilter)
                  call warning("A delayedgroup filter was used on a total &
                       &nuclide tally. Cross section libraries are not &
                       &guaranteed to have the same delayed group structure &
                       &across all isotopes. In particular, ENDF/B-VII.1 does &
                       &not have a consistent delayed group structure across &
                       &all isotopes while the JEFF 3.1.1 library has the same &
                       &delayed group structure across all isotopes. Use with &
                       &caution!")
                end select
              end do

              t % nuclide_bins(j) = -1
              cycle
            end if

            ! If a specific nuclide was specified
            word = to_lower(sarray(j))

            ! Search through nuclides
            pair_list => nuclide_dict % keys()
            do while (associated(pair_list))
              if (starts_with(pair_list % key, word)) then
                word = pair_list % key(1:150)
                exit
              end if

              ! Advance to next
              pair_list => pair_list % next
            end do

            ! Check if no nuclide was found
            if (.not. associated(pair_list)) then
              call fatal_error("Could not find the nuclide " &
                   // trim(word) // " specified in tally " &
                   // trim(to_str(t % id)) // " in any material.")
            end if
            deallocate(pair_list)

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

        ! Check if a delayedgroup filter is present for this tally
        do l = 1, size(t % filters)
          select type(filt => t % filters(l) % obj)
          type is (DelayedGroupFilter)
            call warning("A delayedgroup filter was used on a total nuclide &
                 &tally. Cross section libraries are not guaranteed to have the&
                 & same delayed group structure across all isotopes. In &
                 &particular, ENDF/B-VII.1 does not have a consistent delayed &
                 &group structure across all isotopes while the JEFF 3.1.1 &
                 &library has the same delayed group structure across all &
                 &isotopes. Use with caution!")
          end select
        end do
      end if

      ! =======================================================================
      ! READ DATA FOR SCORES

      if (check_for_node(node_tal, "scores")) then
        n_words = get_arraysize_string(node_tal, "scores")
        allocate(sarray(n_words))
        call get_node_array(node_tal, "scores", sarray)

        ! Before we can allocate storage for scores, we must determine the
        ! number of additional scores required due to the moment scores
        ! (i.e., scatter-p#, flux-y#)
        n_new = 0
        do j = 1, n_words
          sarray(j) = to_lower(sarray(j))
          ! Find if scores(j) is of the form 'moment-p' or 'moment-y' present in
          ! MOMENT_STRS(:)
          ! If so, check the order, store if OK, then reset the number to 'n'
          score_name = trim(sarray(j))

          ! Append the score to the list of possible trigger scores
          if (trigger_on) call trigger_scores % add_key(trim(score_name), j)

          do imomstr = 1, size(MOMENT_STRS)
            if (starts_with(score_name,trim(MOMENT_STRS(imomstr)))) then
              n_order_pos = scan(score_name,'0123456789')
              n_order = int(str_to_int( &
                   score_name(n_order_pos:(len_trim(score_name)))),4)
              if (n_order > MAX_ANG_ORDER) then
                ! User requested too many orders; throw a warning and set to the
                ! maximum order.
                ! The above scheme will essentially take the absolute value
                if (master) call warning("Invalid scattering order of " &
                     // trim(to_str(n_order)) // " requested. Setting to the &
                     &maximum permissible value, " &
                     // trim(to_str(MAX_ANG_ORDER)))
                n_order = MAX_ANG_ORDER
                sarray(j) = trim(MOMENT_STRS(imomstr)) &
                     // trim(to_str(MAX_ANG_ORDER))
              end if
              ! Find total number of bins for this case
              if (imomstr >= YN_LOC) then
                n_bins = (n_order + 1)**2
              else
                n_bins = n_order + 1
              end if
              ! We subtract one since n_words already included
              n_new = n_new + n_bins - 1
              exit
            end if
          end do
        end do
        n_scores = n_words + n_new

        ! Allocate score storage accordingly
        allocate(t % score_bins(n_scores))
        allocate(t % moment_order(n_scores))
        t % moment_order = 0
        j = 0
        do l = 1, n_words
          j = j + 1
          ! Get the input string in scores(l) but if score is one of the moment
          ! scores then strip off the n and store it as an integer to be used
          ! later. Then perform the select case on this modified (number
          ! removed) string
          n_order = -1
          score_name = sarray(l)
          do imomstr = 1, size(MOMENT_STRS)
            if (starts_with(score_name,trim(MOMENT_STRS(imomstr)))) then
              n_order_pos = scan(score_name,'0123456789')
              n_order = int(str_to_int( &
                   score_name(n_order_pos:(len_trim(score_name)))),4)
              if (n_order > MAX_ANG_ORDER) then
                ! User requested too many orders; throw a warning and set to the
                ! maximum order.
                ! The above scheme will essentially take the absolute value
                n_order = MAX_ANG_ORDER
              end if
              score_name = trim(MOMENT_STRS(imomstr)) // "n"
              ! Find total number of bins for this case
              if (imomstr >= YN_LOC) then
                n_bins = (n_order + 1)**2
              else
                n_bins = n_order + 1
              end if
              exit
            end if
          end do
          ! Now check the Moment_N_Strs, but only if we werent successful above
          if (imomstr > size(MOMENT_STRS)) then
            do imomstr = 1, size(MOMENT_N_STRS)
              if (starts_with(score_name,trim(MOMENT_N_STRS(imomstr)))) then
                n_order_pos = scan(score_name,'0123456789')
                n_order = int(str_to_int( &
                     score_name(n_order_pos:(len_trim(score_name)))),4)
                if (n_order > MAX_ANG_ORDER) then
                  ! User requested too many orders; throw a warning and set to the
                  ! maximum order.
                  ! The above scheme will essentially take the absolute value
                  if (master) call warning("Invalid scattering order of " &
                       // trim(to_str(n_order)) // " requested. Setting to &
                       &the maximum permissible value, " &
                       // trim(to_str(MAX_ANG_ORDER)))
                  n_order = MAX_ANG_ORDER
                end if
                score_name = trim(MOMENT_N_STRS(imomstr)) // "n"
                exit
              end if
            end do
          end if

          ! Check if delayed group filter is used with any score besides
          ! delayed-nu-fission or decay-rate
          if ((score_name /= 'delayed-nu-fission' .and. &
               score_name /= 'decay-rate') .and. &
               t % find_filter(FILTER_DELAYEDGROUP) > 0) then
            call fatal_error("Cannot tally " // trim(score_name) // " with a &
                 &delayedgroup filter.")
          end if

          ! Check to see if the mu filter is applied and if that makes sense.
          if ((.not. starts_with(score_name,'scatter')) .and. &
               (.not. starts_with(score_name,'nu-scatter'))) then
            if (t % find_filter(FILTER_MU) > 0) then
              call fatal_error("Cannot tally " // trim(score_name) //" with a &
                               &change of angle (mu) filter.")
            end if
          ! Also check to see if this is a legendre expansion or not.
          ! If so, we can accept this score and filter combo for p0, but not
          ! elsewhere.
          else if (n_order > 0) then
            if (t % find_filter(FILTER_MU) > 0) then
              call fatal_error("Cannot tally " // trim(score_name) //" with a &
                               &change of angle (mu) filter unless order is 0.")
            end if
          end if

          select case (trim(score_name))
          case ('flux')
            ! Prohibit user from tallying flux for an individual nuclide
            if (.not. (t % n_nuclide_bins == 1 .and. &
                 t % nuclide_bins(1) == -1)) then
              call fatal_error("Cannot tally flux for an individual nuclide.")
            end if

            t % score_bins(j) = SCORE_FLUX
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally flux with an outgoing energy &
                   &filter.")
            end if
          case ('flux-yn')
            ! Prohibit user from tallying flux for an individual nuclide
            if (.not. (t % n_nuclide_bins == 1 .and. &
                 t % nuclide_bins(1) == -1)) then
              call fatal_error("Cannot tally flux for an individual nuclide.")
            end if

            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally flux with an outgoing energy &
                   &filter.")
            end if

            t % score_bins(j : j + n_bins - 1) = SCORE_FLUX_YN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins  - 1

          case ('total', '(n,total)')
            t % score_bins(j) = SCORE_TOTAL
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally total reaction rate with an &
                   &outgoing energy filter.")
            end if

          case ('total-yn')
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally total reaction rate with an &
                   &outgoing energy filter.")
            end if

            t % score_bins(j : j + n_bins - 1) = SCORE_TOTAL_YN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins - 1

          case ('scatter')
            t % score_bins(j) = SCORE_SCATTER

          case ('nu-scatter')
            t % score_bins(j) = SCORE_NU_SCATTER

            ! Set tally estimator to analog for CE mode
            ! (MG mode has all data available without a collision being
            ! necessary)
            if (run_CE) then
              t % estimator = ESTIMATOR_ANALOG
            end if

          case ('scatter-n')
            t % score_bins(j) = SCORE_SCATTER_N
            t % moment_order(j) = n_order
            t % estimator = ESTIMATOR_ANALOG

          case ('nu-scatter-n')
            t % score_bins(j) = SCORE_NU_SCATTER_N
            t % moment_order(j) = n_order
            t % estimator = ESTIMATOR_ANALOG

          case ('scatter-pn')
            t % estimator = ESTIMATOR_ANALOG
            ! Setup P0:Pn
            t % score_bins(j : j + n_bins - 1) = SCORE_SCATTER_PN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins - 1

          case ('nu-scatter-pn')
            t % estimator = ESTIMATOR_ANALOG
            ! Setup P0:Pn
            t % score_bins(j : j + n_bins - 1) = SCORE_NU_SCATTER_PN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins - 1

          case ('scatter-yn')
            t % estimator = ESTIMATOR_ANALOG
            ! Setup P0:Pn
            t % score_bins(j : j + n_bins - 1) = SCORE_SCATTER_YN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins - 1

          case ('nu-scatter-yn')
            t % estimator = ESTIMATOR_ANALOG
            ! Setup P0:Pn
            t % score_bins(j : j + n_bins - 1) = SCORE_NU_SCATTER_YN
            t % moment_order(j : j + n_bins - 1) = n_order
            j = j + n_bins - 1

          case('transport')
            call fatal_error("Transport score no longer supported for tallies, &
                 &please remove")

          case ('n1n')
            call fatal_error("n1n score no longer supported for tallies, &
                 &please remove")
          case ('n2n', '(n,2n)')
            t % score_bins(j) = N_2N

          case ('n3n', '(n,3n)')
            t % score_bins(j) = N_3N

          case ('n4n', '(n,4n)')
            t % score_bins(j) = N_4N

          case ('absorption')
            t % score_bins(j) = SCORE_ABSORPTION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally absorption rate with an outgoing &
                   &energy filter.")
            end if
          case ('fission', '18')
            t % score_bins(j) = SCORE_FISSION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally fission rate with an outgoing &
                   &energy filter.")
            end if
          case ('nu-fission')
            t % score_bins(j) = SCORE_NU_FISSION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if
          case ('decay-rate')
            t % score_bins(j) = SCORE_DECAY_RATE
            t % estimator = ESTIMATOR_ANALOG
          case ('delayed-nu-fission')
            t % score_bins(j) = SCORE_DELAYED_NU_FISSION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if
          case ('prompt-nu-fission')
            t % score_bins(j) = SCORE_PROMPT_NU_FISSION
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if

            ! Disallow for MG mode since data not present
            if (.not. run_CE) then
              call fatal_error("Cannot tally delayed nu-fission rate in &
                               &multi-group mode")
            end if
          case ('kappa-fission')
            t % score_bins(j) = SCORE_KAPPA_FISSION
          case ('inverse-velocity')
            t % score_bins(j) = SCORE_INVERSE_VELOCITY
          case ('fission-q-prompt')
            t % score_bins(j) = SCORE_FISS_Q_PROMPT
          case ('fission-q-recoverable')
            t % score_bins(j) = SCORE_FISS_Q_RECOV
          case ('current')
            t % score_bins(j) = SCORE_CURRENT
            t % type = TALLY_SURFACE_CURRENT

            ! Check to make sure that current is the only desired response
            ! for this tally
            if (n_words > 1) then
              call fatal_error("Cannot tally other scores in the &
                   &same tally as surface currents")
            end if

            ! Get index of mesh filter
            k = t % find_filter(FILTER_MESH)

            ! Check to make sure mesh filter was specified
            if (k == 0) then
              call fatal_error("Cannot tally surface current without a mesh &
                   &filter.")
            end if

            ! Copy filters to temporary array
            allocate(filters(size(t % filters) + 1))
            filters(1:size(t % filters)) = t % filters

            ! Move allocation back -- filters becomes deallocated during
            ! this call
            call move_alloc(FROM=filters, TO=t%filters)

            ! Add surface filter
            n_filters = size(t % filters)
            allocate(SurfaceFilter :: t % filters(n_filters) % obj)
            select type (filt => t % filters(size(t % filters)) % obj)
            type is (SurfaceFilter)
              filt % n_bins = 4 * m % n_dimension
              allocate(filt % surfaces(4 * m % n_dimension))
              if (m % n_dimension == 2) then
                filt % surfaces = (/ OUT_LEFT, OUT_RIGHT, OUT_BACK, OUT_FRONT, &
                     IN_LEFT, IN_RIGHT, IN_BACK, IN_FRONT /)
              elseif (m % n_dimension == 3) then
                filt % surfaces = (/ OUT_LEFT, OUT_RIGHT, OUT_BACK, OUT_FRONT, &
                     OUT_BOTTOM, OUT_TOP, IN_LEFT, IN_RIGHT, IN_BACK, &
                     IN_FRONT, IN_BOTTOM, IN_TOP /)
              end if
            end select
            t % find_filter(FILTER_SURFACE) = size(t % filters)

          case ('events')
            t % score_bins(j) = SCORE_EVENTS

          case ('elastic', '(n,elastic)')
            t % score_bins(j) = ELASTIC
          case ('(n,2nd)')
            t % score_bins(j) = N_2ND
          case ('(n,na)')
            t % score_bins(j) = N_2NA
          case ('(n,n3a)')
            t % score_bins(j) = N_N3A
          case ('(n,2na)')
            t % score_bins(j) = N_2NA
          case ('(n,3na)')
            t % score_bins(j) = N_3NA
          case ('(n,np)')
            t % score_bins(j) = N_NP
          case ('(n,n2a)')
            t % score_bins(j) = N_N2A
          case ('(n,2n2a)')
            t % score_bins(j) = N_2N2A
          case ('(n,nd)')
            t % score_bins(j) = N_ND
          case ('(n,nt)')
            t % score_bins(j) = N_NT
          case ('(n,nHe-3)')
            t % score_bins(j) = N_N3HE
          case ('(n,nd2a)')
            t % score_bins(j) = N_ND2A
          case ('(n,nt2a)')
            t % score_bins(j) = N_NT2A
          case ('(n,3nf)')
            t % score_bins(j) = N_3NF
          case ('(n,2np)')
            t % score_bins(j) = N_2NP
          case ('(n,3np)')
            t % score_bins(j) = N_3NP
          case ('(n,n2p)')
            t % score_bins(j) = N_N2P
          case ('(n,npa)')
            t % score_bins(j) = N_NPA
          case ('(n,n1)')
            t % score_bins(j) = N_N1
          case ('(n,nc)')
            t % score_bins(j) = N_NC
          case ('(n,gamma)')
            t % score_bins(j) = N_GAMMA
          case ('(n,p)')
            t % score_bins(j) = N_P
          case ('(n,d)')
            t % score_bins(j) = N_D
          case ('(n,t)')
            t % score_bins(j) = N_T
          case ('(n,3He)')
            t % score_bins(j) = N_3HE
          case ('(n,a)')
            t % score_bins(j) = N_A
          case ('(n,2a)')
            t % score_bins(j) = N_2A
          case ('(n,3a)')
            t % score_bins(j) = N_3A
          case ('(n,2p)')
            t % score_bins(j) = N_2P
          case ('(n,pa)')
            t % score_bins(j) = N_PA
          case ('(n,t2a)')
            t % score_bins(j) = N_T2A
          case ('(n,d2a)')
            t % score_bins(j) = N_D2A
          case ('(n,pd)')
            t % score_bins(j) = N_PD
          case ('(n,pt)')
            t % score_bins(j) = N_PT
          case ('(n,da)')
            t % score_bins(j) = N_DA

          case default
            ! Assume that user has specified an MT number
            MT = int(str_to_int(score_name))

            if (MT /= ERROR_INT) then
              ! Specified score was an integer
              if (MT > 1) then
                t % score_bins(j) = MT
              else
                call fatal_error("Invalid MT on <scores>: " &
                     // trim(sarray(l)))
              end if

            else
              ! Specified score was not an integer
              call fatal_error("Unknown scoring function: " &
                   // trim(sarray(l)))
            end if

          end select

          ! Do a check at the end (instead of for every case) to make sure
          ! the tallies are compatible with MG mode where we have less detailed
          ! nuclear data
          if (.not. run_CE .and. t % score_bins(j) > 0) then
            call fatal_error("Cannot tally " // trim(score_name) // &
                             " reaction rate in multi-group mode")
          end if
        end do

        t % n_score_bins = n_scores
        t % n_user_score_bins = n_words

        ! Deallocate temporary string array of scores
        deallocate(sarray)

        ! Check that no duplicate scores exist
        j = 1
        do while (j < n_scores)
          ! Determine number of bins for scores with expansions
          n_order = t % moment_order(j)
          select case (t % score_bins(j))
          case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
            n_bins = n_order + 1
          case (SCORE_FLUX_YN, SCORE_TOTAL_YN, SCORE_SCATTER_YN, &
               SCORE_NU_SCATTER_YN)
            n_bins = (n_order + 1)**2
          case default
            n_bins = 1
          end select

          do k = j + n_bins, n_scores
            if (t % score_bins(j) == t % score_bins(k) .and. &
                 t % moment_order(j) == t % moment_order(k)) then
              call fatal_error("Duplicate score of type '" // trim(&
                   reaction_name(t % score_bins(j))) // "' found in tally " &
                   // trim(to_str(t % id)))
            end if
          end do
          j = j + n_bins
        end do
      else
        call fatal_error("No <scores> specified on tally " &
             // trim(to_str(t % id)) // ".")
      end if

      ! If settings.xml trigger is turned on, create tally triggers
      if (trigger_on) then

        ! Get list of trigger nodes for this tally
        call get_node_list(node_tal, "trigger", node_trigger_list)

        ! Initialize the number of triggers
        n_user_trig = get_list_size(node_trigger_list)

        ! Count the number of triggers needed for all scores including "all"
        t % n_triggers = 0
        COUNT_TRIGGERS: do user_trig_ind = 1, n_user_trig

          ! Get pointer to trigger node
          call get_list_item(node_trigger_list, user_trig_ind, node_trigger)

          ! Get scores for this trigger
          if (check_for_node(node_trigger, "scores")) then
            n_words = get_arraysize_string(node_trigger, "scores")
            allocate(sarray(n_words))
            call get_node_array(node_trigger, "scores", sarray)
          else
            n_words = 1
            allocate(sarray(n_words))
            sarray(1) = "all"
          end if

          ! Count the number of scores for this trigger
          do j = 1, n_words
            score_name = trim(to_lower(sarray(j)))

            if (score_name == "all") then
              scores => trigger_scores % keys()

              do while (associated(scores))
                next => scores % next
                deallocate(scores)
                scores => next
                t % n_triggers = t % n_triggers + 1
              end do

            else
              t % n_triggers = t % n_triggers + 1
            end if

          end do

          deallocate(sarray)

        end do COUNT_TRIGGERS

        ! Allocate array of triggers for this tally
        if (t % n_triggers > 0) then
          allocate(t % triggers(t % n_triggers))
        end if

        ! Initialize overall trigger index for this tally to zero
        trig_ind = 1

        ! Create triggers for all scores specified on each trigger
        TRIGGER_LOOP: do user_trig_ind = 1, n_user_trig

          ! Get pointer to trigger node
          call get_list_item(node_trigger_list, user_trig_ind, node_trigger)

          ! Get the trigger type - "variance", "std_dev" or "rel_err"
          if (check_for_node(node_trigger, "type")) then
            call get_node_value(node_trigger, "type", temp_str)
            temp_str = to_lower(temp_str)
          else
            call fatal_error("Must specify trigger type for tally " // &
                 trim(to_str(t % id)) // " in tally XML file.")
          end if

          ! Get the convergence threshold for the trigger
          if (check_for_node(node_trigger, "threshold")) then
            call get_node_value(node_trigger, "threshold", threshold)
          else
            call fatal_error("Must specify trigger threshold for tally " // &
                 trim(to_str(t % id)) // " in tally XML file.")
          end if

          ! Get list scores for this trigger
          if (check_for_node(node_trigger, "scores")) then
            n_words = get_arraysize_string(node_trigger, "scores")
            allocate(sarray(n_words))
            call get_node_array(node_trigger, "scores", sarray)
          else
            n_words = 1
            allocate(sarray(n_words))
            sarray(1) = "all"
          end if

          ! Create a trigger for each score
          SCORE_LOOP: do j = 1, n_words
            score_name = trim(to_lower(sarray(j)))

            ! Expand "all" to include TriggerObjects for each score in tally
            if (score_name == "all") then
              scores => trigger_scores % keys()

              ! Loop over all tally scores
              do while (associated(scores))
                score_name = trim(scores % key)

                ! Store the score name and index in the trigger
                t % triggers(trig_ind) % score_name = trim(score_name)
                t % triggers(trig_ind) % score_index = &
                     trigger_scores % get_key(trim(score_name))

                ! Set the trigger convergence threshold type
                select case (temp_str)
                case ('std_dev')
                  t % triggers(trig_ind) % type = STANDARD_DEVIATION
                case ('variance')
                  t % triggers(trig_ind) % type = VARIANCE
                case ('rel_err')
                  t % triggers(trig_ind) % type = RELATIVE_ERROR
                case default
                  call fatal_error("Unknown trigger type " // &
                       trim(temp_str) // " in tally " // trim(to_str(t % id)))
                end select

                ! Store the trigger convergence threshold
                t % triggers(trig_ind) % threshold = threshold

                ! Move to next score
                next => scores % next
                deallocate(scores)
                scores => next

                ! Increment the overall trigger index
                trig_ind = trig_ind + 1
              end do

            ! Scores other than the "all" placeholder
            else

              ! Store the score name and index
              t % triggers(trig_ind) % score_name = trim(score_name)
              t % triggers(trig_ind) % score_index = &
                   trigger_scores % get_key(trim(score_name))

              ! Check if an invalid score was set for the trigger
              if (t % triggers(trig_ind) % score_index == 0) then
                call fatal_error("The trigger score " // trim(score_name) // &
                     " is not set for tally " // trim(to_str(t % id)))
              end if

              ! Store the trigger convergence threshold
              t % triggers(trig_ind) % threshold = threshold

              ! Set the trigger convergence threshold type
              select case (temp_str)
              case ('std_dev')
                t % triggers(trig_ind) % type = STANDARD_DEVIATION
              case ('variance')
                t % triggers(trig_ind) % type = VARIANCE
              case ('rel_err')
                t % triggers(trig_ind) % type = RELATIVE_ERROR
              case default
                call fatal_error("Unknown trigger type " // trim(temp_str) // &
                     " in tally " // trim(to_str(t % id)))
              end select

              ! Increment the overall trigger index
              trig_ind = trig_ind + 1
            end if
          end do SCORE_LOOP

          ! Deallocate the list of tally scores used to create triggers
          deallocate(sarray)
        end do TRIGGER_LOOP

        ! Deallocate dictionary of scores/indices used to populate triggers
        call trigger_scores % clear()
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
            call fatal_error("Cannot use track-length estimator for tally " &
                 // to_str(t % id))
          end if

          ! Set estimator to track-length estimator
          t % estimator = ESTIMATOR_TRACKLENGTH

        case ('collision')
          ! If the estimator was set to an analog estimator, this means the
          ! tally needs post-collision information
          if (t % estimator == ESTIMATOR_ANALOG) then
            call fatal_error("Cannot use collision estimator for tally " &
                 // to_str(t % id))
          end if

          ! Set estimator to collision estimator
          t % estimator = ESTIMATOR_COLLISION

        case default
          call fatal_error("Invalid estimator '" // trim(temp_str) &
               // "' on tally " // to_str(t % id))
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

    integer :: i, j
    integer :: n_cols, col_id, n_comp, n_masks, n_meshlines
    integer :: meshid
    integer :: i_mesh
    integer, allocatable :: iarray(:)
    logical :: file_exists              ! does plots.xml file exist?
    character(MAX_LINE_LEN) :: filename ! absolute path to plots.xml
    character(MAX_LINE_LEN) :: temp_str
    character(MAX_WORD_LEN) :: meshtype
    type(ObjectPlot), pointer :: pl => null()
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_plot => null()
    type(Node), pointer :: node_col => null()
    type(Node), pointer :: node_mask => null()
    type(Node), pointer :: node_meshlines => null()
    type(NodeList), pointer :: node_plot_list => null()
    type(NodeList), pointer :: node_col_list => null()
    type(NodeList), pointer :: node_mask_list => null()
    type(NodeList), pointer :: node_meshline_list => null()

    ! Check if plots.xml exists
    filename = trim(path_input) // "plots.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Plots XML file '" // trim(filename) &
           // "' does not exist!")
    end if

    ! Display output message
    call write_message("Reading plot XML file...", 5)

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
        call fatal_error("Must specify plot id in plots XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (plot_dict % has_key(pl % id)) then
        call fatal_error("Two or more plots use the same unique ID: " &
             // to_str(pl % id))
      end if

      ! Copy plot type
      temp_str = 'slice'
      if (check_for_node(node_plot, "type")) &
           call get_node_value(node_plot, "type", temp_str)
      temp_str = to_lower(temp_str)
      select case (trim(temp_str))
      case ("slice")
        pl % type = PLOT_TYPE_SLICE
      case ("voxel")
        pl % type = PLOT_TYPE_VOXEL
      case default
        call fatal_error("Unsupported plot type '" // trim(temp_str) &
             // "' in plot " // trim(to_str(pl % id)))
      end select

      ! Set output file path
      filename = trim(to_str(pl % id)) // "_plot"
      if (check_for_node(node_plot, "filename")) &
           call get_node_value(node_plot, "filename", filename)
      select case (pl % type)
      case (PLOT_TYPE_SLICE)
        pl % path_plot = trim(path_input) // trim(filename) // ".ppm"
      case (PLOT_TYPE_VOXEL)
        pl % path_plot = trim(path_input) // trim(filename) // ".voxel"
      end select

      ! Copy plot pixel size
      if (pl % type == PLOT_TYPE_SLICE) then
        if (get_arraysize_integer(node_plot, "pixels") == 2) then
          call get_node_array(node_plot, "pixels", pl % pixels(1:2))
        else
          call fatal_error("<pixels> must be length 2 in slice plot " &
               // trim(to_str(pl % id)))
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (get_arraysize_integer(node_plot, "pixels") == 3) then
          call get_node_array(node_plot, "pixels", pl % pixels(1:3))
        else
          call fatal_error("<pixels> must be length 3 in voxel plot " &
               // trim(to_str(pl % id)))
        end if
      end if

      ! Copy plot background color
      if (check_for_node(node_plot, "background")) then
        if (pl % type == PLOT_TYPE_VOXEL) then
          if (master) call warning("Background color ignored in voxel plot " &
               // trim(to_str(pl % id)))
        end if
        if (get_arraysize_integer(node_plot, "background") == 3) then
          call get_node_array(node_plot, "background", pl % not_found % rgb)
        else
          call fatal_error("Bad background RGB in plot " &
               // trim(to_str(pl % id)))
        end if
      else
        pl % not_found % rgb = (/ 255, 255, 255 /)
      end if

      ! Copy plot basis
      if (pl % type == PLOT_TYPE_SLICE) then
        temp_str = 'xy'
        if (check_for_node(node_plot, "basis")) &
             call get_node_value(node_plot, "basis", temp_str)
        temp_str = to_lower(temp_str)
        select case (trim(temp_str))
        case ("xy")
          pl % basis = PLOT_BASIS_XY
        case ("xz")
          pl % basis = PLOT_BASIS_XZ
        case ("yz")
          pl % basis = PLOT_BASIS_YZ
        case default
          call fatal_error("Unsupported plot basis '" // trim(temp_str) &
               // "' in plot " // trim(to_str(pl % id)))
        end select
      end if

      ! Copy plotting origin
      if (get_arraysize_double(node_plot, "origin") == 3) then
        call get_node_array(node_plot, "origin", pl % origin)
      else
        call fatal_error("Origin must be length 3 in plot " &
             // trim(to_str(pl % id)))
      end if

      ! Copy plotting width
      if (pl % type == PLOT_TYPE_SLICE) then
        if (get_arraysize_double(node_plot, "width") == 2) then
          call get_node_array(node_plot, "width", pl % width(1:2))
        else
          call fatal_error("<width> must be length 2 in slice plot " &
               // trim(to_str(pl % id)))
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (get_arraysize_double(node_plot, "width") == 3) then
          call get_node_array(node_plot, "width", pl % width(1:3))
        else
          call fatal_error("<width> must be length 3 in voxel plot " &
               // trim(to_str(pl % id)))
        end if
      end if

      ! Copy plot cell universe level
      if (check_for_node(node_plot, "level")) then
        call get_node_value(node_plot, "level", pl % level)

        if (pl % level < 0) then
          call fatal_error("Bad universe level in plot " &
               // trim(to_str(pl % id)))
        end if
      else
        pl % level = PLOT_LEVEL_LOWEST
      end if

      ! Copy plot color type and initialize all colors randomly
      temp_str = "cell"
      if (check_for_node(node_plot, "color")) &
           call get_node_value(node_plot, "color", temp_str)
      temp_str = to_lower(temp_str)
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
        call fatal_error("Unsupported plot color type '" // trim(temp_str) &
             // "' in plot " // trim(to_str(pl % id)))
      end select

      ! Get the number of <col_spec> nodes and get a list of them
      call get_node_list(node_plot, "col_spec", node_col_list)
      n_cols = get_list_size(node_col_list)

      ! Copy user specified colors
      if (n_cols /= 0) then

        if (pl % type == PLOT_TYPE_VOXEL) then
          if (master) call warning("Color specifications ignored in voxel &
               &plot " // trim(to_str(pl % id)))
        end if

        do j = 1, n_cols

          ! Get pointer to color spec XML node
          call get_list_item(node_col_list, j, node_col)

          ! Check and make sure 3 values are specified for RGB
          if (get_arraysize_double(node_col, "rgb") /= 3) then
            call fatal_error("Bad RGB in plot " &
                 // trim(to_str(pl % id)))
          end if

          ! Ensure that there is an id for this color specification
          if (check_for_node(node_col, "id")) then
            call get_node_value(node_col, "id", col_id)
          else
            call fatal_error("Must specify id for color specification in &
                 &plot " // trim(to_str(pl % id)))
          end if

          ! Add RGB
          if (pl % color_by == PLOT_COLOR_CELLS) then

            if (cell_dict % has_key(col_id)) then
              col_id = cell_dict % get_key(col_id)
              call get_node_array(node_col, "rgb", pl % colors(col_id) % rgb)
            else
              call fatal_error("Could not find cell " // trim(to_str(col_id)) &
                   // " specified in plot " // trim(to_str(pl % id)))
            end if

          else if (pl % color_by == PLOT_COLOR_MATS) then

            if (material_dict % has_key(col_id)) then
              col_id = material_dict % get_key(col_id)
              call get_node_array(node_col, "rgb", pl % colors(col_id) % rgb)
            else
              call fatal_error("Could not find material " &
                   // trim(to_str(col_id)) // " specified in plot " &
                   // trim(to_str(pl % id)))
            end if

          end if
        end do
      end if

      ! Deal with meshlines
      call get_node_list(node_plot, "meshlines", node_meshline_list)
      n_meshlines = get_list_size(node_meshline_list)
      if (n_meshlines /= 0) then

        if (pl % type == PLOT_TYPE_VOXEL) then
          call warning("Meshlines ignored in voxel plot " &
               // trim(to_str(pl % id)))
        end if

        select case(n_meshlines)
          case (0)
            ! Skip if no meshlines are specified
          case (1)

            ! Get pointer to meshlines
            call get_list_item(node_meshline_list, 1, node_meshlines)

            ! Check mesh type
            if (check_for_node(node_meshlines, "meshtype")) then
              call get_node_value(node_meshlines, "meshtype", meshtype)
            else
              call fatal_error("Must specify a meshtype for meshlines &
                   &specification in plot " // trim(to_str(pl % id)))
            end if

            ! Ensure that there is a linewidth for this meshlines specification
            if (check_for_node(node_meshlines, "linewidth")) then
              call get_node_value(node_meshlines, "linewidth", &
                   pl % meshlines_width)
            else
              call fatal_error("Must specify a linewidth for meshlines &
                   &specification in plot " // trim(to_str(pl % id)))
            end if

            ! Check for color
            if (check_for_node(node_meshlines, "color")) then

              ! Check and make sure 3 values are specified for RGB
              if (get_arraysize_double(node_meshlines, "color") /= 3) then
                call fatal_error("Bad RGB for meshlines color in plot " &
                     // trim(to_str(pl % id)))
              end if

              call get_node_array(node_meshlines, "color", &
                   pl % meshlines_color % rgb)
            else

              pl % meshlines_color % rgb = (/ 0, 0, 0 /)

            end if

            ! Set mesh based on type
            select case (trim(meshtype))
            case ('ufs')

              if (.not. associated(ufs_mesh)) then
                call fatal_error("No UFS mesh for meshlines on plot " &
                     // trim(to_str(pl % id)))
              end if

              pl % meshlines_mesh => ufs_mesh

            case ('cmfd')

              if (.not. cmfd_run) then
                call fatal_error("Need CMFD run to plot CMFD mesh for &
                     &meshlines on plot " // trim(to_str(pl % id)))
              end if

              select type(filt => cmfd_tallies(1) % &
                   filters(cmfd_tallies(1) % find_filter(FILTER_MESH)) % obj)
              type is (MeshFilter)
                i_mesh = filt % mesh
              end select
              pl % meshlines_mesh => meshes(i_mesh)

            case ('entropy')

              if (.not. associated(entropy_mesh)) then
                call fatal_error("No entropy mesh for meshlines on plot " &
                     // trim(to_str(pl % id)))
              end if

              if (.not. allocated(entropy_mesh % dimension)) then
                call fatal_error("No dimension specified on entropy mesh &
                     &for meshlines on plot " // trim(to_str(pl % id)))
              end if

              pl % meshlines_mesh => entropy_mesh

            case ('tally')

              ! Ensure that there is a mesh id if the type is tally
              if (check_for_node(node_meshlines, "id")) then
                call get_node_value(node_meshlines, "id", meshid)
              else
                call fatal_error("Must specify a mesh id for meshlines tally &
                     &mesh specification in plot " // trim(to_str(pl % id)))
              end if

              ! Check if the specified tally mesh exists
              if (mesh_dict % has_key(meshid)) then
                pl % meshlines_mesh => meshes(mesh_dict % get_key(meshid))
                if (meshes(meshid) % type /= LATTICE_RECT) then
                  call fatal_error("Non-rectangular mesh specified in &
                       &meshlines for plot " // trim(to_str(pl % id)))
                end if
              else
                call fatal_error("Could not find mesh " &
                     // trim(to_str(meshid)) // " specified in meshlines for &
                     &plot " // trim(to_str(pl % id)))
              end if

            case default
              call fatal_error("Invalid type for meshlines on plot " &
                    // trim(to_str(pl % id)) // ": " // trim(meshtype))
            end select

          case default
            call fatal_error("Mutliple meshlines specified in plot " &
                 // trim(to_str(pl % id)))
        end select

      end if

      ! Deal with masks
      call get_node_list(node_plot, "mask", node_mask_list)
      n_masks = get_list_size(node_mask_list)
      if (n_masks /= 0) then

        if (pl % type == PLOT_TYPE_VOXEL) then
          if (master) call warning("Mask ignored in voxel plot " &
               // trim(to_str(pl % id)))
        end if

        select case(n_masks)
          case default
            call fatal_error("Mutliple masks specified in plot " &
                 // trim(to_str(pl % id)))
          case (1)

            ! Get pointer to mask
            call get_list_item(node_mask_list, 1, node_mask)

            ! Determine how many components there are and allocate
            n_comp = 0
            n_comp = get_arraysize_integer(node_mask, "components")
            if (n_comp == 0) then
              call fatal_error("Missing <components> in mask of plot " &
                   // trim(to_str(pl % id)))
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
                  call fatal_error("Could not find cell " &
                       // trim(to_str(col_id)) // " specified in the mask in &
                       &plot " // trim(to_str(pl % id)))
                end if

              else if (pl % color_by == PLOT_COLOR_MATS) then

                if (material_dict % has_key(col_id)) then
                  iarray(j) = material_dict % get_key(col_id)
                else
                  call fatal_error("Could not find material " &
                       // trim(to_str(col_id)) // " specified in the mask in &
                       &plot " // trim(to_str(pl % id)))
                end if

              end if
            end do

            ! Alter colors based on mask information
            do j=1,size(pl % colors)
              if (.not. any(j .eq. iarray)) then
                if (check_for_node(node_mask, "background")) then
                  call get_node_array(node_mask, "background", pl % colors(j) % rgb)
                else
                  call fatal_error("Missing <background> in mask of plot " &
                       // trim(to_str(pl % id)))
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
! READ_*_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the CE and MG cross sections that may be used.
!===============================================================================

  subroutine read_ce_cross_sections_xml(libraries)
    type(Library), allocatable, intent(out) :: libraries(:)

    integer :: i           ! loop index
    integer :: n
    integer :: n_libraries
    logical :: file_exists ! does cross_sections.xml exist?
    character(MAX_WORD_LEN) :: directory ! directory with cross sections
    character(MAX_WORD_LEN) :: words(MAX_WORDS)
    character(10000) :: temp_str
    type(Node), pointer :: doc
    type(Node), pointer :: node_library
    type(NodeList), pointer :: node_library_list

    ! Check if cross_sections.xml exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find cross_sections.xml file
      call fatal_error("Cross sections XML file '" &
           // trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Reading cross sections XML file...", 5)

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

    ! Get node list of all <library>
    call get_node_list(doc, "library", node_library_list)
    n_libraries = get_list_size(node_library_list)

    ! Allocate xs_listings array
    if (n_libraries == 0) then
      call fatal_error("No cross section libraries present in cross_sections.xml &
           &file!")
    else
      allocate(libraries(n_libraries))
    end if

    do i = 1, n_libraries
      ! Get pointer to ace table XML node
      call get_list_item(node_library_list, i, node_library)

      ! Get list of materials
      if (check_for_node(node_library, "materials")) then
        call get_node_value(node_library, "materials", temp_str)
        call split_string(temp_str, words, n)
        allocate(libraries(i) % materials(n))
        libraries(i) % materials(:) = words(1:n)
      end if

      ! Get type of library
      if (check_for_node(node_library, "type")) then
        call get_node_value(node_library, "type", temp_str)
        select case(to_lower(temp_str))
        case ('neutron')
          libraries(i) % type = LIBRARY_NEUTRON
        case ('thermal')
          libraries(i) % type = LIBRARY_THERMAL
        end select
      else
        call fatal_error("Missing library type")
      end if

      ! determine path of cross section table
      if (check_for_node(node_library, "path")) then
        call get_node_value(node_library, "path", temp_str)
      else
        call fatal_error("Missing library path")
      end if

      if (starts_with(temp_str, '/')) then
        libraries(i) % path = trim(temp_str)
      else
        if (ends_with(directory,'/')) then
          libraries(i) % path = trim(directory) // trim(temp_str)
        else
          libraries(i) % path = trim(directory) // '/' // trim(temp_str)
        end if
      end if

      inquire(FILE=libraries(i) % path, EXIST=file_exists)
      if (.not. file_exists) then
        call warning("Cross section library " // trim(libraries(i) % path) // &
             " does not exist.")
      end if
    end do

    ! Close cross sections XML file
    call close_xmldoc(doc)

  end subroutine read_ce_cross_sections_xml

  subroutine read_mg_cross_sections_xml(libraries)
    type(Library), allocatable, intent(out) :: libraries(:)

    integer :: i           ! loop index
    integer :: n_libraries
    logical :: file_exists ! does cross_sections.xml exist?
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_xsdata => null()
    type(NodeList), pointer :: node_xsdata_list => null()
    real(8), allocatable :: rev_energy_bins(:)

    ! Check if mgxs.xml exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find mgxs.xml file
      call fatal_error("Cross sections XML file '" &
           // trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Reading cross sections XML file...", 5)

    ! Parse mgxs.xml file
    call open_xmldoc(doc, path_cross_sections)

    if (check_for_node(doc, "groups")) then
      ! Get neutron group count
      call get_node_value(doc, "groups", energy_groups)
    else
      call fatal_error("groups element must exist!")
    end if

    allocate(rev_energy_bins(energy_groups + 1))
    allocate(energy_bins(energy_groups + 1))
    if (check_for_node(doc, "group_structure")) then
      ! Get neutron group structure
      call get_node_array(doc, "group_structure", energy_bins)
    else
      call fatal_error("group_structures element must exist!")
    end if

    ! First reverse the order of energy_groups
    energy_bins = energy_bins(energy_groups + 1:1:-1)

    allocate(energy_bin_avg(energy_groups))
    do i = 1, energy_groups
      energy_bin_avg(i) = HALF * (energy_bins(i) + energy_bins(i + 1))
    end do

    allocate(inverse_velocities(energy_groups))
    if (check_for_node(doc, "inverse_velocities")) then
      ! Get inverse velocities
      call get_node_array(doc, "inverse_velocities", inverse_velocities)
    else
      ! If not given, estimate them by using average energy in group which is
      ! assumed to be the midpoint
      do i = 1, energy_groups
        inverse_velocities(i) = ONE / &
             (sqrt(TWO * energy_bin_avg(i) / (MASS_NEUTRON_MEV)) * &
              C_LIGHT * 100.0_8)
      end do
    end if

    ! Get node list of all <xsdata>
    call get_node_list(doc, "xsdata", node_xsdata_list)
    n_libraries = get_list_size(node_xsdata_list)

    ! Allocate xs_listings array
    if (n_libraries == 0) then
      call fatal_error("At least one <xsdata> element must be present in &
                       &mgxs.xml file!")
    else
      allocate(libraries(n_libraries))
    end if

    do i = 1, n_libraries
      ! Get pointer to xsdata table XML node
      call get_list_item(node_xsdata_list, i, node_xsdata)

      ! Get name of material
      allocate(libraries(i) % materials(1))
      call get_node_value(node_xsdata, "name", libraries(i) % materials(1))
    end do

    ! Close cross sections XML file
    call close_xmldoc(doc)

  end subroutine read_mg_cross_sections_xml

!===============================================================================
! EXPAND_NATURAL_ELEMENT converts natural elements specified using an <element>
! tag within a material into individual isotopes based on IUPAC Isotopic
! Compositions of the Elements 2009 (doi:10.1351/PAC-REP-10-06-02). In some
! cases, modifications have been made to work with ENDF/B-VII.1 where
! evaluations of particular isotopes don't exist.
!===============================================================================

  subroutine expand_natural_element(name, xs, density, names, densities)
    character(*),   intent(in)   :: name
    character(*),   intent(in)   :: xs
    real(8),        intent(in)   :: density
    type(VectorChar), intent(inout) :: names
    type(VectorReal), intent(inout) :: densities

    character(2) :: element_name

    element_name = name(1:2)

    select case (to_lower(element_name))
    case ('h')
      call names % push_back('H1.' // xs)
      call densities % push_back(density * 0.999885_8)
      call names % push_back('H2.' // xs)
      call densities % push_back(density * 0.000115_8)
    case ('he')
      call names % push_back('He3.' // xs)
      call densities % push_back(density * 0.00000134_8)
      call names % push_back('He4.' // xs)
      call densities % push_back(density * 0.99999866_8)

    case ('li')
      call names % push_back('Li6.' // xs)
      call densities % push_back(density * 0.0759_8)
      call names % push_back('Li7.' // xs)
      call densities % push_back(density * 0.9241_8)

    case ('be')
      call names % push_back('Be9.' // xs)
      call densities % push_back(density)

    case ('b')
      call names % push_back('B10.' // xs)
      call densities % push_back(density * 0.199_8)
      call names % push_back('B11.' // xs)
      call densities % push_back(density * 0.801_8)

    case ('c')
      ! No evaluations split up Carbon into isotopes yet
      call names % push_back('C0.' // xs)
      call densities % push_back(density)

    case ('n')
      call names % push_back('N14.' // xs)
      call densities % push_back(density * 0.99636_8)
      call names % push_back('N15.' // xs)
      call densities % push_back(density * 0.00364_8)

    case ('o')
      if (default_expand == JEFF_32) then
        call names % push_back('O16.' // xs)
        call densities % push_back(density * 0.99757_8)
        call names % push_back('O17.' // xs)
        call densities % push_back(density * 0.00038_8)
        call names % push_back('O18.' // xs)
        call densities % push_back(density * 0.00205_8)
      elseif (default_expand >= JENDL_32 .and. default_expand <= JENDL_40) then
        call names % push_back('O16.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('O16.' // xs)
        call densities % push_back(density * 0.99962_8)
        call names % push_back('O17.' // xs)
        call densities % push_back(density * 0.00038_8)
      end if

    case ('f')
      call names % push_back('F19.' // xs)
      call densities % push_back(density)

    case ('ne')
      call names % push_back('Ne20.' // xs)
      call densities % push_back(density * 0.9048_8)
      call names % push_back('Ne21.' // xs)
      call densities % push_back(density * 0.0027_8)
      call names % push_back('Ne22.' // xs)
      call densities % push_back(density * 0.0925_8)

    case ('na')
      call names % push_back('Na23.' // xs)
      call densities % push_back(density)

    case ('mg')
      call names % push_back('Mg24.' // xs)
      call densities % push_back(density * 0.7899_8)
      call names % push_back('Mg25.' // xs)
      call densities % push_back(density * 0.1000_8)
      call names % push_back('Mg26.' // xs)
      call densities % push_back(density * 0.1101_8)

    case ('al')
      call names % push_back('Al27.' // xs)
      call densities % push_back(density)

    case ('si')
      call names % push_back('Si28.' // xs)
      call densities % push_back(density * 0.92223_8)
      call names % push_back('Si29.' // xs)
      call densities % push_back(density * 0.04685_8)
      call names % push_back('Si30.' // xs)
      call densities % push_back(density * 0.03092_8)

    case ('p')
      call names % push_back('P31.' // xs)
      call densities % push_back(density)

    case ('s')
      call names % push_back('S32.' // xs)
      call densities % push_back(density * 0.9499_8)
      call names % push_back('S33.' // xs)
      call densities % push_back(density * 0.0075_8)
      call names % push_back('S34.' // xs)
      call densities % push_back(density * 0.0425_8)
      call names % push_back('S36.' // xs)
      call densities % push_back(density * 0.0001_8)

    case ('cl')
      call names % push_back('Cl35.' // xs)
      call densities % push_back(density * 0.7576_8)
      call names % push_back('Cl37.' // xs)
      call densities % push_back(density * 0.2424_8)

    case ('ar')
      call names % push_back('Ar36.' // xs)
      call densities % push_back(density * 0.003336_8)
      call names % push_back('Ar38.' // xs)
      call densities % push_back(density * 0.000629_8)
      call names % push_back('Ar40.' // xs)
      call densities % push_back(density * 0.996035_8)

    case ('k')
      call names % push_back('K39.' // xs)
      call densities % push_back(density * 0.932581_8)
      call names % push_back('K40.' // xs)
      call densities % push_back(density * 0.000117_8)
      call names % push_back('K41.' // xs)
      call densities % push_back(density * 0.067302_8)

    case ('ca')
      call names % push_back('Ca40.' // xs)
      call densities % push_back(density * 0.96941_8)
      call names % push_back('Ca42.' // xs)
      call densities % push_back(density * 0.00647_8)
      call names % push_back('Ca43.' // xs)
      call densities % push_back(density * 0.00135_8)
      call names % push_back('Ca44.' // xs)
      call densities % push_back(density * 0.02086_8)
      call names % push_back('Ca46.' // xs)
      call densities % push_back(density * 0.00004_8)
      call names % push_back('Ca48.' // xs)
      call densities % push_back(density * 0.00187_8)

    case ('sc')
      call names % push_back('Sc45.' // xs)
      call densities % push_back(density)

    case ('ti')
      call names % push_back('Ti46.' // xs)
      call densities % push_back(density * 0.0825_8)
      call names % push_back('Ti47.' // xs)
      call densities % push_back(density * 0.0744_8)
      call names % push_back('Ti48.' // xs)
      call densities % push_back(density * 0.7372_8)
      call names % push_back('Ti49.' // xs)
      call densities % push_back(density * 0.0541_8)
      call names % push_back('Ti50.' // xs)
      call densities % push_back(density * 0.0518_8)

    case ('v')
      if (default_expand == ENDF_BVII0 .or. default_expand == JEFF_311 &
           .or. default_expand == JEFF_32 .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_33)) then
        call names % push_back('V0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('V50.' // xs)
        call densities % push_back(density * 0.0025_8)
        call names % push_back('V51.' // xs)
        call densities % push_back(density * 0.9975_8)
      end if

    case ('cr')
      call names % push_back('Cr50.' // xs)
      call densities % push_back(density * 0.04345_8)
      call names % push_back('Cr52.' // xs)
      call densities % push_back(density * 0.83789_8)
      call names % push_back('Cr53.' // xs)
      call densities % push_back(density * 0.09501_8)
      call names % push_back('Cr54.' // xs)
      call densities % push_back(density * 0.02365_8)

    case ('mn')
      call names % push_back('Mn55.' // xs)
      call densities % push_back(density)

    case ('fe')
      call names % push_back('Fe54.' // xs)
      call densities % push_back(density * 0.05845_8)
      call names % push_back('Fe56.' // xs)
      call densities % push_back(density * 0.91754_8)
      call names % push_back('Fe57.' // xs)
      call densities % push_back(density * 0.02119_8)
      call names % push_back('Fe58.' // xs)
      call densities % push_back(density * 0.00282_8)

    case ('co')
      call names % push_back('Co59.' // xs)
      call densities % push_back(density)

    case ('ni')
      call names % push_back('Ni58.' // xs)
      call densities % push_back(density * 0.68077_8)
      call names % push_back('Ni60.' // xs)
      call densities % push_back(density * 0.26223_8)
      call names % push_back('Ni61.' // xs)
      call densities % push_back(density * 0.011399_8)
      call names % push_back('Ni62.' // xs)
      call densities % push_back(density * 0.036346_8)
      call names % push_back('Ni64.' // xs)
      call densities % push_back(density * 0.009255_8)

    case ('cu')
      call names % push_back('Cu63.' // xs)
      call densities % push_back(density * 0.6915_8)
      call names % push_back('Cu65.' // xs)
      call densities % push_back(density * 0.3085_8)

    case ('zn')
      if (default_expand == ENDF_BVII0 .or. default_expand == &
           JEFF_311 .or. default_expand == JEFF_312) then
        call names % push_back('Zn0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Zn64.' // xs)
        call densities % push_back(density * 0.4917_8)
        call names % push_back('Zn66.' // xs)
        call densities % push_back(density * 0.2773_8)
        call names % push_back('Zn67.' // xs)
        call densities % push_back(density * 0.0404_8)
        call names % push_back('Zn68.' // xs)
        call densities % push_back(density * 0.1845_8)
        call names % push_back('Zn70.' // xs)
        call densities % push_back(density * 0.0061_8)
      end if

    case ('ga')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call names % push_back('Ga0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Ha69.' // xs)
        call densities % push_back(density * 0.60108_8)
        call names % push_back('Ga71.' // xs)
        call densities % push_back(density * 0.39892_8)
      end if

    case ('ge')
      call names % push_back('Ge70.' // xs)
      call densities % push_back(density * 0.2057_8)
      call names % push_back('Ge72.' // xs)
      call densities % push_back(density * 0.2745_8)
      call names % push_back('Ge73.' // xs)
      call densities % push_back(density * 0.0775_8)
      call names % push_back('Ge74.' // xs)
      call densities % push_back(density * 0.3650_8)
      call names % push_back('Ge76.' // xs)
      call densities % push_back(density * 0.0773_8)

    case ('as')
      call names % push_back('As75.' // xs)
      call densities % push_back(density)

    case ('se')
      call names % push_back('Se74.' // xs)
      call densities % push_back(density * 0.0089_8)
      call names % push_back('Se76.' // xs)
      call densities % push_back(density * 0.0937_8)
      call names % push_back('Se77.' // xs)
      call densities % push_back(density * 0.0763_8)
      call names % push_back('Se78.' // xs)
      call densities % push_back(density * 0.2377_8)
      call names % push_back('Se80.' // xs)
      call densities % push_back(density * 0.4961_8)
      call names % push_back('Se82.' // xs)
      call densities % push_back(density * 0.0873_8)

    case ('br')
      call names % push_back('Br79.' // xs)
      call densities % push_back(density * 0.5069_8)
      call names % push_back('Br81.' // xs)
      call densities % push_back(density * 0.4931_8)

    case ('kr')
      call names % push_back('Kr78.' // xs)
      call densities % push_back(density * 0.00355_8)
      call names % push_back('Kr80.' // xs)
      call densities % push_back(density * 0.02286_8)
      call names % push_back('Kr82.' // xs)
      call densities % push_back(density * 0.11593_8)
      call names % push_back('Kr83.' // xs)
      call densities % push_back(density * 0.11500_8)
      call names % push_back('Kr84.' // xs)
      call densities % push_back(density * 0.56987_8)
      call names % push_back('Kr86.' // xs)
      call densities % push_back(density * 0.17279_8)

    case ('rb')
      call names % push_back('Rb85.' // xs)
      call densities % push_back(density * 0.7217_8)
      call names % push_back('Rb87.' // xs)
      call densities % push_back(density * 0.2783_8)

    case ('sr')
      call names % push_back('Sr84.' // xs)
      call densities % push_back(density * 0.0056_8)
      call names % push_back('Sr86.' // xs)
      call densities % push_back(density * 0.0986_8)
      call names % push_back('Sr87.' // xs)
      call densities % push_back(density * 0.0700_8)
      call names % push_back('Sr88.' // xs)
      call densities % push_back(density * 0.8258_8)

    case ('y')
      call names % push_back('Y89.' // xs)
      call densities % push_back(density)

    case ('zr')
      call names % push_back('Zr90.' // xs)
      call densities % push_back(density * 0.5145_8)
      call names % push_back('Zr91.' // xs)
      call densities % push_back(density * 0.1122_8)
      call names % push_back('Zr92.' // xs)
      call densities % push_back(density * 0.1715_8)
      call names % push_back('Zr94.' // xs)
      call densities % push_back(density * 0.1738_8)
      call names % push_back('Zr96.' // xs)
      call densities % push_back(density * 0.0280_8)

    case ('nb')
      call names % push_back('Nb93.' // xs)
      call densities % push_back(density)

    case ('mo')
      call names % push_back('Mo92.' // xs)
      call densities % push_back(density * 0.1453_8)
      call names % push_back('Mo94.' // xs)
      call densities % push_back(density * 0.0915_8)
      call names % push_back('Mo95.' // xs)
      call densities % push_back(density * 0.1584_8)
      call names % push_back('Mo96.' // xs)
      call densities % push_back(density * 0.1667_8)
      call names % push_back('Mo97.' // xs)
      call densities % push_back(density * 0.0960_8)
      call names % push_back('Mo98.' // xs)
      call densities % push_back(density * 0.2439_8)
      call names % push_back('Mo100.' // xs)
      call densities % push_back(density * 0.0982_8)

    case ('ru')
      call names % push_back('Ru96.' // xs)
      call densities % push_back(density * 0.0554_8)
      call names % push_back('Ru98.' // xs)
      call densities % push_back(density * 0.0187_8)
      call names % push_back('Ru99.' // xs)
      call densities % push_back(density * 0.1276_8)
      call names % push_back('Ru100.' // xs)
      call densities % push_back(density * 0.1260_8)
      call names % push_back('Ru101.' // xs)
      call densities % push_back(density * 0.1706_8)
      call names % push_back('Ru102.' // xs)
      call densities % push_back(density * 0.3155_8)
      call names % push_back('Ru104.' // xs)
      call densities % push_back(density * 0.1862_8)

    case ('rh')
      call names % push_back('Rh103.' // xs)
      call densities % push_back(density)

    case ('pd')
      call names % push_back('Pd102.' // xs)
      call densities % push_back(density * 0.0102_8)
      call names % push_back('Pd104.' // xs)
      call densities % push_back(density * 0.1114_8)
      call names % push_back('Pd105.' // xs)
      call densities % push_back(density * 0.2233_8)
      call names % push_back('Pd106.' // xs)
      call densities % push_back(density * 0.2733_8)
      call names % push_back('Pd108.' // xs)
      call densities % push_back(density * 0.2646_8)
      call names % push_back('Pd110.' // xs)
      call densities % push_back(density * 0.1172_8)

    case ('ag')
      call names % push_back('Ag107.' // xs)
      call densities % push_back(density * 0.51839_8)
      call names % push_back('Ag109.' // xs)
      call densities % push_back(density * 0.48161_8)

    case ('cd')
      call names % push_back('Cd106.' // xs)
      call densities % push_back(density * 0.0125_8)
      call names % push_back('Cd108.' // xs)
      call densities % push_back(density * 0.0089_8)
      call names % push_back('Cd110.' // xs)
      call densities % push_back(density * 0.1249_8)
      call names % push_back('Cd111.' // xs)
      call densities % push_back(density * 0.1280_8)
      call names % push_back('Cd112.' // xs)
      call densities % push_back(density * 0.2413_8)
      call names % push_back('Cd113.' // xs)
      call densities % push_back(density * 0.1222_8)
      call names % push_back('Cd114.' // xs)
      call densities % push_back(density * 0.2873_8)
      call names % push_back('Cd116.' // xs)
      call densities % push_back(density * 0.0749_8)

    case ('in')
      call names % push_back('In113.' // xs)
      call densities % push_back(density * 0.0429_8)
      call names % push_back('In115.' // xs)
      call densities % push_back(density * 0.9571_8)

    case ('sn')
      call names % push_back('Sn112.' // xs)
      call densities % push_back(density * 0.0097_8)
      call names % push_back('Sn114.' // xs)
      call densities % push_back(density * 0.0066_8)
      call names % push_back('Sn115.' // xs)
      call densities % push_back(density * 0.0034_8)
      call names % push_back('Sn116.' // xs)
      call densities % push_back(density * 0.1454_8)
      call names % push_back('Sn117.' // xs)
      call densities % push_back(density * 0.0768_8)
      call names % push_back('Sn118.' // xs)
      call densities % push_back(density * 0.2422_8)
      call names % push_back('Sn119.' // xs)
      call densities % push_back(density * 0.0859_8)
      call names % push_back('Sn120.' // xs)
      call densities % push_back(density * 0.3258_8)
      call names % push_back('Sn122.' // xs)
      call densities % push_back(density * 0.0463_8)
      call names % push_back('Sn124.' // xs)
      call densities % push_back(density * 0.0579_8)

    case ('sb')
      call names % push_back('Sb121.' // xs)
      call densities % push_back(density * 0.5721_8)
      call names % push_back('Sb123.' // xs)
      call densities % push_back(density * 0.4279_8)

    case ('te')
      call names % push_back('Te120.' // xs)
      call densities % push_back(density * 0.0009_8)
      call names % push_back('Te122.' // xs)
      call densities % push_back(density * 0.0255_8)
      call names % push_back('Te123.' // xs)
      call densities % push_back(density * 0.0089_8)
      call names % push_back('Te124.' // xs)
      call densities % push_back(density * 0.0474_8)
      call names % push_back('Te125.' // xs)
      call densities % push_back(density * 0.0707_8)
      call names % push_back('Te126.' // xs)
      call densities % push_back(density * 0.1884_8)
      call names % push_back('Te128.' // xs)
      call densities % push_back(density * 0.3174_8)
      call names % push_back('Te130.' // xs)
      call densities % push_back(density * 0.3408_8)

    case ('i')
      call names % push_back('I127.' // xs)
      call densities % push_back(density)

    case ('xe')
      call names % push_back('Xe124.' // xs)
      call densities % push_back(density * 0.000952_8)
      call names % push_back('Xe126.' // xs)
      call densities % push_back(density * 0.000890_8)
      call names % push_back('Xe128.' // xs)
      call densities % push_back(density * 0.019102_8)
      call names % push_back('Xe129.' // xs)
      call densities % push_back(density * 0.264006_8)
      call names % push_back('Xe130.' // xs)
      call densities % push_back(density * 0.040710_8)
      call names % push_back('Xe131.' // xs)
      call densities % push_back(density * 0.212324_8)
      call names % push_back('Xe132.' // xs)
      call densities % push_back(density * 0.269086_8)
      call names % push_back('Xe134.' // xs)
      call densities % push_back(density * 0.104357_8)
      call names % push_back('Xe136.' // xs)
      call densities % push_back(density * 0.088573_8)

    case ('cs')
      call names % push_back('Cs133.' // xs)
      call densities % push_back(density)

    case ('ba')
      call names % push_back('Ba130.' // xs)
      call densities % push_back(density * 0.00106_8)
      call names % push_back('Ba132.' // xs)
      call densities % push_back(density * 0.00101_8)
      call names % push_back('Ba134.' // xs)
      call densities % push_back(density * 0.02417_8)
      call names % push_back('Ba135.' // xs)
      call densities % push_back(density * 0.06592_8)
      call names % push_back('Ba136.' // xs)
      call densities % push_back(density * 0.07854_8)
      call names % push_back('Ba137.' // xs)
      call densities % push_back(density * 0.11232_8)
      call names % push_back('Ba138.' // xs)
      call densities % push_back(density * 0.71698_8)

    case ('la')
      call names % push_back('La138.' // xs)
      call densities % push_back(density * 0.0008881_8)
      call names % push_back('La139.' // xs)
      call densities % push_back(density * 0.9991119_8)

    case ('ce')
      call names % push_back('Ce136.' // xs)
      call densities % push_back(density * 0.00185_8)
      call names % push_back('Ce138.' // xs)
      call densities % push_back(density * 0.00251_8)
      call names % push_back('Ce140.' // xs)
      call densities % push_back(density * 0.88450_8)
      call names % push_back('Ce142.' // xs)
      call densities % push_back(density * 0.11114_8)

    case ('pr')
      call names % push_back('Pr141.' // xs)
      call densities % push_back(density)

    case ('nd')
      call names % push_back('Nd142.' // xs)
      call densities % push_back(density * 0.27152_8)
      call names % push_back('Nd143.' // xs)
      call densities % push_back(density * 0.12174_8)
      call names % push_back('Nd144.' // xs)
      call densities % push_back(density * 0.23798_8)
      call names % push_back('Nd145.' // xs)
      call densities % push_back(density * 0.08293_8)
      call names % push_back('Nd146.' // xs)
      call densities % push_back(density * 0.17189_8)
      call names % push_back('Nd148.' // xs)
      call densities % push_back(density * 0.05756_8)
      call names % push_back('Nd150.' // xs)
      call densities % push_back(density * 0.05638_8)

    case ('sm')
      call names % push_back('Sm144.' // xs)
      call densities % push_back(density * 0.0307_8)
      call names % push_back('Sm147.' // xs)
      call densities % push_back(density * 0.1499_8)
      call names % push_back('Sm148.' // xs)
      call densities % push_back(density * 0.1124_8)
      call names % push_back('Sm149.' // xs)
      call densities % push_back(density * 0.1382_8)
      call names % push_back('Sm150.' // xs)
      call densities % push_back(density * 0.0738_8)
      call names % push_back('Sm152.' // xs)
      call densities % push_back(density * 0.2675_8)
      call names % push_back('Sm154.' // xs)
      call densities % push_back(density * 0.2275_8)

    case ('eu')
      call names % push_back('Eu151.' // xs)
      call densities % push_back(density * 0.4781_8)
      call names % push_back('Eu153.' // xs)
      call densities % push_back(density * 0.5219_8)

    case ('gd')
      call names % push_back('Gd152.' // xs)
      call densities % push_back(density * 0.0020_8)
      call names % push_back('Gd154.' // xs)
      call densities % push_back(density * 0.0218_8)
      call names % push_back('Gd155.' // xs)
      call densities % push_back(density * 0.1480_8)
      call names % push_back('Gd156.' // xs)
      call densities % push_back(density * 0.2047_8)
      call names % push_back('Gd157.' // xs)
      call densities % push_back(density * 0.1565_8)
      call names % push_back('Gd158.' // xs)
      call densities % push_back(density * 0.2484_8)
      call names % push_back('Gd160.' // xs)
      call densities % push_back(density * 0.2186_8)

    case ('tb')
      call names % push_back('Tb159.' // xs)
      call densities % push_back(density)

    case ('dy')
      call names % push_back('Dy156.' // xs)
      call densities % push_back(density * 0.00056_8)
      call names % push_back('Dy158.' // xs)
      call densities % push_back(density * 0.00095_8)
      call names % push_back('Dy160.' // xs)
      call densities % push_back(density * 0.02329_8)
      call names % push_back('Dy161.' // xs)
      call densities % push_back(density * 0.18889_8)
      call names % push_back('Dy162.' // xs)
      call densities % push_back(density * 0.25475_8)
      call names % push_back('Dy163.' // xs)
      call densities % push_back(density * 0.24896_8)
      call names % push_back('Dy164.' // xs)
      call densities % push_back(density * 0.28260_8)

    case ('ho')
      call names % push_back('Ho165.' // xs)
      call densities % push_back(density)

    case ('er')
      call names % push_back('Er162.' // xs)
      call densities % push_back(density * 0.00139_8)
      call names % push_back('Er164.' // xs)
      call densities % push_back(density * 0.01601_8)
      call names % push_back('Er166.' // xs)
      call densities % push_back(density * 0.33503_8)
      call names % push_back('Er167.' // xs)
      call densities % push_back(density * 0.22869_8)
      call names % push_back('Er168.' // xs)
      call densities % push_back(density * 0.26978_8)
      call names % push_back('Er170.' // xs)
      call densities % push_back(density * 0.14910_8)

    case ('tm')
      call names % push_back('Tm169.' // xs)
      call densities % push_back(density)

    case ('yb')
      call names % push_back('Yb168.' // xs)
      call densities % push_back(density * 0.00123_8)
      call names % push_back('Yb170.' // xs)
      call densities % push_back(density * 0.02982_8)
      call names % push_back('Yb171.' // xs)
      call densities % push_back(density * 0.1409_8)
      call names % push_back('Yb172.' // xs)
      call densities % push_back(density * 0.2168_8)
      call names % push_back('Yb173.' // xs)
      call densities % push_back(density * 0.16103_8)
      call names % push_back('Yb174.' // xs)
      call densities % push_back(density * 0.32026_8)
      call names % push_back('Yb176.' // xs)
      call densities % push_back(density * 0.12996_8)

    case ('lu')
      call names % push_back('Lu175.' // xs)
      call densities % push_back(density * 0.97401_8)
      call names % push_back('Lu176.' // xs)
      call densities % push_back(density * 0.02599_8)

    case ('hf')
      call names % push_back('Hf174.' // xs)
      call densities % push_back(density * 0.0016_8)
      call names % push_back('Hf176.' // xs)
      call densities % push_back(density * 0.0526_8)
      call names % push_back('Hf177.' // xs)
      call densities % push_back(density * 0.1860_8)
      call names % push_back('Hf178.' // xs)
      call densities % push_back(density * 0.2728_8)
      call names % push_back('Hf179.' // xs)
      call densities % push_back(density * 0.1362_8)
      call names % push_back('Hf180.' // xs)
      call densities % push_back(density * 0.3508_8)

    case ('ta')
      if (default_expand == ENDF_BVII0 .or. &
           (default_expand >= JEFF_311 .and. default_expand <= JEFF_312) .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_40)) then
        call names % push_back('Ta181.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Ta180.' // xs)
        call densities % push_back(density * 0.0001201_8)
        call names % push_back('Ta181.' // xs)
        call densities % push_back(density * 0.9998799_8)
      end if

    case ('w')
      if (default_expand == ENDF_BVII0 .or. default_expand == JEFF_311 &
           .or. default_expand == JEFF_312 .or. &
           (default_expand >= JENDL_32 .and. default_expand <= JENDL_33)) then
        ! Combine W-180 with W-182
        call names % push_back('W182.' // xs)
        call densities % push_back(density * 0.2662_8)
        call names % push_back('W183.' // xs)
        call densities % push_back(density * 0.1431_8)
        call names % push_back('W184.' // xs)
        call densities % push_back(density * 0.3064_8)
        call names % push_back('W186.' // xs)
        call densities % push_back(density * 0.2843_8)
      else
        call names % push_back('W180.' // xs)
        call densities % push_back(density * 0.0012_8)
        call names % push_back('W182.' // xs)
        call densities % push_back(density * 0.2650_8)
        call names % push_back('W183.' // xs)
        call densities % push_back(density * 0.1431_8)
        call names % push_back('W184.' // xs)
        call densities % push_back(density * 0.3064_8)
        call names % push_back('W186.' // xs)
        call densities % push_back(density * 0.2843_8)
      end if

    case ('re')
      call names % push_back('Re185.' // xs)
      call densities % push_back(density * 0.3740_8)
      call names % push_back('Re187.' // xs)
      call densities % push_back(density * 0.6260_8)

    case ('os')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call names % push_back('Os0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Os184.' // xs)
        call densities % push_back(density * 0.0002_8)
        call names % push_back('Os186.' // xs)
        call densities % push_back(density * 0.0159_8)
        call names % push_back('Os187.' // xs)
        call densities % push_back(density * 0.0196_8)
        call names % push_back('Os188.' // xs)
        call densities % push_back(density * 0.1324_8)
        call names % push_back('Os189.' // xs)
        call densities % push_back(density * 0.1615_8)
        call names % push_back('Os190.' // xs)
        call densities % push_back(density * 0.2626_8)
        call names % push_back('Os192.' // xs)
        call densities % push_back(density * 0.4078_8)
      end if

    case ('ir')
      call names % push_back('Ir191.' // xs)
      call densities % push_back(density * 0.373_8)
      call names % push_back('Ir193.' // xs)
      call densities % push_back(density * 0.627_8)

    case ('pt')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call names % push_back('Pt0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Pt190.' // xs)
        call densities % push_back(density * 0.00012_8)
        call names % push_back('Pt192.' // xs)
        call densities % push_back(density * 0.00782_8)
        call names % push_back('Pt194.' // xs)
        call densities % push_back(density * 0.3286_8)
        call names % push_back('Pt195.' // xs)
        call densities % push_back(density * 0.3378_8)
        call names % push_back('Pt196.' // xs)
        call densities % push_back(density * 0.2521_8)
        call names % push_back('Pt198.' // xs)
        call densities % push_back(density * 0.07356_8)
      end if

    case ('au')
      call names % push_back('Au197.' // xs)
      call densities % push_back(density)

    case ('hg')
      call names % push_back('Hg196.' // xs)
      call densities % push_back(density * 0.0015_8)
      call names % push_back('Hg198.' // xs)
      call densities % push_back(density * 0.0997_8)
      call names % push_back('Hg199.' // xs)
      call densities % push_back(density * 0.1687_8)
      call names % push_back('Hg200.' // xs)
      call densities % push_back(density * 0.2310_8)
      call names % push_back('Hg201.' // xs)
      call densities % push_back(density * 0.1318_8)
      call names % push_back('Hg202.' // xs)
      call densities % push_back(density * 0.2986_8)
      call names % push_back('Hg204.' // xs)
      call densities % push_back(density * 0.0687_8)

    case ('tl')
      if (default_expand == JEFF_311 .or. default_expand == JEFF_312) then
        call names % push_back('Tl0.' // xs)
        call densities % push_back(density)
      else
        call names % push_back('Tl203.' // xs)
        call densities % push_back(density * 0.2952_8)
        call names % push_back('Tl205.' // xs)
        call densities % push_back(density * 0.7048_8)
      end if

    case ('pb')
      call names % push_back('Pb204.' // xs)
      call densities % push_back(density * 0.014_8)
      call names % push_back('Pb206.' // xs)
      call densities % push_back(density * 0.241_8)
      call names % push_back('Pb207.' // xs)
      call densities % push_back(density * 0.221_8)
      call names % push_back('Pb208.' // xs)
      call densities % push_back(density * 0.524_8)

    case ('bi')
      call names % push_back('Bi209.' // xs)
      call densities % push_back(density)

    case ('th')
      call names % push_back('Th232.' // xs)
      call densities % push_back(density)

    case ('pa')
      call names % push_back('Pa231.' // xs)
      call densities % push_back(density)

    case ('u')
      call names % push_back('U234.' // xs)
      call densities % push_back(density * 0.000054_8)
      call names % push_back('U235.' // xs)
      call densities % push_back(density * 0.007204_8)
      call names % push_back('U238.' // xs)
      call densities % push_back(density * 0.992742_8)

    case default
      call fatal_error("Cannot expand element: " // name)

    end select

  end subroutine expand_natural_element

!===============================================================================
! GENERATE_RPN implements the shunting-yard algorithm to generate a Reverse
! Polish notation (RPN) expression for the region specification of a cell given
! the infix notation.
!===============================================================================

  subroutine generate_rpn(cell_id, tokens, output)
    integer, intent(in) :: cell_id
    type(VectorInt), intent(in) :: tokens    ! infix notation
    type(VectorInt), intent(inout) :: output ! RPN notation

    integer :: i
    integer :: token
    integer :: op
    type(VectorInt) :: stack

    do i = 1, tokens%size()
      token = tokens%data(i)

      if (token < OP_UNION) then
        ! If token is not an operator, add it to output
        call output%push_back(token)

      elseif (token < OP_RIGHT_PAREN) then
        ! Regular operators union, intersection, complement
        do while (stack%size() > 0)
          op = stack%data(stack%size())

          if (op < OP_RIGHT_PAREN .and. &
               ((token == OP_COMPLEMENT .and. token < op) .or. &
               (token /= OP_COMPLEMENT .and. token <= op))) then
            ! While there is an operator, op, on top of the stack, if the token
            ! is left-associative and its precedence is less than or equal to
            ! that of op or if the token is right-associative and its precedence
            ! is less than that of op, move op to the output queue and push the
            ! token on to the stack. Note that only complement is
            ! right-associative.
            call output%push_back(op)
            call stack%pop_back()
          else
            exit
          end if
        end do

        call stack%push_back(token)

      elseif (token == OP_LEFT_PAREN) then
        ! If the token is a left parenthesis, push it onto the stack
        call stack%push_back(token)

      else
        ! If the token is a right parenthesis, move operators from the stack to
        ! the output queue until reaching the left parenthesis.
        do
          ! If we run out of operators without finding a left parenthesis, it
          ! means there are mismatched parentheses.
          if (stack%size() == 0) then
            call fatal_error('Mimatched parentheses in region specification &
                 &for cell ' // trim(to_str(cell_id)) // '.')
          end if

          op = stack%data(stack%size())
          if (op == OP_LEFT_PAREN) exit
          call output%push_back(op)
          call stack%pop_back()
        end do

        ! Pop the left parenthesis.
        call stack%pop_back()
      end if
    end do

    ! While there are operators on the stack, move them to the output queue
    do while (stack%size() > 0)
      op = stack%data(stack%size())

      ! If the operator is a parenthesis, it is mismatched
      if (op >= OP_RIGHT_PAREN) then
        call fatal_error('Mimatched parentheses in region specification &
             &for cell ' // trim(to_str(cell_id)) // '.')
      end if

      call output%push_back(op)
      call stack%pop_back()
    end do
  end subroutine generate_rpn

!===============================================================================
! NORMALIZE_AO normalizes the atom or weight percentages for each material
!===============================================================================

  subroutine normalize_ao()
    integer        :: i               ! index in materials array
    integer        :: j               ! index over nuclides in material
    real(8)        :: sum_percent     ! summation
    real(8)        :: awr             ! atomic weight ratio
    real(8)        :: x               ! atom percent
    logical        :: percent_in_atom ! nuclides specified in atom percent?
    logical        :: density_in_atom ! density specified in atom/b-cm?

    do i = 1, size(materials)
      associate (mat => materials(i))
        percent_in_atom = (mat % atom_density(1) > ZERO)
        density_in_atom = (mat % density > ZERO)

        sum_percent = ZERO
        do j = 1, size(mat % nuclide)
          ! determine atomic weight ratio
          if (run_CE) then
            awr = nuclides(mat % nuclide(j)) % awr
          else
            awr = ONE
          end if

          ! if given weight percent, convert all values so that they are divided
          ! by awr. thus, when a sum is done over the values, it's actually
          ! sum(w/awr)
          if (.not. percent_in_atom) then
            mat % atom_density(j) = -mat % atom_density(j) / awr
          end if
        end do

        ! determine normalized atom percents. if given atom percents, this is
        ! straightforward. if given weight percents, the value is w/awr and is
        ! divided by sum(w/awr)
        sum_percent = sum(mat % atom_density)
        mat % atom_density = mat % atom_density / sum_percent

        ! Change density in g/cm^3 to atom/b-cm. Since all values are now in atom
        ! percent, the sum needs to be re-evaluated as 1/sum(x*awr)
        if (.not. density_in_atom) then
          sum_percent = ZERO
          do j = 1, mat % n_nuclides
            if (run_CE) then
              awr = nuclides(mat % nuclide(j)) % awr
            else
              awr = ONE
            end if
            x = mat % atom_density(j)
            sum_percent = sum_percent + x*awr
          end do
          sum_percent = ONE / sum_percent
          mat%density = -mat % density * N_AVOGADRO &
               / MASS_NEUTRON * sum_percent
        end if

        ! Calculate nuclide atom densities
        mat % atom_density = mat % density * mat % atom_density
      end associate
    end do

  end subroutine normalize_ao

!===============================================================================
! ASSIGN_SAB_TABLES assigns S(alpha,beta) tables to specific nuclides within
! materials so the code knows when to apply bound thermal scattering data
!===============================================================================

  subroutine assign_sab_tables()
    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: k            ! index over S(a,b) tables in material
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting

    do i = 1, size(materials)
      ! Skip materials with no S(a,b) tables
      if (.not. allocated(materials(i) % i_sab_tables)) cycle

      associate (mat => materials(i))

        ASSIGN_SAB: do k = 1, size(mat % i_sab_tables)
          ! In order to know which nuclide the S(a,b) table applies to, we need
          ! to search through the list of nuclides for one which has a matching
          ! zaid
          associate (sab => sab_tables(mat % i_sab_tables(k)))
            FIND_NUCLIDE: do j = 1, size(mat % nuclide)
              if (any(sab % zaid == nuclides(mat % nuclide(j)) % zaid)) then
                mat % i_sab_nuclides(k) = j
                exit FIND_NUCLIDE
              end if
            end do FIND_NUCLIDE
          end associate

          ! Check to make sure S(a,b) table matched a nuclide
          if (mat % i_sab_nuclides(k) == NONE) then
            call fatal_error("S(a,b) table " // trim(mat % &
                 sab_names(k)) // " did not match any nuclide on material " &
                 // trim(to_str(mat % id)))
          end if
        end do ASSIGN_SAB

        ! If there are multiple S(a,b) tables, we need to make sure that the
        ! entries in i_sab_nuclides are sorted or else they won't be applied
        ! correctly in the cross_section module. The algorithm here is a simple
        ! insertion sort -- don't need anything fancy!

        if (size(mat % i_sab_tables) > 1) then
          SORT_SAB: do k = 2, size(mat % i_sab_tables)
            ! Save value to move
            m = k
            temp_nuclide = mat % i_sab_nuclides(k)
            temp_table   = mat % i_sab_tables(k)

            MOVE_OVER: do
              ! Check if insertion value is greater than (m-1)th value
              if (temp_nuclide >= mat % i_sab_nuclides(m-1)) exit

              ! Move values over until hitting one that's not larger
              mat % i_sab_nuclides(m) = mat % i_sab_nuclides(m-1)
              mat % i_sab_tables(m)   = mat % i_sab_tables(m-1)
              m = m - 1

              ! Exit if we've reached the beginning of the list
              if (m == 1) exit
            end do MOVE_OVER

            ! Put the original value into its new position
            mat % i_sab_nuclides(m) = temp_nuclide
            mat % i_sab_tables(m)   = temp_table
          end do SORT_SAB
        end if

        ! Deallocate temporary arrays for names of nuclides and S(a,b) tables
        if (allocated(mat % names)) deallocate(mat % names)
      end associate
    end do
  end subroutine assign_sab_tables

  subroutine read_ce_cross_sections(libraries, library_dict)
    type(Library),   intent(in)      :: libraries(:)
    type(DictCharInt), intent(inout) :: library_dict

    integer :: i, j
    integer :: i_library
    integer :: i_nuclide
    integer :: i_sab
    integer :: index_nuc_zaid   ! index in nuclide ZAID
    integer :: zaid             ! ZAID of nuclide
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    logical :: mp_found     ! if windowed multipole libraries were found
    character(MAX_WORD_LEN) :: name
    type(SetChar) :: already_read

    allocate(nuclides(n_nuclides_total))
    allocate(sab_tables(n_sab_tables))
!$omp parallel
    allocate(micro_xs(n_nuclides_total))
!$omp end parallel

    index_nuc_zaid = 0

    ! Read cross sections
    do i = 1, size(materials)
      do j = 1, size(materials(i) % names)
        name = materials(i) % names(j)

        if (.not. already_read % contains(name)) then
          i_library = library_dict % get_key(to_lower(name))
          i_nuclide = nuclide_dict % get_key(to_lower(name))

          call write_message('Reading ' // trim(name) // ' from ' // &
               trim(libraries(i_library) % path), 6)

          ! Read nuclide data from HDF5
          file_id = file_open(libraries(i_library) % path, 'r')
          group_id = open_group(file_id, name)
          call nuclides(i_nuclide) % from_hdf5(group_id)
          call close_group(group_id)
          call file_close(file_id)

          ! Assign resonant scattering data
          if (treat_res_scat) call read_0K_elastic_scattering(&
               nuclides(i_nuclide), libraries, library_dict)

          ! Determine if minimum/maximum energy for this nuclide is greater/less
          ! than the previous
          energy_min_neutron = max(energy_min_neutron, nuclides(i_nuclide) % energy(1))
          energy_max_neutron = min(energy_max_neutron, nuclides(i_nuclide) % energy(&
               size(nuclides(i_nuclide) % energy)))

          ! Add name and alias to dictionary
          call already_read % add(name)

          ! Construct dictionary mapping nuclide zaids to [1,N] -- used for
          ! unresolved resonance probability tables
          zaid = nuclides(i_nuclide) % zaid
          if (.not. nuc_zaid_dict % has_key(zaid)) then
            index_nuc_zaid  = index_nuc_zaid + 1
            call nuc_zaid_dict % add_key(zaid, index_nuc_zaid)
          end if

          ! Read multipole file into the appropriate entry on the nuclides array
          if (multipole_active) call read_multipole_data(i_nuclide)
        end if

        ! Check if material is fissionable
        if (nuclides(materials(i) % nuclide(j)) % fissionable) then
          materials(i) % fissionable = .true.
        end if
      end do
    end do

    do i = 1, size(materials)
      ! Skip materials with no S(a,b) tables
      if (.not. allocated(materials(i) % sab_names)) cycle

      do j = 1, size(materials(i) % sab_names)
        ! Get name of S(a,b) table
        name = materials(i) % sab_names(j)

        if (.not. already_read % contains(name)) then
          i_library = library_dict % get_key(to_lower(name))
          i_sab  = sab_dict % get_key(to_lower(name))

          call write_message('Reading ' // trim(name) // ' from ' // &
               trim(libraries(i_library) % path), 6)

          ! Read S(a,b) data from HDF5
          file_id = file_open(libraries(i_library) % path, 'r')
          group_id = open_group(file_id, name)
          call sab_tables(i_sab) % from_hdf5(group_id)
          call close_group(group_id)
          call file_close(file_id)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do
    end do

    n_nuc_zaid_total = index_nuc_zaid

    ! Associate S(a,b) tables with specific nuclides
    call assign_sab_tables()

    ! Show which nuclide results in lowest energy for neutron transport
    do i = 1, size(nuclides)
      if (nuclides(i) % energy(nuclides(i) % n_grid) == energy_max_neutron) then
        call write_message("Maximum neutron transport energy: " // &
             trim(to_str(energy_max_neutron)) // " MeV for " // &
             trim(adjustl(nuclides(i) % name)), 6)
        exit
      end if
    end do

    ! If the user wants multipole, make sure we found a multipole library.
    if (multipole_active) then
      mp_found = .false.
      do i = 1, size(nuclides)
        if (nuclides(i) % mp_present) then
          mp_found = .true.
          exit
        end if
      end do
      if (.not. mp_found) call warning("Windowed multipole functionality is &
           &turned on, but no multipole libraries were found.  Set the &
           &<multipole_library> element in settings.xml or the &
           &OPENMC_MULTIPOLE_LIBRARY environment variable.")
    end if

  end subroutine read_ce_cross_sections

!===============================================================================
! READ_0K_ELASTIC_SCATTERING
!===============================================================================

  subroutine read_0K_elastic_scattering(nuc, libraries, library_dict)
    type(Nuclide), intent(inout)   :: nuc
    type(Library),   intent(in)      :: libraries(:)
    type(DictCharInt), intent(inout) :: library_dict

    integer :: i, j
    integer :: i_library
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    real(8) :: xs_cdf_sum
    character(MAX_WORD_LEN) :: name
    type(Nuclide) :: resonant_nuc

    do i = 1, size(nuclides_0K)
      if (nuc % name == nuclides_0K(i) % name) then
        ! Copy basic information from settings.xml
        nuc % resonant = .true.
        nuc % name_0K = trim(nuclides_0K(i) % name_0K)
        nuc % scheme = trim(nuclides_0K(i) % scheme)
        nuc % E_min = nuclides_0K(i) % E_min
        nuc % E_max = nuclides_0K(i) % E_max

        ! Get index in libraries array
        name = nuc % name_0K
        i_library = library_dict % get_key(to_lower(name))

        call write_message('Reading ' // trim(name) // ' 0K data from ' // &
             trim(libraries(i_library) % path), 6)

        ! Read nuclide data from HDF5
        file_id = file_open(libraries(i_library) % path, 'r')
        group_id = open_group(file_id, name)
        call resonant_nuc % from_hdf5(group_id)
        call close_group(group_id)
        call file_close(file_id)

        ! Copy 0K energy grid and elastic scattering cross section
        call move_alloc(TO=nuc % energy_0K, FROM=resonant_nuc % energy)
        call move_alloc(TO=nuc % elastic_0K, FROM=resonant_nuc % elastic)
        nuc % n_grid_0K = size(nuc % energy_0K)

        ! Build CDF for 0K elastic scattering
        xs_cdf_sum = ZERO
        allocate(nuc % xs_cdf(size(nuc % energy_0K)))

        do j = 1, size(nuc % energy_0K) - 1
          ! Negative cross sections result in a CDF that is not monotonically
          ! increasing. Set all negative xs values to ZERO.
          if (nuc % elastic_0K(j) < ZERO) nuc % elastic_0K(j) = ZERO

          ! build xs cdf
          xs_cdf_sum = xs_cdf_sum &
               + (sqrt(nuc % energy_0K(j)) * nuc % elastic_0K(j) &
               + sqrt(nuc % energy_0K(j+1)) * nuc % elastic_0K(j+1)) / TWO &
               * (nuc % energy_0K(j+1) - nuc % energy_0K(j))
          nuc % xs_cdf(j) = xs_cdf_sum
        end do

        exit
      end if
    end do

  end subroutine read_0K_elastic_scattering

!===============================================================================
! READ_MULTIPOLE_DATA checks for the existence of a multipole library in the
! directory and loads it using multipole_read
!===============================================================================

  subroutine read_multipole_data(i_table)

    integer, intent(in) :: i_table  ! index in nuclides/sab_tables

    integer :: i
    logical :: file_exists                 ! Does multipole library exist?
    character(7) :: readable               ! Is multipole library readable?
    character(6) :: zaid_string            ! String of the ZAID
    character(MAX_FILE_LEN+9) :: filename  ! Path to multipole xs library

    ! For the time being, and I know this is a bit hacky, we just assume
    ! that the file will be zaid.h5.
    associate (nuc => nuclides(i_table))

      write(zaid_string, '(I6.6)') nuc % zaid
      filename = trim(path_multipole) // zaid_string // ".h5"

      ! Check if Multipole library exists and is readable
      inquire(FILE=filename, EXIST=file_exists, READ=readable)
      if (.not. file_exists) then
        nuc % mp_present = .false.
        return
      elseif (readable(1:3) == 'NO') then
        call fatal_error("Multipole library '" // trim(filename) // "' is not &
             &readable! Change file permissions with chmod command.")
      end if

      ! Display message
      call write_message("Loading Multipole XS table: " // filename, 6)

      allocate(nuc % multipole)

      ! Call the read routine
      call multipole_read(filename, nuc % multipole, i_table)
      nuc % mp_present = .true.

      ! Recreate nu-fission cross section
      if (nuc % fissionable) then
        do i = 1, size(nuc % energy)
          nuc % nu_fission(i) = nuc % nu(nuc % energy(i), EMISSION_TOTAL) * &
               nuc % fission(i)
        end do
      else
        nuc % nu_fission(:) = ZERO
      end if

    end associate

  end subroutine read_multipole_data

end module input_xml
