module input_xml

  use algorithm,        only: find
  use cmfd_input,       only: configure_cmfd
  use constants
  use dict_header,      only: DictIntInt, DictCharInt, ElemKeyValueCI
  use distribution_multivariate
  use distribution_univariate
  use endf,             only: reaction_name
  use energy_grid,      only: grid_method, n_log_bins
  use error,            only: fatal_error, warning
  use geometry_header,  only: Cell, Lattice, RectLattice, HexLattice, &
                              get_temperatures, root_universe
  use global
  use hdf5_interface
  use list_header,      only: ListChar, ListInt, ListReal
  use mesh_header,      only: RegularMesh
  use message_passing
  use mgxs_data,        only: create_macro_xs, read_mgxs
  use multipole,        only: multipole_read
  use output,           only: write_message, title, header
  use plot_header
  use random_lcg,       only: prn, seed
  use surface_header
  use set_header,       only: SetChar
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with, tokenize, split_string, &
                              zero_padded
  use tally_header,     only: TallyObject
  use tally_filter_header, only: TallyFilterContainer
  use tally_filter
  use tally_initialize, only: add_tallies
  use xml_interface

#ifdef PURXS
  use purxs_api, only:&
       URR_EVENT,&
       URR_HISTORY,&
       URR_BATCH,&
       URR_SIMULATION,&
       URR_E_NEUTRON,&
       URR_E_RESONANCE,&
       URR_SLBW,&
       URR_MLBW,&
       URR_REICH_MOORE,&
       URR_MNBW,&
       URR_ENDF6,&
       URR_USER,&
       URR_ENDFFILE,&
       URR_FALSE,&
       URR_INF,&
       URR_LINEAR,&
       URR_LOGARITHMIC,&
       URR_HISTOGRAM,&
       URR_LINEAR_LINEAR,&
       URR_LINEAR_LOG,&
       URR_LOG_LINEAR,&
       URR_LOG_LOG,&
       URR_SQRT_LINEAR,&
       URR_SQRT_LOG,&
       URR_STATISTICAL,&
       URR_LOW_NEIGHBOR,&
       URR_ON_THE_FLY,&
       URR_PROB_BANDS,&
       URR_POINTWISE,&
       URR_MIT_W,&
       URR_QUICK_W,&
       URR_isotopes,&
       URR_use_urr,&
       URR_num_isotopes,&
       URR_endf_filenames,&
       URR_formalism,&
       URR_T_grid_prob_tables,&
       URR_E_grid_prob_tables,&
       URR_num_l_waves,&
       URR_faddeeva_method,&
       URR_write_avg_xs,&
       URR_write_prob_tables,&
       URR_background_xs_treatment,&
       URR_competitive_structure,&
       URR_E_grid_scheme_prob_tables,&
       URR_i_realization_user,&
       URR_temperature_interp_scheme,&
       URR_pregenerated_point_xs,&
       URR_pregenerated_prob_tables,&
       URR_num_energies_prob_tables,&
       URR_num_temperatures_prob_tables,&
       URR_num_bands_prob_tables,&
       URR_num_histories_prob_tables,&
       URR_rel_err_tolerance_avg_xs,&
       URR_realization_frequency,&
       URR_rel_err_tolerance_pointwise,&
       URR_xs_representation,&
       URR_min_num_batches_prob_tables,&
       URR_max_num_batches_prob_tables,&
       URR_min_delta_E_pointwise,&
       URR_num_urr_realizations,&
       URR_path_avg_xs,&
       URR_path_point_xs,&
       URR_path_prob_tables,&
       URR_path_endf_files,&
       URR_parameter_energy_dependence,&
       URR_xs_source_pointwise,&
       URR_temperature_pointwise,&
       URR_RECONSTRUCTION,&
       URR_HDF5
#endif

  implicit none
  save

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, tallies, and unresolved resonance region (URR) treatment.
!===============================================================================

  subroutine read_input_xml()

    call read_urr_xml()
    if (run_mode == MODE_PURXS) return
    call read_settings_xml()
    call read_geometry_xml()
    call read_materials()
    call read_tallies_xml()
    if (cmfd_run) call configure_cmfd()

    if (.not. run_CE) then
      ! Create material macroscopic data for MGXS
      call time_read_xs % start()
      call read_mgxs()
      call create_macro_xs()
      call time_read_xs % stop()
    end if

    ! Normalize atom/weight percents
    if (run_mode /= MODE_PLOTTING) call normalize_ao()

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
    real(8), allocatable :: temp_real(:)
    integer :: n_tracks
    logical :: file_exists
    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: filename
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_mode
    type(XMLNode) :: node_source
    type(XMLNode) :: node_space
    type(XMLNode) :: node_angle
    type(XMLNode) :: node_dist
    type(XMLNode) :: node_cutoff
    type(XMLNode) :: node_entropy
    type(XMLNode) :: node_ufs
    type(XMLNode) :: node_sp
    type(XMLNode) :: node_output
    type(XMLNode) :: node_res_scat
    type(XMLNode) :: node_trigger
    type(XMLNode) :: node_vol
    type(XMLNode) :: node_tab_leg
    type(XMLNode), allocatable :: node_source_list(:)
    type(XMLNode), allocatable :: node_vol_list(:)

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
    end if

    ! Parse settings.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Verbosity
    if (check_for_node(root, "verbosity")) then
      call get_node_value(root, "verbosity", verbosity)
    end if

    ! To this point, we haven't displayed any output since we didn't know what
    ! the verbosity is. Now that we checked for it, show the title if necessary
    if (master) then
      if (verbosity >= 2) call title()
    end if
    call write_message("Reading settings XML file...", 5)

    ! Find if a multi-group or continuous-energy simulation is desired
    if (check_for_node(root, "energy_mode")) then
      call get_node_value(root, "energy_mode", temp_str)
      temp_str = trim(to_lower(temp_str))
      if (temp_str == "mg" .or. temp_str == "multi-group") then
        run_CE = .false.
      else if (temp_str == "ce" .or. temp_str == "continuous-energy") then
        run_CE = .true.
      end if
    end if

    ! Look for deprecated cross_sections.xml file in settings.xml
    if (check_for_node(root, "cross_sections")) then
      call warning("Setting cross_sections in settings.xml has been deprecated.&
           & The cross_sections are now set in materials.xml and the &
           &cross_sections input to materials.xml and the OPENMC_CROSS_SECTIONS&
           & environment variable will take precendent over setting &
           &cross_sections in settings.xml.")
      call get_node_value(root, "cross_sections", path_cross_sections)
    end if

    ! Look for deprecated windowed_multipole file in settings.xml
    if (run_mode /= MODE_PLOTTING) then
      if (check_for_node(root, "multipole_library")) then
        call warning("Setting multipole_library in settings.xml has been &
             &deprecated. The multipole_library is now set in materials.xml and&
             & the multipole_library input to materials.xml and the &
             &OPENMC_MULTIPOLE_LIBRARY environment variable will take &
             &precendent over setting multipole_library in settings.xml.")
        call get_node_value(root, "multipole_library", path_multipole)
      end if
      if (.not. ends_with(path_multipole, "/")) &
           path_multipole = trim(path_multipole) // "/"
    end if

    if (.not. run_CE) then
      ! Scattering Treatments
      if (check_for_node(root, "max_order")) then
        call get_node_value(root, "max_order", max_order)
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

    ! Check for a trigger node and get trigger information
    if (check_for_node(root, "trigger")) then
      node_trigger = root % child("trigger")

      ! Check if trigger(s) are to be turned on
      call get_node_value(node_trigger, "active", trigger_on)

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

    ! Check run mode if it hasn't been set from the command line
    if (run_mode == NONE) then
      if (check_for_node(root, "run_mode")) then
        call get_node_value(root, "run_mode", temp_str)
        select case (to_lower(temp_str))
        case ("eigenvalue")
          run_mode = MODE_EIGENVALUE
        case ("fixed source")
          run_mode = MODE_FIXEDSOURCE
        case ("plot")
          run_mode = MODE_PLOTTING
        case ("particle restart")
          run_mode = MODE_PARTICLE
        case ("volume")
          run_mode = MODE_VOLUME
        case default
          call fatal_error("Unrecognized run mode: " // &
               trim(temp_str) // ".")
        end select

        ! Assume XML specifics <particles>, <batches>, etc. directly
        node_mode = root
      else
        call warning("<run_mode> should be specified.")

        ! Make sure that either eigenvalue or fixed source was specified
        node_mode = root % child("eigenvalue")
        if (node_mode % associated()) then
          if (run_mode == NONE) run_mode = MODE_EIGENVALUE
        else
          node_mode = root % child("fixed_source")
          if (node_mode % associated()) then
            if (run_mode == NONE) run_mode = MODE_FIXEDSOURCE
          else
            call fatal_error("<eigenvalue> or <fixed_source> not specified.")
          end if
        end if
      end if
    end if

    if (run_mode == MODE_EIGENVALUE .or. run_mode == MODE_FIXEDSOURCE) then
      ! Read run parameters
      call get_run_parameters(node_mode)

      ! Check number of active batches, inactive batches, and particles
      if (n_active <= 0) then
        call fatal_error("Number of active batches must be greater than zero.")
      elseif (n_inactive < 0) then
        call fatal_error("Number of inactive batches must be non-negative.")
      elseif (n_particles <= 0) then
        call fatal_error("Number of particles must be greater than zero.")
      end if
    end if

    ! Copy random number seed if specified
    if (check_for_node(root, "seed")) call get_node_value(root, "seed", seed)

    ! Number of bins for logarithmic grid
    if (check_for_node(root, "log_grid_bins")) then
      call get_node_value(root, "log_grid_bins", n_log_bins)
      if (n_log_bins < 1) then
        call fatal_error("Number of bins for logarithmic grid must be &
             &greater than zero.")
      end if
    else
      n_log_bins = 8000
    end if

    ! Number of OpenMP threads
    if (check_for_node(root, "threads")) then
#ifdef _OPENMP
      if (n_threads == NONE) then
        call get_node_value(root, "threads", n_threads)
        if (n_threads < 1) then
          call fatal_error("Invalid number of threads: " // to_str(n_threads))
        end if
        call omp_set_num_threads(n_threads)
      end if
#else
      call warning("Ignoring number of threads.")
#endif
    end if

    ! ==========================================================================
    ! EXTERNAL SOURCE

    ! Get point to list of <source> elements and make sure there is at least one
    call get_node_list(root, "source", node_source_list)
    n = size(node_source_list)

    if (n == 0) then
      ! Default source is isotropic point source at origin with Watt spectrum
      allocate(external_source(1))
      external_source % strength = ONE

      allocate(SpatialPoint :: external_source(1) % space)
      select type (space => external_source(1) % space)
      type is (SpatialPoint)
        space % xyz(:) = [ZERO, ZERO, ZERO]
      end select

      allocate(Isotropic :: external_source(1) % angle)
      external_source(1) % angle % reference_uvw(:) = [ZERO, ZERO, ONE]

      allocate(Watt :: external_source(1) % energy)
      select type(energy => external_source(1) % energy)
      type is (Watt)
        energy % a = 0.988e6_8
        energy % b = 2.249e-6_8
      end select
    else
      ! Allocate array for sources
      allocate(external_source(n))
    end if

    ! Read each source
    do i = 1, n
      ! Get pointer to source
      node_source = node_source_list(i)

      ! Check if we want to write out source
      if (check_for_node(node_source, "write_initial")) then
        call get_node_value(node_source, "write_initial", write_initial_source)
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
          call fatal_error("Source file '" // trim(path_source) &
               // "' does not exist!")
        end if

      else

        ! Spatial distribution for external source
        if (check_for_node(node_source, "space")) then

          ! Get pointer to spatial distribution
          node_space = node_source % child("space")

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
              node_dist = node_space % child("x")
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
              node_dist = node_space % child("y")
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
              node_dist = node_space % child("z")
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
            if (node_word_count(node_space, "parameters") /= 6) then
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
            if (node_word_count(node_space, "parameters") /= 3) then
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
          ! If no spatial distribution specified, make it a point source
          allocate(SpatialPoint :: external_source(i) % space)
          select type (space => external_source(i) % space)
          type is (SpatialPoint)
            space % xyz(:) = [ZERO, ZERO, ZERO]
          end select
        end if

        ! Determine external source angular distribution
        if (check_for_node(node_source, "angle")) then

          ! Get pointer to angular distribution
          node_angle = node_source % child("angle")

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
            n = node_word_count(node_angle, "reference_uvw")
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
              node_dist = node_angle % child("mu")
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
              node_dist = node_angle % child("phi")
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
          node_dist = node_source % child("energy")
          call distribution_from_xml(external_source(i)%energy, node_dist)
        else
          ! Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
          allocate(Watt :: external_source(i)%energy)
          select type(energy => external_source(i)%energy)
          type is (Watt)
            energy%a = 0.988e6_8
            energy%b = 2.249e-6_8
          end select
        end if
      end if
    end do

    ! Survival biasing
    if (check_for_node(root, "survival_biasing")) then
      call get_node_value(root, "survival_biasing", survival_biasing)
    end if

    ! Probability tables
    if (check_for_node(root, "ptables")) then
      call get_node_value(root, "ptables", urr_ptables_on)
    end if

    ! Cutoffs
    if (check_for_node(root, "cutoff")) then
      node_cutoff = root % child("cutoff")
      if (check_for_node(node_cutoff, "weight")) then
        call get_node_value(node_cutoff, "weight", weight_cutoff)
      end if
      if (check_for_node(node_cutoff, "weight_avg")) then
        call get_node_value(node_cutoff, "weight_avg", weight_survive)
      end if
      if (check_for_node(node_cutoff, "energy")) then
        call get_node_value(node_cutoff, "energy", energy_cutoff)
      end if
    end if

    ! Particle trace
    if (check_for_node(root, "trace")) then
      call get_node_array(root, "trace", temp_int_array3)
      trace_batch    = temp_int_array3(1)
      trace_gen      = temp_int_array3(2)
      trace_particle = int(temp_int_array3(3), 8)
    end if

    ! Particle tracks
    if (check_for_node(root, "track")) then
      ! Make sure that there are three values per particle
      n_tracks = node_word_count(root, "track")
      if (mod(n_tracks, 3) /= 0) then
        call fatal_error("Number of integers specified in 'track' is not divisible &
             &by 3.  Please provide 3 integers per particle to be tracked.")
      end if

      ! Allocate space and get list of tracks
      allocate(temp_int_array(n_tracks))
      call get_node_array(root, "track", temp_int_array)

      ! Reshape into track_identifiers
      allocate(track_identifiers(3, n_tracks/3))
      track_identifiers = reshape(temp_int_array, [3, n_tracks/3])
    end if

    ! Shannon Entropy mesh
    if (check_for_node(root, "entropy")) then

      ! Get pointer to entropy node
      node_entropy = root % child("entropy")

      ! Check to make sure enough values were supplied
      if (node_word_count(node_entropy, "lower_left") /= 3) then
        call fatal_error("Need to specify (x,y,z) coordinates of lower-left &
             &corner of Shannon entropy mesh.")
      elseif (node_word_count(node_entropy, "upper_right") /= 3) then
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
      if (.not. all(entropy_mesh % upper_right > entropy_mesh % lower_left)) then
        call fatal_error("Upper-right coordinate must be greater than lower-left &
             &coordinate for Shannon entropy mesh.")
      end if

      ! Check if dimensions were specified -- if not, they will be calculated
      ! automatically upon first entry into shannon_entropy
      if (check_for_node(node_entropy, "dimension")) then

        ! If so, make sure proper number of values were given
        if (node_word_count(node_entropy, "dimension") /= 3) then
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
    if (check_for_node(root, "uniform_fs")) then

      ! Get pointer to ufs node
      node_ufs = root % child("uniform_fs")

      ! Check to make sure enough values were supplied
      if (node_word_count(node_ufs, "lower_left") /= 3) then
        call fatal_error("Need to specify (x,y,z) coordinates of lower-left &
             &corner of UFS mesh.")
      elseif (node_word_count(node_ufs, "upper_right") /= 3) then
        call fatal_error("Need to specify (x,y,z) coordinates of upper-right &
             &corner of UFS mesh.")
      elseif (node_word_count(node_ufs, "dimension") /= 3) then
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
        call fatal_error("Upper-right coordinate must be greater than lower-left &
             &coordinate for UFS mesh.")
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
    if (check_for_node(root, "state_point")) then

      ! Get pointer to state_point node
      node_sp = root % child("state_point")

      ! Determine number of batches at which to store state points
      if (check_for_node(node_sp, "batches")) then
        n_state_points = node_word_count(node_sp, "batches")
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
    if (check_for_node(root, "source_point")) then

      ! Get pointer to source_point node
      node_sp = root % child("source_point")

      ! Determine number of batches at which to store source points
      if (check_for_node(node_sp, "batches")) then
        n_source_points = node_word_count(node_sp, "batches")
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
      else
        ! If neither were specified, write source points with state points
        n_source_points = n_state_points
        do i = 1, n_state_points
          call sourcepoint_batch % add(statepoint_batch % get_item(i))
        end do
      end if

      ! Check if the user has specified to write binary source file
      if (check_for_node(node_sp, "separate")) then
        call get_node_value(node_sp, "separate", source_separate)
      end if
      if (check_for_node(node_sp, "write")) then
        call get_node_value(node_sp, "write", source_write)
      end if
      if (check_for_node(node_sp, "overwrite_latest")) then
        call get_node_value(node_sp, "overwrite_latest", source_latest)
        source_separate = source_latest
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
    if (check_for_node(root, "no_reduce")) then
      call get_node_value(root, "no_reduce", reduce_tallies)
    end if

    ! Check if the user has specified to use confidence intervals for
    ! uncertainties rather than standard deviations
    if (check_for_node(root, "confidence_intervals")) then
      call get_node_value(root, "confidence_intervals", confidence_intervals)
    end if

    ! Check for output options
    if (check_for_node(root, "output")) then

      ! Get pointer to output node
      node_output = root % child("output")

      ! Check for summary option
      if (check_for_node(node_output, "summary")) then
        call get_node_value(node_output, "summary", output_summary)
      end if

      ! Check for ASCII tallies output option
      if (check_for_node(node_output, "tallies")) then
        call get_node_value(node_output, "tallies", output_tallies)
      end if

      ! Set output directory if a path has been specified
      if (check_for_node(node_output, "path")) then
        call get_node_value(node_output, "path", path_output)
        if (.not. ends_with(path_output, "/")) &
             path_output = trim(path_output) // "/"
      end if
    end if

    ! Check for cmfd run
    if (check_for_node(root, "run_cmfd")) then
      call get_node_value(root, "run_cmfd", cmfd_run)
    end if

    ! Resonance scattering parameters
    if (check_for_node(root, "resonance_scattering")) then
      node_res_scat = root % child("resonance_scattering")

      ! See if resonance scattering is enabled
      if (check_for_node(node_res_scat, "enable")) then
        call get_node_value(node_res_scat, "enable", res_scat_on)
      else
        res_scat_on = .true.
      end if

      ! Determine what method is used
      if (check_for_node(node_res_scat, "method")) then
        call get_node_value(node_res_scat, "method", temp_str)
        select case(to_lower(temp_str))
        case ('ares')
          res_scat_method = RES_SCAT_ARES
        case ('dbrc')
          res_scat_method = RES_SCAT_DBRC
        case ('wcm')
          res_scat_method = RES_SCAT_WCM
        case default
          call fatal_error("Unrecognized resonance elastic scattering method: " &
               // trim(temp_str) // ".")
        end select
      end if

      ! Minimum energy for resonance scattering
      if (check_for_node(node_res_scat, "energy_min")) then
        call get_node_value(node_res_scat, "energy_min", res_scat_energy_min)
      end if
      if (res_scat_energy_min < ZERO) then
        call fatal_error("Lower resonance scattering energy bound is negative")
      end if

      ! Maximum energy for resonance scattering
      if (check_for_node(node_res_scat, "energy_max")) then
        call get_node_value(node_res_scat, "energy_max", res_scat_energy_max)
      end if
      if (res_scat_energy_max < res_scat_energy_min) then
        call fatal_error("Upper resonance scattering energy bound is below the &
             &lower resonance scattering energy bound.")
      end if

      ! Get nuclides that resonance scattering should be applied to
      if (check_for_node(node_res_scat, "nuclides")) then
        n = node_word_count(node_res_scat, "nuclides")
        allocate(res_scat_nuclides(n))
        if (n > 0) then
          call get_node_array(node_res_scat, "nuclides", res_scat_nuclides)
        end if
      end if
    end if

    call get_node_list(root, "volume_calc", node_vol_list)
    n = size(node_vol_list)
    allocate(volume_calcs(n))
    do i = 1, n
      node_vol = node_vol_list(i)
      call volume_calcs(i) % from_xml(node_vol)
    end do

    ! Get temperature settings
    if (check_for_node(root, "temperature_default")) then
      call get_node_value(root, "temperature_default", temperature_default)
    end if
    if (check_for_node(root, "temperature_method")) then
      call get_node_value(root, "temperature_method", temp_str)
      select case (to_lower(temp_str))
      case ('nearest')
        temperature_method = TEMPERATURE_NEAREST
      case ('interpolation')
        temperature_method = TEMPERATURE_INTERPOLATION
      case default
        call fatal_error("Unknown temperature method: " // trim(temp_str))
      end select
    end if
    if (check_for_node(root, "temperature_tolerance")) then
      call get_node_value(root, "temperature_tolerance", temperature_tolerance)
    end if
    if (check_for_node(root, "temperature_multipole")) then
      call get_node_value(root, "temperature_multipole", temperature_multipole)
    end if

    ! Check for tabular_legendre options
    if (check_for_node(root, "tabular_legendre")) then

      ! Get pointer to tabular_legendre node
      node_tab_leg = root % child("tabular_legendre")

      ! Check for enable option
      if (check_for_node(node_tab_leg, "enable")) then
        call get_node_value(node_tab_leg, "enable", legendre_to_tabular)
      end if

      ! Check for the number of points
      if (check_for_node(node_tab_leg, "num_points")) then
        call get_node_value(node_tab_leg, "num_points", &
             legendre_to_tabular_points)
        if (legendre_to_tabular_points <= 1 .and. (.not. run_CE)) then
          call fatal_error("The 'num_points' subelement/attribute of the &
               &'tabular_legendre' element must contain a value greater than 1")
        end if
      end if
    end if

    ! Check whether create fission sites
    if (run_mode == MODE_FIXEDSOURCE) then
      if (check_for_node(root, "create_fission_neutrons")) then
        call get_node_value(root, "create_fission_neutrons", &
             create_fission_neutrons)
      end if
    end if

    ! Close settings XML file
    call doc % clear()

  end subroutine read_settings_xml

!===============================================================================
! GET_RUN_PARAMETERS
!===============================================================================

  subroutine get_run_parameters(node_base)
    type(XMLNode), intent(in) :: node_base

    character(MAX_LINE_LEN) :: temp_str
    type(XMLNode) :: node_keff_trigger

    ! Check number of particles
    if (.not. check_for_node(node_base, "particles")) then
      call fatal_error("Need to specify number of particles.")
    end if

    ! Get number of particles if it wasn't specified as a command-line argument
    if (n_particles == 0) then
      call get_node_value(node_base, "particles", n_particles)
    end if

    ! Get number of basic batches
    call get_node_value(node_base, "batches", n_batches)
    if (.not. trigger_on) then
      n_max_batches = n_batches
    end if
    n_inactive = 0
    gen_per_batch = 1

    ! Get number of inactive batches
    if (run_mode == MODE_EIGENVALUE) then
      call get_node_value(node_base, "inactive", n_inactive)
      if (check_for_node(node_base, "generations_per_batch")) then
        call get_node_value(node_base, "generations_per_batch", gen_per_batch)
      end if

      ! Allocate array for batch keff and entropy
      allocate(k_generation(n_max_batches*gen_per_batch))
      allocate(entropy(n_max_batches*gen_per_batch))
      entropy = ZERO

      ! Get the trigger information for keff
      if (check_for_node(node_base, "keff_trigger")) then
        node_keff_trigger = node_base % child("keff_trigger")

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

    ! Determine number of active batches
    n_active = n_batches - n_inactive

  end subroutine get_run_parameters

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i, j, k, m, i_x, i_a, input_index
    integer :: n, n_mats, n_x, n_y, n_z, n_rings, n_rlats, n_hlats
    integer :: univ_id
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
    character(:), allocatable :: region_spec
    type(Cell),     pointer :: c
    class(Surface), pointer :: s
    class(Lattice), pointer :: lat
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_cell
    type(XMLNode) :: node_surf
    type(XMLNode) :: node_lat
    type(XMLNode), allocatable :: node_cell_list(:)
    type(XMLNode), allocatable :: node_surf_list(:)
    type(XMLNode), allocatable :: node_rlat_list(:)
    type(XMLNode), allocatable :: node_hlat_list(:)
    type(VectorInt) :: tokens
    type(VectorInt) :: rpn
    type(VectorInt) :: fill_univ_ids ! List of fill universe IDs
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains

    call write_message("Reading geometry XML file...")

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Geometry XML file '" // trim(filename) // "' does not exist!")
    end if

    ! Parse geometry.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Get pointer to list of XML <cell>
    call get_node_list(root, "cell", node_cell_list)

    ! Get number of <cell> tags
    n_cells = size(node_cell_list)

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
      node_cell = node_cell_list(i)

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
        c % universe = 0
      end if
      if (check_for_node(node_cell, "fill")) then
        call get_node_value(node_cell, "fill", c % fill)
        if (find(fill_univ_ids, c % fill) == -1) &
             call fill_univ_ids % push_back(c % fill)
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
        n_mats = node_word_count(node_cell, "material")

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
      if (check_for_node(node_cell, "surfaces")) then
        call warning("The use of 'surfaces' is deprecated and will be &
             &disallowed in a future release.  Use 'region' instead. The &
             &openmc-update-inputs utility can be used to automatically &
             &update geometry.xml files.")
        region_spec = node_value_string(node_cell, "surfaces")
        call get_node_value(node_cell, "surfaces", region_spec)
      elseif (check_for_node(node_cell, "region")) then
        region_spec = node_value_string(node_cell, "region")
      else
        region_spec = ''
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
               c % id)) // " because it is not filled with another universe")
        end if

        ! Read number of rotation parameters
        n = node_word_count(node_cell, "rotation")
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
        n = node_word_count(node_cell, "translation")
        if (n /= 3) then
          call fatal_error("Incorrect number of translation parameters on cell " &
               // to_str(c % id))
        end if

        ! Copy translation vector
        allocate(c % translation(3))
        call get_node_array(node_cell, "translation", c % translation)
      end if

      ! Read cell temperatures.  If the temperature is not specified, set it to
      ! ERROR_REAL for now.  During initialization we'll replace ERROR_REAL with
      ! the temperature from the material data.
      if (check_for_node(node_cell, "temperature")) then
        n = node_word_count(node_cell, "temperature")
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
      univ_id = c % universe
      if (.not. cells_in_univ_dict % has_key(univ_id)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % add_key(univ_id, n_universes)
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get_key(univ_id)
      end if
      call cells_in_univ_dict % add_key(univ_id, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! This variable is used to check whether at least one boundary condition was
    ! applied to a surface
    boundary_exists = .false.

    ! get pointer to list of xml <surface>
    call get_node_list(root, "surface", node_surf_list)

    ! Get number of <surface> tags
    n_surfaces = size(node_surf_list)

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
      node_surf = node_surf_list(i)

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

      n = node_word_count(node_surf, "coeffs")
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
    if (run_mode /= MODE_PLOTTING) then
      if (.not. boundary_exists) then
        call fatal_error("No boundary conditions were applied to any surfaces!")
      end if
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

        type is (SurfacePlane)
          if (surf % i_periodic == NONE) then
            call fatal_error("No matching periodic surface specified for &
                 &periodic boundary condition on surface " // &
                 trim(to_str(surf % id)) // ".")
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
    call get_node_list(root, "lattice", node_rlat_list)
    call get_node_list(root, "hex_lattice", node_hlat_list)

    ! Allocate lattices array
    n_rlats = size(node_rlat_list)
    n_hlats = size(node_hlat_list)
    n_lattices = n_rlats + n_hlats
    allocate(lattices(n_lattices))

    RECT_LATTICES: do i = 1, n_rlats
      allocate(RectLattice::lattices(i) % obj)
      lat => lattices(i) % obj
      select type(lat)
      type is (RectLattice)

      ! Get pointer to i-th lattice
      node_lat = node_rlat_list(i)

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
      n = node_word_count(node_lat, "dimension")
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
      if (node_word_count(node_lat, "lower_left") /= n) then
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
        if (node_word_count(node_lat, "width") /= n) then
          call fatal_error("Number of entries on <pitch> must be the same as &
               &the number of entries on <dimension>.")
        end if

      else if (node_word_count(node_lat, "pitch") /= n) then
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
      n = node_word_count(node_lat, "universes")
      if (n /= n_x*n_y*n_z) then
        call fatal_error("Number of universes on <universes> does not match size of &
             &lattice " // trim(to_str(lat % id)) // ".")
      end if

      allocate(temp_int_array(n))
      call get_node_array(node_lat, "universes", temp_int_array)

      ! Read universes
      do m = 1, n_z
        do k = 0, n_y - 1
          do j = 1, n_x
            lat % universes(j, n_y - k, m) = &
                 temp_int_array(j + n_x*k + n_x*n_y*(m-1))
            if (find(fill_univ_ids, lat % universes(j, n_y - k, m)) == -1) &
                 call fill_univ_ids % push_back(lat % universes(j, n_y - k, m))
          end do
        end do
      end do
      deallocate(temp_int_array)

      ! Read outer universe for area outside lattice.
      lat % outer = NO_OUTER_UNIVERSE
      if (check_for_node(node_lat, "outer")) then
        call get_node_value(node_lat, "outer", lat % outer)
        if (find(fill_univ_ids, lat % outer) == -1) &
             call fill_univ_ids % push_back(lat % outer)
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
      node_lat = node_hlat_list(i)

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
      n = node_word_count(node_lat, "center")
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
      n = node_word_count(node_lat, "pitch")
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
      n = node_word_count(node_lat, "universes")
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
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
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
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
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
                 temp_int_array(input_index)
            if (find(fill_univ_ids, temp_int_array(input_index)) == -1) &
                 call fill_univ_ids % push_back(temp_int_array(input_index))
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
        if (find(fill_univ_ids, lat % outer) == -1) &
             call fill_univ_ids % push_back(lat % outer)
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

    ! ==========================================================================
    ! SETUP UNIVERSES

    ! Allocate universes, universe cell arrays, and assign base universe
    allocate(universes(n_universes))
    do i = 1, n_universes
      associate (u => universes(i))
        u % id = univ_ids % data(i)

        ! Allocate cell list
        n_cells_in_univ = cells_in_univ_dict % get_key(u % id)
        allocate(u % cells(n_cells_in_univ))
        u % cells(:) = 0

        ! Check whether universe is a fill universe
        if (find(fill_univ_ids, u % id) == -1) then
          if (root_universe > 0) then
            call fatal_error("Two or more universes are not used as fill &
                 &universes, so it is not possible to distinguish which one &
                 &is the root universe.")
          else
            root_universe = i
          end if
        end if
      end associate
    end do

    do i = 1, n_cells
      ! Get index in universes array
      j = universe_dict % get_key(cells(i) % universe)

      ! Set the first zero entry in the universe cells array to the index in the
      ! global cells array
      associate (u => universes(j))
        u % cells(find(u % cells, 0)) = i
      end associate
    end do

    ! Clear dictionary
    call cells_in_univ_dict%clear()

    ! Close geometry XML file
    call doc % clear()

  end subroutine read_geometry_xml

  subroutine read_materials()
    integer :: i, j
    type(DictCharInt) :: library_dict
    type(Library), allocatable :: libraries(:)
    type(VectorReal), allocatable :: nuc_temps(:) ! List of T to read for each nuclide
    type(VectorReal), allocatable :: sab_temps(:) ! List of T to read for each S(a,b)
    real(8), allocatable    :: material_temps(:)
    logical                 :: file_exists
    character(MAX_FILE_LEN) :: env_variable
    character(MAX_LINE_LEN) :: filename
    type(XMLDocument)       :: doc
    type(XMLNode)           :: root

    call write_message("Reading materials XML file...")

    ! Check is materials.xml exists
    filename = trim(path_input) // "materials.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Material XML file '" // trim(filename) // "' does not exist!")
    end if

    ! Parse materials.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Find cross_sections.xml file -- the first place to look is the
    ! materials.xml file. If no file is found there, then we check the
    ! OPENMC_CROSS_SECTIONS environment variable
    if (.not. check_for_node(root, "cross_sections")) then
      ! No cross_sections.xml file specified in settings.xml, check
      ! environment variable
      if (run_CE) then
        call get_environment_variable("OPENMC_CROSS_SECTIONS", env_variable)
        if (len_trim(env_variable) == 0) then
          call get_environment_variable("CROSS_SECTIONS", env_variable)
          ! FIXME: When deprecated option of setting the cross sections in
          ! settings.xml is removed, remove ".and. path_cross_sections == ''"
          if (len_trim(env_variable) == 0 .and. path_cross_sections == '') then
            call fatal_error("No cross_sections.xml file was specified in &
                 &materials.xml, settings.xml,  or in the OPENMC_CROSS_SECTIONS&
                 & environment variable. OpenMC needs such a file to identify &
                 &where to find ACE cross section libraries. Please consult the&
                 & user's guide at http://mit-crpg.github.io/openmc for &
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
          ! FIXME: When deprecated option of setting the mg cross sections in
          ! settings.xml is removed, remove ".and. path_cross_sections == ''"
        if (len_trim(env_variable) == 0 .and. path_cross_sections == '') then
          call fatal_error("No mgxs.h5 file was specified in &
               &materials.xml or in the OPENMC_MG_CROSS_SECTIONS environment &
               &variable. OpenMC needs such a file to identify where to &
               &find MG cross section libraries. Please consult the user's &
               &guide at http://mit-crpg.github.io/openmc for information on &
               &how to set up MG cross section libraries.")
        else if (len_trim(env_variable) /= 0) then
          path_cross_sections = trim(env_variable)
        end if
      end if
    else
      call get_node_value(root, "cross_sections", path_cross_sections)
    end if

    ! Find the windowed multipole library
    if (run_mode /= MODE_PLOTTING) then
      if (.not. check_for_node(root, "multipole_library")) then
        ! No library location specified in materials.xml, check
        ! environment variable
        call get_environment_variable("OPENMC_MULTIPOLE_LIBRARY", env_variable)
        path_multipole = trim(env_variable)
      else
        call get_node_value(root, "multipole_library", path_multipole)
      end if
      if (.not. ends_with(path_multipole, "/")) &
           path_multipole = trim(path_multipole) // "/"
    end if

    ! Close materials XML file
    call doc % clear()

    ! Now that the cross_sections.xml or mgxs.h5 has been located, read it in
    if (run_CE) then
      call read_ce_cross_sections_xml(libraries)
    else
      call read_mg_cross_sections_header(libraries)
    end if

    ! Creating dictionary that maps the name of the material to the entry
    do i = 1, size(libraries)
      do j = 1, size(libraries(i) % materials)
        call library_dict % add_key(to_lower(libraries(i) % materials(j)), i)
      end do
    end do

    ! Check that 0K nuclides are listed in the cross_sections.xml file
    if (allocated(res_scat_nuclides)) then
      do i = 1, size(res_scat_nuclides)
        if (.not. library_dict % has_key(to_lower(res_scat_nuclides(i)))) then
          call fatal_error("Could not find resonant scatterer " &
               // trim(res_scat_nuclides(i)) // " in cross_sections.xml file!")
        end if
      end do
    end if

    ! Parse data from materials.xml
    call read_materials_xml(libraries, library_dict, material_temps)

    ! Assign temperatures to cells that don't have temperatures already assigned
    call assign_temperatures(material_temps)

    ! Determine desired temperatures for each nuclide and S(a,b) table
    call get_temperatures(cells, materials, material_dict, nuclide_dict, &
                          n_nuclides_total, nuc_temps, sab_dict, &
                          n_sab_tables, sab_temps)

    ! Read continuous-energy cross sections
    if (run_CE .and. run_mode /= MODE_PLOTTING) then
      call time_read_xs % start()
      call read_ce_cross_sections(libraries, library_dict, nuc_temps, sab_temps)
      call time_read_xs % stop()
    end if

    ! Clear dictionary
    call library_dict % clear()
  end subroutine read_materials

!===============================================================================
! READ_MATERIAL_XML reads data from a materials.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_materials_xml(libraries, library_dict, material_temps)
    type(Library), intent(in) :: libraries(:)
    type(DictCharInt), intent(inout) :: library_dict
    real(8), allocatable, intent(out) :: material_temps(:)

    integer :: i              ! loop index for materials
    integer :: j              ! loop index for nuclides
    integer :: n              ! number of nuclides
    integer :: n_sab          ! number of sab tables for a material
    integer :: i_library      ! index in libraries array
    integer :: index_nuclide  ! index in nuclides
    integer :: index_sab      ! index in sab_tables
    logical :: file_exists    ! does materials.xml exist?
    character(20)           :: name         ! name of nuclide, e.g. 92235.03c
    character(MAX_WORD_LEN) :: units        ! units on density
    character(MAX_LINE_LEN) :: filename     ! absolute path to materials.xml
    character(MAX_LINE_LEN) :: temp_str     ! temporary string when reading
    real(8)                 :: val          ! value entered for density
    real(8)                 :: temp_dble    ! temporary double prec. real
    logical                 :: sum_density  ! density is sum of nuclide densities
    type(VectorChar)        :: names        ! temporary list of nuclide names
    type(VectorInt)         :: list_iso_lab ! temporary list of isotropic lab scatterers
    type(VectorReal)        :: densities    ! temporary list of nuclide densities
    type(Material), pointer :: mat => null()
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_mat
    type(XMLNode) :: node_dens
    type(XMLNode) :: node_nuc
    type(XMLNode) :: node_sab
    type(XMLNode), allocatable :: node_mat_list(:)
    type(XMLNode), allocatable :: node_nuc_list(:)
    type(XMLNode), allocatable :: node_ele_list(:)
    type(XMLNode), allocatable :: node_macro_list(:)
    type(XMLNode), allocatable :: node_sab_list(:)

    ! Check is materials.xml exists
    filename = trim(path_input) // "materials.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Material XML file '" // trim(filename) // "' does not &
           &exist!")
    end if

    ! Parse materials.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Get pointer to list of XML <material>
    call get_node_list(root, "material", node_mat_list)

    ! Allocate cells array
    n_materials = size(node_mat_list)
    allocate(materials(n_materials))
    allocate(material_temps(n_materials))

    ! Initialize count for number of nuclides/S(a,b) tables
    index_nuclide = 0
    index_sab = 0

    do i = 1, n_materials
      mat => materials(i)

      ! Get pointer to i-th material node
      node_mat = node_mat_list(i)

      ! Copy material id
      if (check_for_node(node_mat, "id")) then
        call get_node_value(node_mat, "id", mat % id)
      else
        call fatal_error("Must specify id of material in materials XML file")
      end if

      ! Check if material is depletable
      if (check_for_node(node_mat, "depletable")) then
        call get_node_value(node_mat, "depletable", temp_str)
        if (to_lower(temp_str) == "true" .or. temp_str == "1") &
             mat % depletable = .true.
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

      ! Get material default temperature
      if (check_for_node(node_mat, "temperature")) then
        call get_node_value(node_mat, "temperature", material_temps(i))
      else
        material_temps(i) = ERROR_REAL
      end if

      ! Get pointer to density element
      if (check_for_node(node_mat, "density")) then
        node_dens = node_mat % child("density")
      else
        call fatal_error("Must specify density element in material " &
             // trim(to_str(mat % id)))
      end if

      ! Copy units
      call get_node_value(node_dens, "units", units)

      ! If the units is 'sum', then the total density of the material is taken
      ! to be the sum of the atom fractions listed on the nuclides
      if (units == 'sum') then
        sum_density = .true.

      else if (units == 'macro') then
        if (check_for_node(node_dens, "value")) then
          call get_node_value(node_dens, "value", val)
        else
          val = ONE
        end if

        ! Set density
        mat % density = val

        sum_density = .false.

      else
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

      ! Issue error if elements are provided
      call get_node_list(node_mat, "element", node_ele_list)

      if (size(node_ele_list) > 0) then
        call fatal_error("Unable to add an element to material " &
             // trim(to_str(mat % id)) // " since the element option has &
             &been removed from the xml input. Elements can only be added via &
             &the Python API, which will expand elements into their natural &
             &nuclides.")
      end if

      ! =======================================================================
      ! READ AND PARSE <nuclide> TAGS

      ! Check to ensure material has at least one nuclide
      if (.not. check_for_node(node_mat, "nuclide") .and. &
           .not. check_for_node(node_mat, "macroscopic")) then
        call fatal_error("No macroscopic data or nuclides specified on &
             &material " // trim(to_str(mat % id)))
      end if

      ! Create list of macroscopic x/s based on those specified, just treat
      ! them as nuclides. This is all really a facade so the user thinks they
      ! are entering in macroscopic data but the code treats them the same
      ! as nuclides internally.
      ! Get pointer list of XML <macroscopic>
      call get_node_list(node_mat, "macroscopic", node_macro_list)
      if (run_CE .and. (size(node_macro_list) > 0)) then
        call fatal_error("Macroscopic can not be used in continuous-energy&
                         & mode!")
      else if (size(node_macro_list) > 1) then
        call fatal_error("Only one macroscopic object permitted per material, " &
             // trim(to_str(mat % id)))
      else if (size(node_macro_list) == 1) then

        node_nuc = node_macro_list(1)

        ! Check for empty name on nuclide
        if (.not. check_for_node(node_nuc, "name")) then
          call fatal_error("No name specified on macroscopic data in material " &
               // trim(to_str(mat % id)))
        end if

        ! store nuclide name
        call get_node_value(node_nuc, "name", name)
        name = trim(name)

        ! save name to list
        call names % push_back(name)

        ! Set density for macroscopic data
        if (units == 'macro') then
          call densities % push_back(ONE)
        else
          call fatal_error("Units can only be macro for macroscopic data " &
               // trim(name))
        end if
      else

        ! Get pointer list of XML <nuclide>
        call get_node_list(node_mat, "nuclide", node_nuc_list)

        ! Create list of nuclides based on those specified
        INDIVIDUAL_NUCLIDES: do j = 1, size(node_nuc_list)
          ! Combine nuclide identifier and cross section and copy into names
          node_nuc = node_nuc_list(j)

          ! Check for empty name on nuclide
          if (.not. check_for_node(node_nuc, "name")) then
            call fatal_error("No name specified on nuclide in material " &
                 // trim(to_str(mat % id)))
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

          ! store nuclide name
          call get_node_value(node_nuc, "name", name)
          name = trim(name)

          ! save name to list
          call names % push_back(name)

          ! Check if no atom/weight percents were specified or if both atom and
          ! weight percents were specified
          if (units == 'macro') then
            call densities % push_back(ONE)
          else
            if (.not. check_for_node(node_nuc, "ao") .and. &
                 .not. check_for_node(node_nuc, "wo")) then
              call fatal_error("No atom or weight percent specified for &
                   &nuclide" // trim(name))
            elseif (check_for_node(node_nuc, "ao") .and. &
                    check_for_node(node_nuc, "wo")) then
              call fatal_error("Cannot specify both atom and weight percents &
                   &for a nuclide: " // trim(name))
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

        n_sab = size(node_sab_list)
        if (n_sab > 0) then
          ! Set number of S(a,b) tables
          mat % n_sab = n_sab

          ! Allocate names and indices for nuclides and tables -- for now we
          ! allocate these as the number of S(a,b) tables listed. Since a single
          ! table might apply to multiple nuclides, they are resized later if a
          ! table is indeed applied to multiple nuclides.
          allocate(mat % sab_names(n_sab))
          allocate(mat % i_sab_tables(n_sab))

          do j = 1, n_sab
            ! Get pointer to S(a,b) table
            node_sab = node_sab_list(j)

            ! Determine name of S(a,b) table
            if (.not. check_for_node(node_sab, "name")) then
              call fatal_error("Need to specify <name> for S(a,b) table.")
            end if
            call get_node_value(node_sab, "name", name)
            name = trim(name)
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
    call doc % clear()

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
    integer :: filter_id     ! user-specified identifier for filter
    integer :: i_mesh        ! index in meshes array
    integer :: i_filt        ! index in filters array
    integer :: i_filter_mesh ! index of mesh filter
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filter      ! number of filters
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
    integer, allocatable :: temp_filter(:) ! temporary filter indices
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    type(DictCharInt) :: trigger_scores
    type(ElemKeyValueCI), pointer :: pair_list
    type(TallyObject), pointer :: t
    type(TallyFilterContainer), pointer :: f
    type(RegularMesh), pointer :: m
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_mesh
    type(XMLNode) :: node_tal
    type(XMLNode) :: node_filt
    type(XMLNode) :: node_trigger
    type(XMLNode) :: node_deriv
    type(XMLNode), allocatable :: node_mesh_list(:)
    type(XMLNode), allocatable :: node_tal_list(:)
    type(XMLNode), allocatable :: node_filt_list(:)
    type(XMLNode), allocatable :: node_trigger_list(:)
    type(XMLNode), allocatable :: node_deriv_list(:)
    type(ElemKeyValueCI), pointer :: scores
    type(ElemKeyValueCI), pointer :: next

    ! Check if tallies.xml exists
    filename = trim(path_input) // "tallies.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! We need to allocate tally_derivs to avoid segfaults.  Also needs to be
      ! done in parallel because tally derivs are threadprivate.
!$omp parallel
      allocate(tally_derivs(0))
!$omp end parallel

      ! Since a tallies.xml file is optional, no error is issued here
      return
    end if

    call write_message("Reading tallies XML file...")

    ! Parse tallies.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! ==========================================================================
    ! DETERMINE SIZE OF ARRAYS AND ALLOCATE

    ! Get pointer list to XML <mesh>
    call get_node_list(root, "mesh", node_mesh_list)

    ! Get pointer list to XML <filter>
    call get_node_list(root, "filter", node_filt_list)

    ! Get pointer list to XML <tally>
    call get_node_list(root, "tally", node_tal_list)

    ! Check for user meshes
    n_user_meshes = size(node_mesh_list)
    if (cmfd_run) then
      n_meshes = n_user_meshes + n_cmfd_meshes
    else
      n_meshes = n_user_meshes
    end if

    ! Allocate mesh array
    if (n_meshes > 0) allocate(meshes(n_meshes))

    ! Check for user filters
    n_user_filters = size(node_filt_list)

    ! Allocate filters array
    if (n_user_filters > 0) call add_filters(n_user_filters)

    ! Check for user tallies
    n_user_tallies = size(node_tal_list)
    if (n_user_tallies == 0) then
      call warning("No tallies present in tallies.xml file!")
    end if

    ! Allocate tally array
    if (n_user_tallies > 0) then
      call add_tallies("user", n_user_tallies)
    end if

    ! Check for <assume_separate> setting
    if (check_for_node(root, "assume_separate")) then
      call get_node_value(root, "assume_separate", assume_separate)
    end if

    ! ==========================================================================
    ! READ MESH DATA

    do i = 1, n_user_meshes
      m => meshes(i)

      ! Get pointer to mesh node
      node_mesh = node_mesh_list(i)

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
      n = node_word_count(node_mesh, "dimension")
      if (n /= 1 .and. n /= 2 .and. n /= 3) then
        call fatal_error("Mesh must be one, two, or three dimensions.")
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
        call fatal_error("All entries on the <dimension> element for a tally mesh &
             &must be positive.")
      end if

      ! Read dimensions in each direction
      m % dimension = iarray3(1:n)

      ! Read mesh lower-left corner location
      if (m % n_dimension /= node_word_count(node_mesh, "lower_left")) then
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
        if (node_word_count(node_mesh, "width") /= &
             node_word_count(node_mesh, "lower_left")) then
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
        if (node_word_count(node_mesh, "upper_right") /= &
             node_word_count(node_mesh, "lower_left")) then
          call fatal_error("Number of entries on <upper_right> must be the &
               &same as the number of entries on <lower_left>.")
        end if

        ! Check that upper-right is above lower-left
        call get_node_array(node_mesh, "upper_right", rarray3(1:n))
        if (any(rarray3(1:n) < m % lower_left)) then
          call fatal_error("The <upper_right> coordinates must be greater than the &
               &<lower_left> coordinates on a tally mesh.")
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
    ! READ DATA FOR DERIVATIVES

    ! Get pointer list to XML <derivative> nodes and allocate global array.
    ! The array is threadprivate so it must be allocated in parallel.
    call get_node_list(root, "derivative", node_deriv_list)
!$omp parallel
    allocate(tally_derivs(size(node_deriv_list)))
!$omp end parallel

    ! Make sure this is not an MG run.
    if (.not. run_CE .and. size(node_deriv_list) > 0) then
      call fatal_error("Differential tallies not supported in multi-group mode")
    end if

    ! Read derivative attributes.
    do i = 1, size(node_deriv_list)
      associate(deriv => tally_derivs(i))
        ! Get pointer to derivative node.
        node_deriv = node_deriv_list(i)

        ! Copy the derivative id.
        if (check_for_node(node_deriv, "id")) then
          call get_node_value(node_deriv, "id", deriv % id)
        else
          call fatal_error("Must specify an ID for <derivative> elements in the&
               & tally XML file")
        end if

        ! Make sure the id is > 0.
        if (deriv % id <= 0) then
          call fatal_error("<derivative> IDs must be an integer greater than &
               &zero")
        end if

        ! Make sure this id has not already been used.
        do j = 1, i-1
          if (tally_derivs(j) % id == deriv % id) then
            call fatal_error("Two or more <derivative>'s use the same unique &
                 &ID: " // trim(to_str(deriv % id)))
          end if
        end do

        ! Read the independent variable name.
        temp_str = ""
        call get_node_value(node_deriv, "variable", temp_str)
        temp_str = to_lower(temp_str)

        select case(temp_str)

        case("density")
          deriv % variable = DIFF_DENSITY
          call get_node_value(node_deriv, "material", deriv % diff_material)

        case("nuclide_density")
          deriv % variable = DIFF_NUCLIDE_DENSITY
          call get_node_value(node_deriv, "material", deriv % diff_material)

          call get_node_value(node_deriv, "nuclide", word)
          word = trim(to_lower(word))
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
                 // trim(word) // " specified in derivative " &
                 // trim(to_str(deriv % id)) // " in any material.")
          end if
          deallocate(pair_list)

          deriv % diff_nuclide = nuclide_dict % get_key(word)

        case("temperature")
          deriv % variable = DIFF_TEMPERATURE
          call get_node_value(node_deriv, "material", deriv % diff_material)
        end select
      end associate
    end do

    ! ==========================================================================
    ! READ FILTER DATA

    READ_FILTERS: do i = 1, n_user_filters
      f => filters(i)

      ! Get pointer to filter xml node
      node_filt = node_filt_list(i)

      ! Copy filter id
      if (check_for_node(node_filt, "id")) then
        call get_node_value(node_filt, "id", filter_id)
      else
        call fatal_error("Must specify id for filter in tally XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (filter_dict % has_key(filter_id)) then
        call fatal_error("Two or more filters use the same unique ID: " &
             // to_str(filter_id))
      end if

      ! Convert filter type to lower case
      temp_str = ''
      if (check_for_node(node_filt, "type")) &
           call get_node_value(node_filt, "type", temp_str)
      temp_str = to_lower(temp_str)

      ! Determine number of bins
      select case(temp_str)
      case ("energy", "energyout", "mu", "polar", "azimuthal")
        if (.not. check_for_node(node_filt, "bins")) then
          call fatal_error("Bins not set in filter " // trim(to_str(filter_id)))
        end if
        n_words = node_word_count(node_filt, "bins")
      case ("mesh", "universe", "material", "cell", "distribcell", &
            "cellborn", "surface", "delayedgroup")
        if (.not. check_for_node(node_filt, "bins")) then
          call fatal_error("Bins not set in filter " // trim(to_str(filter_id)))
        end if
        n_words = node_word_count(node_filt, "bins")
      end select

      ! Determine type of filter
      select case (temp_str)

      case ('distribcell')
        ! Allocate and declare the filter type
        allocate(DistribcellFilter :: f % obj)
        select type (filt => f % obj)
        type is (DistribcellFilter)
          if (n_words /= 1) call fatal_error("Only one cell can be &
               &specified per distribcell filter.")
          ! Store bins
          call get_node_value(node_filt, "bins", filt % cell)
        end select

      case ('cell')
        ! Allocate and declare the filter type
        allocate(CellFilter :: f % obj)
        select type (filt => f % obj)
        type is (CellFilter)
          ! Allocate and store bins
          filt % n_bins = n_words
          allocate(filt % cells(n_words))
          call get_node_array(node_filt, "bins", filt % cells)
        end select

      case ('cellborn')
        ! Allocate and declare the filter type
        allocate(CellbornFilter :: f % obj)
        select type (filt => f % obj)
        type is (CellbornFilter)
          ! Allocate and store bins
          filt % n_bins = n_words
          allocate(filt % cells(n_words))
          call get_node_array(node_filt, "bins", filt % cells)
        end select

      case ('material')
        ! Allocate and declare the filter type
        allocate(MaterialFilter :: f % obj)
        select type (filt => f % obj)
        type is (MaterialFilter)
          ! Allocate and store bins
          filt % n_bins = n_words
          allocate(filt % materials(n_words))
          call get_node_array(node_filt, "bins", filt % materials)
        end select

      case ('universe')
        ! Allocate and declare the filter type
        allocate(UniverseFilter :: f % obj)
        select type (filt => f % obj)
        type is (UniverseFilter)
          ! Allocate and store bins
          filt % n_bins = n_words
          allocate(filt % universes(n_words))
          call get_node_array(node_filt, "bins", filt % universes)
        end select

      case ('surface')
        call fatal_error("Surface filter is not yet supported!")
        ! Allocate and declare the filter type
        allocate(SurfaceFilter :: f % obj)
        select type (filt => f % obj)
        type is (SurfaceFilter)
          ! Allocate and store bins
          filt % n_bins = n_words
          allocate(filt % surfaces(n_words))
          call get_node_array(node_filt, "bins", filt % surfaces)
        end select

      case ('mesh')
        ! Allocate and declare the filter type
        allocate(MeshFilter :: f % obj)
        select type (filt => f % obj)
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
                 // " specified on filter " // trim(to_str(filter_id)))
          end if

          ! Determine number of bins
          filt % n_bins = product(m % dimension)

          ! Store the index of the mesh
          filt % mesh = i_mesh
        end select

      case ('energy')

        ! Allocate and declare the filter type
        allocate(EnergyFilter :: f % obj)
        select type (filt => f % obj)
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
            if (n_words == num_energy_groups + 1) then
              if (all(filt % bins == energy_bins(num_energy_groups + 1:1:-1))) &
                   then
                filt % matches_transport_groups = .true.
              end if
            end if
          end if
        end select

      case ('energyout')
        ! Allocate and declare the filter type
        allocate(EnergyoutFilter :: f % obj)
        select type (filt => f % obj)
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
            if (n_words == num_energy_groups + 1) then
              if (all(filt % bins == energy_bins(num_energy_groups + 1:1:-1))) &
                   then
                filt % matches_transport_groups = .true.
              end if
            end if
          end if
        end select

      case ('delayedgroup')

        ! Allocate and declare the filter type
        allocate(DelayedGroupFilter :: f % obj)
        select type (filt => f % obj)
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

      case ('mu')
        ! Allocate and declare the filter type
        allocate(MuFilter :: f % obj)
        select type (filt => f % obj)
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
                   & greater than 1 on filter " &
                   // trim(to_str(filter_id)) // ".")
            end if
          end if
        end select

      case ('polar')
        ! Allocate and declare the filter type
        allocate(PolarFilter :: f % obj)
        select type (filt => f % obj)
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
              call fatal_error("Number of bins for polar filter must be&
                   & greater than 1 on filter " &
                   // trim(to_str(filter_id)) // ".")
            end if
          end if
        end select

      case ('azimuthal')
        ! Allocate and declare the filter type
        allocate(AzimuthalFilter :: f % obj)
        select type (filt => f % obj)
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
              call fatal_error("Number of bins for azimuthal filter must be&
                   & greater than 1 on filter " &
                   // trim(to_str(filter_id)) // ".")
            end if
          end if
        end select

      case ('energyfunction')
        ! Allocate and declare the filter type.
        allocate(EnergyFunctionFilter :: f % obj)
        select type (filt => f % obj)
        type is (EnergyFunctionFilter)
          filt % n_bins = 1
          ! Make sure this is continuous-energy mode.
          if (.not. run_CE) then
            call fatal_error("EnergyFunction filters are only supported for &
                 &continuous-energy transport calculations")
          end if

          ! Allocate and store energy grid.
          if (.not. check_for_node(node_filt, "energy")) then
            call fatal_error("Energy grid not specified for EnergyFunction &
                 &filter on filter " // trim(to_str(filter_id)))
          end if
          n_words = node_word_count(node_filt, "energy")
          allocate(filt % energy(n_words))
          call get_node_array(node_filt, "energy", filt % energy)

          ! Allocate and store interpolant values.
          if (.not. check_for_node(node_filt, "y")) then
            call fatal_error("y values not specified for EnergyFunction &
                 &filter on filter " // trim(to_str(filter_id)))
          end if
          n_words = node_word_count(node_filt, "y")
          allocate(filt % y(n_words))
          call get_node_array(node_filt, "y", filt % y)
        end select

      case default
        ! Specified filter is invalid, raise error
        call fatal_error("Unknown filter type '" &
             // trim(temp_str) // "' on filter " &
             // trim(to_str(filter_id)) // ".")

      end select

      ! Set filter id
      f % obj % id = filter_id

      ! Add filter to dictionary
      call filter_dict % add_key(filter_id, i)

    end do READ_FILTERS

    ! ==========================================================================
    ! READ TALLY DATA

    READ_TALLIES: do i = 1, n_user_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Get pointer to tally xml node
      node_tal = node_tal_list(i)

      ! Set tally type to volume by default
      t % type = TALLY_VOLUME

      ! It's desirable to use a track-length esimator for tallies since
      ! generally more events will score to the tally, reducing the
      ! variance. However, for tallies that require information on
      ! post-collision parameters (e.g. tally with an energyout filter) the
      ! analog esimator must be used.

      t % estimator = ESTIMATOR_TRACKLENGTH

      ! Copy tally id
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

      ! Determine number of filters
      if (check_for_node(node_tal, "filters")) then
        n_filter = node_word_count(node_tal, "filters")
      else
        n_filter = 0
      end if

      ! Allocate and store filter user ids
      allocate(temp_filter(n_filter))
      if (n_filter > 0) then
        call get_node_array(node_tal, "filters", temp_filter)
      end if

      do j = 1, n_filter
        ! Get pointer to filter
        if (filter_dict % has_key(temp_filter(j))) then
          i_filt = filter_dict % get_key(temp_filter(j))
          f => filters(i_filt)
        else
          call fatal_error("Could not find filter " &
               // trim(to_str(temp_filter(j))) // " specified on tally " &
               // trim(to_str(t % id)))
        end if

        ! Set the filter index in the tally find_filter array
        select type (filt => f % obj)
        type is (DistribcellFilter)
          t % find_filter(FILTER_DISTRIBCELL) = j
        type is (CellFilter)
          t % find_filter(FILTER_CELL) = j
        type is (CellbornFilter)
          t % find_filter(FILTER_CELLBORN) = j
        type is (MaterialFilter)
          t % find_filter(FILTER_MATERIAL) = j
        type is (UniverseFilter)
          t % find_filter(FILTER_UNIVERSE) = j
        type is (SurfaceFilter)
          t % find_filter(FILTER_SURFACE) = j
        type is (MeshFilter)
          t % find_filter(FILTER_MESH) = j
        type is (EnergyFilter)
          t % find_filter(FILTER_ENERGYIN) = j
        type is (EnergyoutFilter)
          t % find_filter(FILTER_ENERGYOUT) = j
          ! Set to analog estimator
          t % estimator = ESTIMATOR_ANALOG
        type is (DelayedGroupFilter)
          t % find_filter(FILTER_DELAYEDGROUP) = j
        type is (MuFilter)
          t % find_filter(FILTER_MU) = j
          ! Set to analog estimator
          t % estimator = ESTIMATOR_ANALOG
        type is (PolarFilter)
          t % find_filter(FILTER_POLAR) = j
        type is (AzimuthalFilter)
          t % find_filter(FILTER_AZIMUTHAL) = j
        type is (EnergyFunctionFilter)
          t % find_filter(FILTER_ENERGYFUNCTION) = j
        end select

        ! Store the index of the filter
        temp_filter(j) = i_filt
      end do

      ! Store the filter indices
      call move_alloc(FROM=temp_filter, TO=t % filter)

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
        allocate(sarray(node_word_count(node_tal, "nuclides")))
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
          n_words = node_word_count(node_tal, "nuclides")
          allocate(t % nuclide_bins(n_words))
          do j = 1, n_words

            ! Check if total material was specified
            if (trim(sarray(j)) == 'total') then
              t % nuclide_bins(j) = -1
              cycle
            end if

            ! If a specific nuclide was specified
            word = to_lower(sarray(j))

            ! Search through nuclides
            pair_list => nuclide_dict % keys()
            do while (associated(pair_list))
              if (trim(pair_list % key) == trim(word)) then
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
      end if

      ! =======================================================================
      ! READ DATA FOR SCORES

      if (check_for_node(node_tal, "scores")) then
        n_words = node_word_count(node_tal, "scores")
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
                call warning("Invalid scattering order of " // trim(to_str(n_order)) // &
                     " requested. Setting to the maximum permissible value, " // &
                     trim(to_str(MAX_ANG_ORDER)))
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
              call fatal_error("Cannot tally flux with an outgoing energy filter.")
            end if

          case ('flux-yn')
            ! Prohibit user from tallying flux for an individual nuclide
            if (.not. (t % n_nuclide_bins == 1 .and. &
                 t % nuclide_bins(1) == -1)) then
              call fatal_error("Cannot tally flux for an individual nuclide.")
            end if

            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally flux with an outgoing energy filter.")
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
            i_filter_mesh = t % filter(t % find_filter(FILTER_MESH))

            ! Check to make sure mesh filter was specified
            if (i_filter_mesh == 0) then
              call fatal_error("Cannot tally surface current without a mesh &
                   &filter.")
            end if

            ! Get pointer to mesh
            select type(filt => filters(i_filter_mesh) % obj)
            type is (MeshFilter)
              i_mesh = filt % mesh
              m => meshes(i_mesh)
            end select

            ! Copy filter indices to temporary array
            allocate(temp_filter(size(t % filter) + 1))
            temp_filter(1:size(t % filter)) = t % filter

            ! Move allocation back -- temp_filter becomes deallocated during
            ! this call
            call move_alloc(FROM=temp_filter, TO=t % filter)
            n_filter = size(t % filter)

            ! Extend the filters array so we can add a surface filter and mesh
            ! filter
            call add_filters(2)

            ! Increment number of user filters
            n_user_filters = n_user_filters + 2

            ! Get index of the new mesh filter
            i_filt = n_user_filters - 1

            ! Duplicate the mesh filter since other tallies might use this
            ! filter and we need to change the dimension
            allocate(MeshFilter :: filters(i_filt) % obj)
            select type(filt => filters(i_filt) % obj)
            type is (MeshFilter)
              filt % id = i_filt
              filt % mesh = i_mesh

              ! We need to increase the dimension by one since we also need
              ! currents coming into and out of the boundary mesh cells.
              filt % n_bins = product(m % dimension + 1)

              ! Add filter to dictionary
              call filter_dict % add_key(filt % id, i_filt)
            end select
            t % filter(t % find_filter(FILTER_MESH)) = i_filt

            ! Get index of the new surface filter
            i_filt = n_user_filters

            ! Add surface filter
            allocate(SurfaceFilter :: filters(i_filt) % obj)
            select type (filt => filters(i_filt) % obj)
            type is (SurfaceFilter)
              filt % id = i_filt
              filt % n_bins = 4 * m % n_dimension
              allocate(filt % surfaces(4 * m % n_dimension))
              if (m % n_dimension == 1) then
                filt % surfaces = (/ OUT_LEFT, IN_LEFT, OUT_RIGHT, IN_RIGHT /)
              elseif (m % n_dimension == 2) then
                filt % surfaces = (/ OUT_LEFT, IN_LEFT, OUT_RIGHT, IN_RIGHT, &
                     OUT_BACK, IN_BACK, OUT_FRONT, IN_FRONT /)
              elseif (m % n_dimension == 3) then
                filt % surfaces = (/ OUT_LEFT, IN_LEFT, OUT_RIGHT, IN_RIGHT, &
                     OUT_BACK, IN_BACK, OUT_FRONT, IN_FRONT, OUT_BOTTOM, &
                     IN_BOTTOM, OUT_TOP, IN_TOP /)
              end if
              filt % current = .true.

              ! Add filter to dictionary
              call filter_dict % add_key(filt % id, i_filt)
            end select
            t % find_filter(FILTER_SURFACE) = n_filter
            t % filter(n_filter) = i_filt

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

      ! Check for a tally derivative.
      if (check_for_node(node_tal, "derivative")) then
        ! Temporarily store the derivative id.
        call get_node_value(node_tal, "derivative", t % deriv)

        ! Find the derivative with the given id, and store it's index.
        do j = 1, size(tally_derivs)
          if (tally_derivs(j) % id == t % deriv) then
            t % deriv = j
            ! Only analog or collision estimators are supported for differential
            ! tallies.
            if (t % estimator == ESTIMATOR_TRACKLENGTH) then
              t % estimator = ESTIMATOR_COLLISION
            end if
            ! We found the derivative we were looking for; exit the do loop.
            exit
          end if
          if (j == size(tally_derivs)) then
            call fatal_error("Could not find derivative " &
                 // trim(to_str(t % deriv)) // " specified on tally " &
                 // trim(to_str(t % id)))
          end if
        end do

        if (tally_derivs(t % deriv) % variable == DIFF_NUCLIDE_DENSITY &
             .or. tally_derivs(t % deriv) % variable == DIFF_TEMPERATURE) then
          if (any(t % nuclide_bins == -1)) then
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Error on tally " // trim(to_str(t % id)) &
                   // ": Cannot use a 'nuclide_density' or 'temperature' &
                   &derivative on a tally with an outgoing energy filter and &
                   &'total' nuclide rate. Instead, tally each nuclide in the &
                   &material individually.")
              ! Note that diff tallies with these characteristics would work
              ! correctly if no tally events occur in the perturbed material
              ! (e.g. pertrubing moderator but only tallying fuel), but this
              ! case would be hard to check for by only reading inputs.
            end if
          end if
        end if
      end if

      ! If settings.xml trigger is turned on, create tally triggers
      if (trigger_on) then

        ! Get list of trigger nodes for this tally
        call get_node_list(node_tal, "trigger", node_trigger_list)

        ! Initialize the number of triggers
        n_user_trig = size(node_trigger_list)

        ! Count the number of triggers needed for all scores including "all"
        t % n_triggers = 0
        COUNT_TRIGGERS: do user_trig_ind = 1, n_user_trig

          ! Get pointer to trigger node
          node_trigger = node_trigger_list(user_trig_ind)

          ! Get scores for this trigger
          if (check_for_node(node_trigger, "scores")) then
            n_words = node_word_count(node_trigger, "scores")
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
          node_trigger = node_trigger_list(user_trig_ind)

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
            n_words = node_word_count(node_trigger, "scores")
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
    call doc % clear()

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
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_plot
    type(XMLNode) :: node_col
    type(XMLNode) :: node_mask
    type(XMLNode) :: node_meshlines
    type(XMLNode), allocatable :: node_plot_list(:)
    type(XMLNode), allocatable :: node_col_list(:)
    type(XMLNode), allocatable :: node_mask_list(:)
    type(XMLNode), allocatable :: node_meshline_list(:)

    ! Check if plots.xml exists
    filename = trim(path_input) // "plots.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Plots XML file '" // trim(filename) &
           // "' does not exist!")
    end if

    call write_message("Reading plot XML file...")

    ! Parse plots.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Get list pointer to XML <plot>
    call get_node_list(root, "plot", node_plot_list)

    ! Allocate plots array
    n_plots = size(node_plot_list)
    allocate(plots(n_plots))

    READ_PLOTS: do i = 1, n_plots
      pl => plots(i)

      ! Get pointer to plot XML node
      node_plot = node_plot_list(i)

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
      filename = "plot_" // trim(to_str(pl % id))
      if (check_for_node(node_plot, "filename")) &
           call get_node_value(node_plot, "filename", filename)
      select case (pl % type)
      case (PLOT_TYPE_SLICE)
        pl % path_plot = trim(path_input) // trim(filename) // ".ppm"
      case (PLOT_TYPE_VOXEL)
        pl % path_plot = trim(path_input) // trim(filename) // ".h5"
      end select

      ! Copy plot pixel size
      if (pl % type == PLOT_TYPE_SLICE) then
        if (node_word_count(node_plot, "pixels") == 2) then
          call get_node_array(node_plot, "pixels", pl % pixels(1:2))
        else
          call fatal_error("<pixels> must be length 2 in slice plot " &
               // trim(to_str(pl % id)))
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (node_word_count(node_plot, "pixels") == 3) then
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
        if (node_word_count(node_plot, "background") == 3) then
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
      if (node_word_count(node_plot, "origin") == 3) then
        call get_node_array(node_plot, "origin", pl % origin)
      else
        call fatal_error("Origin must be length 3 in plot " &
             // trim(to_str(pl % id)))
      end if

      ! Copy plotting width
      if (pl % type == PLOT_TYPE_SLICE) then
        if (node_word_count(node_plot, "width") == 2) then
          call get_node_array(node_plot, "width", pl % width(1:2))
        else
          call fatal_error("<width> must be length 2 in slice plot " &
               // trim(to_str(pl % id)))
        end if
      else if (pl % type == PLOT_TYPE_VOXEL) then
        if (node_word_count(node_plot, "width") == 3) then
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
      if (check_for_node(node_plot, "color_by")) &
           call get_node_value(node_plot, "color_by", temp_str)
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

      case ("material")

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

      ! Get the number of <color> nodes and get a list of them
      call get_node_list(node_plot, "color", node_col_list)
      n_cols = size(node_col_list)

      ! Copy user specified colors
      if (n_cols /= 0) then

        if (pl % type == PLOT_TYPE_VOXEL) then
          call warning("Color specifications ignored in voxel plot " // &
                     trim(to_str(pl % id)))
        end if

        do j = 1, n_cols

          ! Get pointer to color spec XML node
          node_col = node_col_list(j)

          ! Check and make sure 3 values are specified for RGB
          if (node_word_count(node_col, "rgb") /= 3) then
            call fatal_error("Bad RGB in plot " &
                 // trim(to_str(pl % id)))
          end if

          ! Ensure that there is an id for this color specification
          if (check_for_node(node_col, "id")) then
            call get_node_value(node_col, "id", col_id)
          else
            call fatal_error("Must specify id for color specification in plot " // &
                      trim(to_str(pl % id)))
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
      n_meshlines = size(node_meshline_list)
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
            node_meshlines = node_meshline_list(1)

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
              if (node_word_count(node_meshlines, "color") /= 3) then
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

              select type(filt => filters(cmfd_tallies(1) % &
                   filter(cmfd_tallies(1) % find_filter(FILTER_MESH))) % obj)
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
      n_masks = size(node_mask_list)
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
            node_mask = node_mask_list(1)

            ! Determine how many components there are and allocate
            n_comp = 0
            n_comp = node_word_count(node_mask, "components")
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
            do j = 1, size(pl % colors)
              if (.not. any(j == iarray)) then
                if (check_for_node(node_mask, "background")) then
                  call get_node_array(node_mask, "background", pl % colors(j) % rgb)
                else
                  pl % colors(j) % rgb(:) = [255, 255, 255]
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
    call doc % clear()

  end subroutine read_plots_xml

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Read data from a urr.xml file and parse it, placing properly-formatted data in
! the right data structures
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine read_urr_xml()
#ifdef PURXS
    character(len=:), allocatable :: temp_str
    character(MAX_LINE_LEN) :: temp_str_fixed_len
    character(len=:), allocatable :: filename
    integer :: i ! re-usable index
    character(len=:), allocatable :: zaid_str
    integer :: zaid_str_len ! number of characters in zaid_str
    logical :: file_exists ! does urr.xml file exist?
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: settings_node
    type(XMLNode) :: otf_xs_node
    type(XMLNode) :: prob_table_node
    type(XMLNode) :: point_xs_node
    type(XMLNode) :: isotope_node
    type(XMLNode) :: isotopes_node
    type(XMLNode), allocatable :: isotope_node_list(:)

    ! Check if urr.xml exists
    filename = trim(adjustl(path_input)) // "urr.xml"
    print*,filename
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! Since a urr.xml file is optional, no error is issued here
      URR_use_urr = .false.
      return
    else
      URR_use_urr = .true.
      continue
    end if

    call write_message("Reading URR XML file...")

    ! Parse urr.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Check that required settings element is present
    if (check_for_node(root, "settings")) then
      settings_node = root % child("settings")

      ! Check that a path to ENDF data is specified
      if (check_for_node(settings_node, "endf_6_filepath")) then
        call get_node_value(settings_node, "endf_6_filepath", temp_str_fixed_len)
        temp_str = trim(adjustl(temp_str_fixed_len))
        if (temp_str(len(temp_str):) /= '/') temp_str = temp_str // '/'
        URR_path_endf_files = temp_str
      else
        call fatal_error('Specify path to ENDF-6 data files for URR treatment&
             & via endf_6_filepath in urr.xml')
      end if

      ! Check that a path to averaged URR cross section values is specified
      if (check_for_node(settings_node, "avg_xs_filepath")) then
        call get_node_value(settings_node, "avg_xs_filepath", temp_str_fixed_len)
        temp_str = trim(adjustl(temp_str_fixed_len))
        if (temp_str(len(temp_str):) /= '/') temp_str = temp_str // '/'
        URR_path_avg_xs = temp_str
      else
        call fatal_error('Specify path to averaged URR cross section data files&
             & via avg_xs_filepath in urr.xml')
      end if

      ! Check which formalism for the on-the-fly URR treatment is specified
      if (check_for_node(settings_node, 'formalism')) then
        call get_node_value(settings_node, 'formalism', temp_str_fixed_len)

        select case (trim(adjustl(to_lower(temp_str_fixed_len))))
        case ('slbw')
          URR_formalism = URR_SLBW
        case ('mlbw')
          URR_formalism = URR_MLBW
        case ('mnbw')
          URR_formalism = URR_MNBW
          call fatal_error('MNBW formalism not yet supported for the URR')
        case ('reich-moore')
          URR_formalism = URR_REICH_MOORE
          if (run_mode /= MODE_PURXS)&
               call fatal_error('Reich-Moore formalism not yet supported for the URR')
        case default
          call fatal_error('Unrecognized resonance formalism given in urr.xml')
        end select
      else
        call fatal_error('No formalism element given in urr.xml')
      end if

      ! Determine which W function evaluation to use
      if (check_for_node(settings_node, "w_function_implementation")) then
        call get_node_value(settings_node, "w_function_implementation", temp_str_fixed_len)

        select case (trim(adjustl(to_lower(temp_str_fixed_len))))
        case ('mit_w')
          URR_faddeeva_method = URR_MIT_W
        case ('quick_w')
          URR_faddeeva_method = URR_QUICK_W
        case default
          call fatal_error('Unrecognized Faddeeva function evaluation method&
               & in urr.xml')
        end select
      else
        call fatal_error('No Faddeeva function evaluation method specified in urr.xml')
      end if

      ! Assign number of l-wave resonances to be used in xs calculations
      if (check_for_node(settings_node, 'num_s_wave')) then
        call get_node_value(settings_node, 'num_s_wave', URR_num_l_waves(1))
      else
        call fatal_error('Must specify number of s-wave resonances in urr.xml')
      end if
      if (check_for_node(settings_node, 'num_p_wave')) then
        call get_node_value(settings_node, 'num_p_wave', URR_num_l_waves(2))
      else
        call fatal_error('Must specify number of p-wave resonances in urr.xml')
      end if
      if (check_for_node(settings_node, 'num_d_wave')) then
        call get_node_value(settings_node, 'num_d_wave', URR_num_l_waves(3))
      else
        call fatal_error('Must specify number of d-wave resonances in urr.xml')
      end if
      if (check_for_node(settings_node, 'num_f_wave')) then
        call get_node_value(settings_node, 'num_f_wave', URR_num_l_waves(4))
      else
        call fatal_error('Must specify number of f-wave resonances in urr.xml')
      end if

      ! Check for resonance parameter energy dependence treatment
      if (check_for_node(settings_node, 'parameter_energy_dependence')) then
        call get_node_value(settings_node, "parameter_energy_dependence", temp_str_fixed_len)

        select case (trim(adjustl(to_lower(temp_str_fixed_len))))
        case ('neutron')
          URR_parameter_energy_dependence = URR_E_NEUTRON
        case ('resonance')
          URR_parameter_energy_dependence = URR_E_RESONANCE
        case default
          call fatal_error('Unrecognized resonance parameter energy dependence&
               & in urr.xml')
        end select
      else
        call fatal_error('No resonance parameter energy dependence given in&
             & urr.xml')
      end if

      ! Include resonance structure of competitive URR cross sections?
      if (check_for_node(settings_node, "competitive_structure")) then
        call get_node_value(settings_node, "competitive_structure", temp_str_fixed_len)
        if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false' &
             .or. trim(adjustl(temp_str_fixed_len)) == '0') then
          URR_competitive_structure = .false.
        else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'true' &
             .or. trim(adjustl(temp_str_fixed_len)) == '1') then
          URR_competitive_structure = .true.
        else
          call fatal_error('Modeling of competitive reaction cross section &
               &resonance structure must be true or false in urr.xml')
        end if
      else
        call fatal_error('Must specify whether or not to model resonance &
             &structure of competitive reaction cross section in urr.xml')
      end if

    else
      call fatal_error('No settings element present in urr.xml')
    end if

    ! Check for required node containing list of isotopes to apply the treatment to
    if (check_for_node(root, 'isotopes')) then
      isotopes_node = root % child('isotopes')

      URR_num_isotopes = 0
      if (check_for_node(isotopes_node, 'isotope')) then
        call get_node_list(isotopes_node, 'isotope', isotope_node_list)
        URR_num_isotopes = size(isotope_node_list)
        allocate(URR_isotopes(URR_num_isotopes))
        allocate(URR_endf_filenames(URR_num_isotopes))

        ! Read in ZAID values and ENDF filenames for the isotopes(s)
        do i = 1, URR_num_isotopes
          isotope_node = isotope_node_list(i)

          ! Check that an ENDF data file is given
          if (check_for_node(isotope_node, 'endf_6_file')) then
            call get_node_value(isotope_node, 'endf_6_file', temp_str_fixed_len)
            URR_endf_filenames(i) = trim(adjustl(temp_str_fixed_len))
          else
            call fatal_error('No ENDF-6 data file specified for isotope ' &
                 // trim(to_str(i)) // ' in urr.xml file')
          end if

          ! Check if a max energy is given
          if (check_for_node(isotope_node, 'max_energy')) then
            call get_node_value(isotope_node, 'max_energy', URR_isotopes(i) % max_E_urr)
          else
            URR_isotopes(i) % max_E_urr = URR_INF
          end if
        end do

      else
        call fatal_error('No isotope(s) specified in urr.xml')
      end if

    else
      call fatal_error('No isotopes element in urr.xml')
    end if

    ! Check that an on-the-fly URR treatment and/or average xs calc is specified
    if (check_for_node(root, "on_the_fly")&
         .or. check_for_node(root, "probability_tables")&
         .or. check_for_node(root, "pointwise")) then
      continue
    else
      call fatal_error('No on-the-fly, pointwise, or probability table &
           &calculation parameters given in urr.xml')
    end if

    ! Check for pointwise cross section calculation
    if (check_for_node(root, "pointwise")) then
      URR_xs_representation = URR_POINTWISE
      URR_num_urr_realizations = 1
      point_xs_node = root % child("pointwise")
      if (URR_xs_representation == URR_ON_THE_FLY .or. URR_xs_representation == URR_PROB_BANDS) then
        call fatal_error('Cannot represent URR cross sections as pointwise and&
             & on-the-fly and/or as probability bands')
      end if

      ! load previously-generated pointwise xs?
      if (check_for_node(point_xs_node, 'read_xs')) then
        call get_node_value(point_xs_node, 'read_xs', temp_str_fixed_len)
        if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'true'&
             .or. trim(adjustl(temp_str_fixed_len)) == '1') then
          URR_pregenerated_point_xs = .true.
          if (master) then
            call write_message('Loading pointwise URR cross sections from files. Any&
                 & pointwise xs generation parameters specified in urr.xml&
                 & are being ignored.')
          end if
        else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false'&
             .or. trim(adjustl(temp_str_fixed_len)) == '0') then
          URR_pregenerated_point_xs = .false.
        else
          call fatal_error('read_xs must be true or false in urr.xml')
        end if
      else
        call fatal_error('Must specify whether or not to load pointwise &
             &cross sections in urr.xml')
      end if

      if (check_for_node(point_xs_node, 'source')) then
        call get_node_value(point_xs_node, "source", temp_str_fixed_len)
        select case (trim(adjustl(to_lower(temp_str_fixed_len))))
        case ('reconstruction')
          URR_xs_source_pointwise = URR_RECONSTRUCTION
        case ('hdf5')
          URR_xs_source_pointwise = URR_HDF5
        case default
          call fatal_error('Unrecognized pointwise URR xs source.')
        end select
      else
        call fatal_error('No pointwise xs source specified in urr.xml')
      end if
      if (URR_xs_source_pointwise == URR_RECONSTRUCTION) then
        if (check_for_node(point_xs_node, 'filepath')) then
          call get_node_value(point_xs_node, 'filepath', temp_str_fixed_len)
          temp_str = trim(adjustl(temp_str_fixed_len))
          if (temp_str(len(temp_str):) /= '/') temp_str = temp_str // '/'
          URR_path_point_xs = temp_str
        else
          call fatal_error('Must specify path for read/write of pointwise xs files &
               &in urr.xml')
        end if
        if (check_for_node(point_xs_node, 'min_energy_spacing')) then
          call get_node_value(&
               point_xs_node, 'min_energy_spacing', URR_min_delta_E_pointwise)
        else
          call fatal_error('No min energy spacing [eV] for cross section &
               &reconstruction is given in urr.xml')
        end if
        if (check_for_node(point_xs_node, 'tolerance')) then
          call get_node_value(point_xs_node, 'tolerance', URR_rel_err_tolerance_pointwise)
        else
          call fatal_error('No max fractional error tolerance for cross section&
               & reconstruction is given in urr.xml')
        end if
        if (check_for_node(point_xs_node, 'temperature')) then
          call get_node_value(&
               point_xs_node, 'temperature', URR_temperature_pointwise)
        else
          call fatal_error('No temperature [K] for cross section &
               &reconstruction is given in urr.xml')
        end if
      end if
    end if

    ! Check if an on-the-fly treatment is specified
    if (check_for_node(root, "on_the_fly")) then
      if (URR_xs_representation == URR_POINTWISE .or. URR_xs_representation == URR_PROB_BANDS) then
        call fatal_error('Cannot represent URR cross sections on-the-fly and&
             & as pointwise and/or as probability bands')
      end if
      URR_xs_representation = URR_ON_THE_FLY
      otf_xs_node = root % child("on_the_fly")

      ! Check if a xs calculation frequency is specified
      if (check_for_node(otf_xs_node, 'frequency')) then
        call get_node_value(otf_xs_node, 'frequency', temp_str_fixed_len)

        select case (trim(adjustl(to_lower(temp_str_fixed_len))))
        case ('event')
          URR_realization_frequency = URR_EVENT

        case ('history')
          URR_realization_frequency = URR_HISTORY
          call fatal_error('History-based URR realizations not yet supported')

        case ('batch')
          URR_realization_frequency = URR_BATCH
          call fatal_error('Batch-based URR realizations not yet supported')

        case ('simulation')
          URR_realization_frequency = URR_SIMULATION

        case default
          call fatal_error('Unrecognized URR realization frequency in urr.xml')
        end select
      else
        call fatal_error('No URR realization frequency given in urr.xml')
      end if

      if (URR_realization_frequency /= URR_EVENT) then
        ! determine number of independent URR realizations to generate
        if (check_for_node(otf_xs_node, 'num_realizations')) then
          call get_node_value(otf_xs_node, 'num_realizations', URR_num_urr_realizations)
          if (URR_num_urr_realizations < 1) then
            call fatal_error('Number of independent URR realizations&
                 & must be greater than 0 in urr.xml')
          end if
        else
          call fatal_error('No number of independent URR realizations to &
               &generate is given in urr.xml')
        end if

        ! determine which realization to use if multiple are generated
        if (check_for_node(otf_xs_node, "realization_index")) then
          call get_node_value(otf_xs_node, "realization_index", URR_i_realization_user)
          if (URR_i_realization_user < 0) then
            call fatal_error('realization_index must be non-negative in &
                 &urr.xml - 0 for a random realization or a positive integer &
                 &equal to the index of a single, desired realization')
          end if
        else
          call fatal_error('No URR realization index given in urr.xml')
        end if
      else
        if (check_for_node(otf_xs_node, 'realizations')&
             .or. check_for_node(otf_xs_node, 'realization_index')) then
          call fatal_error('a new local URR realization is being generated&
               & at each event; do not specify number or index of realizations')
        end if
      end if
    end if

    ! Check for probability band calculations parameters
    if (check_for_node(root, "probability_tables")) then
      if (URR_xs_representation == URR_POINTWISE .or. URR_xs_representation == URR_ON_THE_FLY) then
        call fatal_error('Cannot represent URR cross sections as probability&
             & bands and as pointwise and/or on-the-fly')
      end if
      URR_xs_representation = URR_PROB_BANDS
      prob_table_node = root % child("probability_tables")

      ! check that a filepath for probability table file I/O is specified
      if (check_for_node(prob_table_node, 'filepath')) then
        call get_node_value(prob_table_node, 'filepath', temp_str_fixed_len)
        temp_str = trim(adjustl(temp_str_fixed_len))
        if (temp_str(len(temp_str):) /= '/') temp_str = temp_str // '/'
        URR_path_prob_tables = temp_str
      else
        call fatal_error('Must specify path to probability table files &
             &in urr.xml')
      end if

      ! load previously-generated tables?
      if (check_for_node(prob_table_node, 'read_tables')) then
        call get_node_value(prob_table_node, 'read_tables', temp_str_fixed_len)
        if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'true'&
             .or. trim(adjustl(temp_str_fixed_len)) == '1') then
          URR_pregenerated_prob_tables = .true.
          if (master) then
            call write_message('Loading probability tables from files. Any&
                 & probability table generation parameters specified in urr.xml&
                 & are being ignored.')
          end if
        else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false'&
             .or. trim(adjustl(temp_str_fixed_len)) == '0') then
          URR_pregenerated_prob_tables = .false.
        else
          call fatal_error('Loading probability tables must be true or false&
               & in urr.xml')
        end if
      else
        call fatal_error('Must specify whether or not to load probability &
             &tables in urr.xml')
      end if

      ! temperature interpolation scheme
      if (check_for_node(prob_table_node, 'temperature_interpolation_method')) then
        call get_node_value(prob_table_node, 'temperature_interpolation_method',&
             temp_str_fixed_len)
        select case(trim(adjustl(to_lower(temp_str_fixed_len))))
        case('histogram')
          URR_temperature_interp_scheme = URR_HISTOGRAM
        case('lin-lin')
          URR_temperature_interp_scheme = URR_LINEAR_LINEAR
        case('lin-log')
          URR_temperature_interp_scheme = URR_LINEAR_LOG
        case('log-lin')
          URR_temperature_interp_scheme = URR_LOG_LINEAR
        case('log-log')
          URR_temperature_interp_scheme = URR_LOG_LOG
        case('sqrt-lin')
          URR_temperature_interp_scheme = URR_SQRT_LINEAR
        case('sqrt-log')
          URR_temperature_interp_scheme = URR_SQRT_LOG
        case('statistical')
          URR_temperature_interp_scheme = URR_STATISTICAL
        case('low-neighbor')
          URR_temperature_interp_scheme = URR_LOW_NEIGHBOR
        case default
          call fatal_error('Unrecognized temperature interpolation scheme')
        end select
      else
        call fatal_error('Must specify temperature interpolation scheme in &
             &urr.xml')
      end if

      if (.not. URR_pregenerated_prob_tables) then
        ! number of bands
        if (check_for_node(prob_table_node, "num_bands")) then
          call get_node_value(prob_table_node, "num_bands", URR_num_bands_prob_tables)
        else
          call fatal_error('No number of probability table cross section bands&
               & given in urr.xml')
        end if

        ! temperatures
        if (check_for_node(prob_table_node, "temperature_grid")) then
          URR_num_temperatures_prob_tables = node_word_count(prob_table_node, "temperature_grid")
          allocate(URR_T_grid_prob_tables(URR_num_temperatures_prob_tables))
          call get_node_array(prob_table_node, "temperature_grid", URR_T_grid_prob_tables)
          if (URR_num_temperatures_prob_tables > 1) then
            do i = 2, URR_num_temperatures_prob_tables
              if (URR_T_grid_prob_tables(i) <= URR_T_grid_prob_tables(i-1)) then
                call fatal_error('Cross section &
                     & reconstruction temperatures must be strictly ascending.')
              end if
            end do
          end if
        else
          call fatal_error('No probability table temperatures given in urr.xml')
        end if

        ! min number of batches
        if (check_for_node(prob_table_node, "min_num_batches")) then
          call get_node_value(prob_table_node, "min_num_batches", URR_min_num_batches_prob_tables)
        else
          call fatal_error('No minimum batches for probability table&
               & calculation given in urr.xml')
        end if

        ! max number of batches
        if (check_for_node(prob_table_node, "max_num_batches")) then
          call get_node_value(prob_table_node, "max_num_batches", URR_max_num_batches_prob_tables)
        else
          call fatal_error('No maximum batches for probability table&
               & calculation given in urr.xml')
        end if

        ! histories per batch
        if (check_for_node(prob_table_node, "num_histories")) then
          call get_node_value(prob_table_node, "num_histories", URR_num_histories_prob_tables)
          if (mod(URR_num_histories_prob_tables, URR_num_bands_prob_tables) /= 0) then
            call fatal_error('Histories &
                 &per batch must be evenly divisible by number of probability bands')
          end if
        else
          call fatal_error('No number of realization histories for probability&
               & table calculation given in urr.xml')
        end if

        ! relative uncertainty tolerance (1sigma / mean)
        if (check_for_node(prob_table_node, "tolerance")) then
          call get_node_value(prob_table_node, "tolerance", URR_rel_err_tolerance_avg_xs)
        else
          call fatal_error('No relative uncertainty tolerance for probability&
               & table calculation given in urr.xml')
        end if

        ! energy spacing scheme
        if (check_for_node(prob_table_node, "energy_spacing")) then
          call get_node_value(prob_table_node, "energy_spacing", temp_str_fixed_len)
          if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'linear') then
            URR_E_grid_scheme_prob_tables = URR_LINEAR
            call fatal_error('Linear probability table energy spacing not yet&
                 & supported in urr.xml')
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'logarithmic') then
            URR_E_grid_scheme_prob_tables = URR_LOGARITHMIC
            if (check_for_node(prob_table_node, "num_energies")) then
              call get_node_value(prob_table_node, "num_energies", URR_num_energies_prob_tables)
            else
              call fatal_error('Must specify number of logarithmic probability&
                   & table energies')
            end if
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'endf') then
            URR_E_grid_scheme_prob_tables = URR_ENDF6
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'user') then
            URR_E_grid_scheme_prob_tables = URR_USER
            if (check_for_node(prob_table_node, "energy_grid")) then
              URR_num_energies_prob_tables = node_word_count(prob_table_node, "energy_grid")
              allocate(URR_E_grid_prob_tables(URR_num_energies_prob_tables))
              call get_node_array(prob_table_node, "energy_grid", URR_E_grid_prob_tables)
            else
              call fatal_error('No probability table energy grid values&
                   & given in urr.xml')
            end if
          else
            call fatal_error('Unrecognized energy spacing scheme for probability&
                 & table calculation given in urr.xml')
          end if
        else
          call fatal_error('No energy spacing scheme for probability&
               & table calculation given in urr.xml')
        end if

        ! background cross section component treatment
        URR_background_xs_treatment = URR_FALSE
        if (check_for_node(prob_table_node, "background_xs")) then
          call get_node_value(prob_table_node, "background_xs", temp_str_fixed_len)
          if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'endf_6') then
            URR_background_xs_treatment = URR_ENDFFILE
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false') then
            URR_background_xs_treatment = URR_FALSE
          else
            call fatal_error('Unrecognized background cross section source&
                 & in probability table calculation given in urr.xml')
          end if
        end if

        if (check_for_node(prob_table_node, "write_tables")) then
          call get_node_value(prob_table_node, "write_tables", temp_str_fixed_len)
          if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false' &
               .or. trim(adjustl(temp_str_fixed_len)) == '0') then
            URR_write_prob_tables = .false.
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'true' &
               .or. trim(adjustl(temp_str_fixed_len)) == '1') then
            URR_write_prob_tables = .true.
          else
            call fatal_error('write_tables must be true or false in urr.xml')
          end if
        else
          call fatal_error('Specify whether or not to write out URR probability&
               & tables to a file via write_tables in urr.xml')
        end if

        if (check_for_node(prob_table_node, "write_avg_xs")) then
          call get_node_value(prob_table_node, "write_avg_xs", temp_str_fixed_len)
          if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'false' &
               .or. trim(adjustl(temp_str_fixed_len)) == '0') then
            URR_write_avg_xs = .false.
          else if (trim(adjustl(to_lower(temp_str_fixed_len))) == 'true' &
               .or. trim(adjustl(temp_str_fixed_len)) == '1') then
            URR_write_avg_xs = .true.
          else
            call fatal_error('write_avg_xs must be true or false in urr.xml')
          end if
        else
          call fatal_error('Specify whether or not to write out averaged URR&
               & cross section values to a file via write_avg_xs in urr.xml')
        end if
      end if
    end if

    ! Close URR XML file
    call doc % clear()
#endif
  end subroutine read_urr_xml

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
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_library
    type(XMLNode), allocatable :: node_library_list(:)

    ! Check if cross_sections.xml exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find cross_sections.xml file
      call fatal_error("Cross sections XML file '" &
           // trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Reading cross sections XML file...")

    ! Parse cross_sections.xml file
    call doc % load_file(path_cross_sections)
    root = doc % document_element()

    if (check_for_node(root, "directory")) then
      ! Copy directory information if present
      call get_node_value(root, "directory", directory)
    else
      ! If no directory is listed in cross_sections.xml, by default select the
      ! directory in which the cross_sections.xml file resides
      i = index(path_cross_sections, "/", BACK=.true.)
      directory = path_cross_sections(1:i)
    end if

    ! Get node list of all <library>
    call get_node_list(root, "library", node_library_list)
    n_libraries = size(node_library_list)

    ! Allocate xs_listings array
    if (n_libraries == 0) then
      call fatal_error("No cross section libraries present in cross_sections.xml &
           &file!")
    else
      allocate(libraries(n_libraries))
    end if

    do i = 1, n_libraries
      ! Get pointer to ace table XML node
      node_library = node_library_list(i)

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
    call doc % clear()

  end subroutine read_ce_cross_sections_xml

  subroutine read_mg_cross_sections_header(libraries)
    type(Library), allocatable, intent(out) :: libraries(:)

    integer :: i           ! loop index
    integer :: n_libraries
    logical :: file_exists ! does mgxs.h5 exist?
    integer(HID_T) :: file_id
    real(8), allocatable :: rev_energy_bins(:)
    character(len=MAX_WORD_LEN), allocatable :: names(:)

    ! Check if MGXS Library exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find MGXS Library file
      call fatal_error("Cross sections HDF5 file '" &
           // trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Reading cross sections HDF5 file...", 5)

    ! Open file for reading
    file_id = file_open(path_cross_sections, 'r', parallel=.true.)

    if (attribute_exists(file_id, "energy_groups")) then
      ! Get neutron energy group count
      call read_attribute(num_energy_groups, file_id, "energy_groups")
    else
      call fatal_error("'energy_groups' attribute must exist!")
    end if

    if (attribute_exists(file_id, "delayed_groups")) then
      ! Get neutron delayed group count
      call read_attribute(num_delayed_groups, file_id, "delayed_groups")
    else
      num_delayed_groups = 0
    end if

    allocate(rev_energy_bins(num_energy_groups + 1))
    allocate(energy_bins(num_energy_groups + 1))

    if (attribute_exists(file_id, "group structure")) then
      ! Get neutron group structure
      call read_attribute(energy_bins, file_id, "group structure")
    else
      call fatal_error("'group structure' attribute must exist!")
    end if

    ! First reverse the order of energy_groups
    energy_bins = energy_bins(num_energy_groups + 1:1:-1)

    allocate(energy_bin_avg(num_energy_groups))
    do i = 1, num_energy_groups
      energy_bin_avg(i) = HALF * (energy_bins(i) + energy_bins(i + 1))
    end do

    ! Get the datasets present in the library
    call get_groups(file_id, names)
    n_libraries = size(names)

    ! Allocate libraries array
    if (n_libraries == 0) then
      call fatal_error("At least one MGXS data set must be present in &
                       &mgxs library file!")
    else
      allocate(libraries(n_libraries))
    end if

    do i = 1, n_libraries
      ! Get name of material
      allocate(libraries(i) % materials(1))
      libraries(i) % materials(1) = names(i)
    end do

    ! Close MGXS HDF5 file
    call file_close(file_id)

  end subroutine read_mg_cross_sections_header

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
! NORMALIZE_AO Normalize the nuclide atom percents
!===============================================================================

  subroutine normalize_ao()
    integer :: i               ! index in materials array
    integer :: j               ! index over nuclides in material
    real(8) :: sum_percent     ! summation
    real(8) :: awr             ! atomic weight ratio
    real(8) :: x               ! atom percent
    logical :: percent_in_atom ! nuclides specified in atom percent?
    logical :: density_in_atom ! density specified in atom/b-cm?

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
            awr = nuclides_MG(mat % nuclide(j)) % obj % awr
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

        ! Change density in g/cm^3 to atom/b-cm. Since all values are now in
        ! atom percent, the sum needs to be re-evaluated as 1/sum(x*awr)
        if (.not. density_in_atom) then
          sum_percent = ZERO
          do j = 1, mat % n_nuclides
            if (run_CE) then
              awr = nuclides(mat % nuclide(j)) % awr
            else
              awr = nuclides_MG(mat % nuclide(j)) % obj % awr
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

        ! Calculate density in g/cm^3.
        mat % density_gpcc = ZERO
        do j = 1, mat % n_nuclides
          if (run_CE) then
            awr = nuclides(mat % nuclide(j)) % awr
          else
            awr = ONE
          end if
          mat % density_gpcc = mat % density_gpcc &
               + mat % atom_density(j) * awr * MASS_NEUTRON / N_AVOGADRO
        end do
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
    logical :: found
    type(VectorInt) :: i_sab_tables
    type(VectorInt) :: i_sab_nuclides

    do i = 1, size(materials)
      ! Skip materials with no S(a,b) tables
      if (.not. allocated(materials(i) % i_sab_tables)) cycle

      associate (mat => materials(i))

        ASSIGN_SAB: do k = 1, size(mat % i_sab_tables)
          ! In order to know which nuclide the S(a,b) table applies to, we need
          ! to search through the list of nuclides for one which has a matching
          ! name
          found = .false.
          associate (sab => sab_tables(mat % i_sab_tables(k)))
            FIND_NUCLIDE: do j = 1, size(mat % nuclide)
              if (any(sab % nuclides == nuclides(mat % nuclide(j)) % name)) then
                call i_sab_tables % push_back(mat % i_sab_tables(k))
                call i_sab_nuclides % push_back(j)
                found = .true.
              end if
            end do FIND_NUCLIDE
          end associate

          ! Check to make sure S(a,b) table matched a nuclide
          if (.not. found) then
            call fatal_error("S(a,b) table " // trim(mat % &
                 sab_names(k)) // " did not match any nuclide on material " &
                 // trim(to_str(mat % id)))
          end if
        end do ASSIGN_SAB

        ! Update i_sab_tables and i_sab_nuclides
        deallocate(mat % i_sab_tables)
        m = i_sab_tables % size()
        allocate(mat % i_sab_tables(m))
        allocate(mat % i_sab_nuclides(m))
        mat % i_sab_tables(:) = i_sab_tables % data(1:m)
        mat % i_sab_nuclides(:) = i_sab_nuclides % data(1:m)

        ! Clear entries in vectors for next material
        call i_sab_tables % clear()
        call i_sab_nuclides % clear()

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

  subroutine read_ce_cross_sections(libraries, library_dict, nuc_temps, sab_temps)
    type(Library),   intent(in)      :: libraries(:)
    type(DictCharInt), intent(inout) :: library_dict
    type(VectorReal), intent(in)     :: nuc_temps(:)
    type(VectorReal), intent(in)     :: sab_temps(:)

    integer :: i, j
    integer :: i_library
    integer :: i_nuclide
    integer :: i_sab
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

    ! Read cross sections
    do i = 1, size(materials)
      do j = 1, size(materials(i) % names)
        name = materials(i) % names(j)

        if (.not. already_read % contains(name)) then
          i_library = library_dict % get_key(to_lower(name))
          i_nuclide = nuclide_dict % get_key(to_lower(name))

          call write_message('Reading ' // trim(name) // ' from ' // &
               trim(libraries(i_library) % path), 6)

          ! Open file and make sure version is sufficient
          file_id = file_open(libraries(i_library) % path, 'r')
          call check_data_version(file_id)

          ! Read nuclide data from HDF5
          group_id = open_group(file_id, name)
          call nuclides(i_nuclide) % from_hdf5(group_id, nuc_temps(i_nuclide), &
               temperature_method, temperature_tolerance, master)
          call close_group(group_id)
          call file_close(file_id)

          ! Assign resonant scattering data
          if (res_scat_on) call assign_0K_elastic_scattering(nuclides(i_nuclide))

          ! Determine if minimum/maximum energy for this nuclide is greater/less
          ! than the previous
          if (size(nuclides(i_nuclide) % grid) >= 1) then
            energy_min_neutron = max(energy_min_neutron, &
                 nuclides(i_nuclide) % grid(1) % energy(1))
            energy_max_neutron = min(energy_max_neutron, nuclides(i_nuclide) % &
                 grid(1) % energy(size(nuclides(i_nuclide) % grid(1) % energy)))
          end if

          ! Add name and alias to dictionary
          call already_read % add(name)

          ! Read multipole file into the appropriate entry on the nuclides array
          if (temperature_multipole) call read_multipole_data(i_nuclide)
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

          ! Open file and make sure version matches
          file_id = file_open(libraries(i_library) % path, 'r')
          call check_data_version(file_id)

          ! Read S(a,b) data from HDF5
          group_id = open_group(file_id, name)
          call sab_tables(i_sab) % from_hdf5(group_id, sab_temps(i_sab), &
               temperature_method, temperature_tolerance)
          call close_group(group_id)
          call file_close(file_id)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do
    end do

    ! Associate S(a,b) tables with specific nuclides
    call assign_sab_tables()

    ! Show which nuclide results in lowest energy for neutron transport
    do i = 1, size(nuclides)
      if (nuclides(i) % grid(1) % energy(size(nuclides(i) % grid(1) % energy)) &
           == energy_max_neutron) then
        call write_message("Maximum neutron transport energy: " // &
             trim(to_str(energy_max_neutron)) // " eV for " // &
             trim(adjustl(nuclides(i) % name)), 7)
        exit
      end if
    end do

    ! If the user wants multipole, make sure we found a multipole library.
    if (temperature_multipole) then
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
! ASSIGN_TEMPERATURES If any cells have undefined temperatures, try to find
! their temperatures from material or global default temperatures
!===============================================================================

  subroutine assign_temperatures(material_temps)
    real(8), intent(in) :: material_temps(:)

    integer :: i, j
    integer :: i_material

    do i = 1, n_cells
      ! Ignore non-normal cells and cells with defined temperature.
      if (cells(i) % material(1) == NONE) cycle
      if (cells(i) % sqrtkT(1) /= ERROR_REAL) cycle

      ! Set the number of temperatures equal to the number of materials.
      deallocate(cells(i) % sqrtkT)
      allocate(cells(i) % sqrtkT(size(cells(i) % material)))

      ! Check each of the cell materials for temperature data.
      do j = 1, size(cells(i) % material)
        ! Arbitrarily set void regions to 0K.
        if (cells(i) % material(j) == MATERIAL_VOID) then
          cells(i) % sqrtkT(j) = ZERO
          cycle
        end if

        ! Use material default or global default temperature
        i_material = material_dict % get_key(cells(i) % material(j))
        if (material_temps(i_material) /= ERROR_REAL) then
          cells(i) % sqrtkT(j) = sqrt(K_BOLTZMANN * &
               material_temps(i_material))
        else
          cells(i) % sqrtkT(j) = sqrt(K_BOLTZMANN * temperature_default)
        end if
      end do
    end do
  end subroutine assign_temperatures

!===============================================================================
! ASSIGN_0K_ELASTIC_SCATTERING
!===============================================================================

  subroutine assign_0K_elastic_scattering(nuc)
    type(Nuclide), intent(inout) :: nuc

    integer :: i
    real(8) :: xs_cdf_sum

    nuc % resonant = .false.
    if (allocated(res_scat_nuclides)) then
      ! If resonant nuclides were specified, check the list explicitly
      do i = 1, size(res_scat_nuclides)
        if (nuc % name == res_scat_nuclides(i)) then
          nuc % resonant = .true.

          ! Make sure nuclide has 0K data
          if (.not. allocated(nuc % energy_0K)) then
            call fatal_error("Cannot treat " // trim(nuc % name) // " as a &
                 &resonant scatterer because 0 K elastic scattering data is &
                 &not present.")
          end if

          exit
        end if
      end do
    else
      ! Otherwise, assume that any that have 0 K elastic scattering data are
      ! resonant
      nuc % resonant = allocated(nuc % energy_0K)
    end if

    if (nuc % resonant) then
      ! Build CDF for 0K elastic scattering
      xs_cdf_sum = ZERO
      allocate(nuc % xs_cdf(0:size(nuc % energy_0K)))
      nuc % xs_cdf(0) = ZERO

      associate (E => nuc % energy_0K, xs => nuc % elastic_0K)
        do i = 1, size(E) - 1
          ! Negative cross sections result in a CDF that is not monotonically
          ! increasing. Set all negative xs values to zero.
          if (xs(i) < ZERO) xs(i) = ZERO

          ! build xs cdf
          xs_cdf_sum = xs_cdf_sum + (sqrt(E(i))*xs(i) + sqrt(E(i+1))*xs(i+1))&
               / TWO * (E(i+1) - E(i))
          nuc % xs_cdf(i) = xs_cdf_sum
        end do
      end associate
    end if
  end subroutine assign_0K_elastic_scattering

!===============================================================================
! READ_MULTIPOLE_DATA checks for the existence of a multipole library in the
! directory and loads it using multipole_read
!===============================================================================

  subroutine read_multipole_data(i_table)

    integer, intent(in) :: i_table  ! index in nuclides/sab_tables

    logical :: file_exists                 ! Does multipole library exist?
    character(7) :: readable               ! Is multipole library readable?
    character(MAX_FILE_LEN) :: filename  ! Path to multipole xs library

    ! For the time being, and I know this is a bit hacky, we just assume
    ! that the file will be ZZZAAAmM.h5.
    associate (nuc => nuclides(i_table))

      if (nuc % metastable > 0) then
        filename = trim(path_multipole) // trim(zero_padded(nuc % Z, 3)) // &
             trim(zero_padded(nuc % A, 3)) // 'm' // &
             trim(to_str(nuc % metastable)) // ".h5"
      else
        filename = trim(path_multipole) // trim(zero_padded(nuc % Z, 3)) // &
             trim(zero_padded(nuc % A, 3)) // ".h5"
      end if

      ! Check if Multipole library exists and is readable
      inquire(FILE=filename, EXIST=file_exists, READ=readable)
      if (.not. file_exists) then
        nuc % mp_present = .false.
        return
      elseif (readable(1:3) == 'NO') then
        call fatal_error("Multipole library '" // trim(filename) // "' is not &
             &readable! Change file permissions with chmod command.")
      end if

      call write_message("Loading Multipole XS table: " // filename, 6)

      allocate(nuc % multipole)

      call multipole_read(filename, nuc % multipole, i_table)
      nuc % mp_present = .true.

    end associate

  end subroutine read_multipole_data

!===============================================================================
! CHECK_DATA_VERSION checks for the right version of nuclear data within HDF5
! files
!===============================================================================

  subroutine check_data_version(file_id)
    integer(HID_T), intent(in) :: file_id

    integer, allocatable :: version(:)

    if (attribute_exists(file_id, 'version')) then
      call read_attribute(version, file_id, 'version')
      if (version(1) /= HDF5_VERSION(1)) then
        call fatal_error("HDF5 data format uses version " // trim(to_str(&
             version(1))) // "." // trim(to_str(version(2))) // " whereas &
             &your installation of OpenMC expects version " // trim(to_str(&
             HDF5_VERSION(1))) // ".x data.")
      end if
    else
      call fatal_error("HDF5 data does not indicate a version. Your &
           &installation of OpenMC expects version " // trim(to_str(&
           HDF5_VERSION(1))) // ".x data.")
    end if
  end subroutine check_data_version

end module input_xml
