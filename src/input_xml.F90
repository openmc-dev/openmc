module input_xml

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: find
  use cmfd_input,       only: configure_cmfd
  use cmfd_header,      only: cmfd_mesh
  use constants
  use dict_header,      only: DictIntInt, DictCharInt, DictEntryCI
  use distribution_multivariate
  use distribution_univariate
  use endf,             only: reaction_name
  use error,            only: fatal_error, warning, write_message
  use geometry,         only: calc_offsets, maximum_levels, count_instance, &
                              neighbor_lists
  use geometry_header
  use hdf5_interface
  use list_header,      only: ListChar, ListInt, ListReal
  use material_header
  use mesh_header
  use message_passing
  use mgxs_data,        only: create_macro_xs, read_mgxs
  use mgxs_header
  use multipole,        only: multipole_read
  use nuclide_header
  use output,           only: title, header, print_plot
  use plot_header
  use random_lcg,       only: prn, openmc_set_seed
  use surface_header
  use set_header,       only: SetChar
  use settings
  use source_header
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with, tokenize, split_string, &
                              zero_padded, to_c_string
  use summary,          only: write_summary
  use tally
  use tally_header,     only: openmc_extend_tallies
  use tally_derivative_header
  use tally_filter_header
  use tally_filter
  use timer_header,     only: time_read_xs
  use trigger_header
  use volume_header
  use xml_interface

  implicit none
  save

  interface
    subroutine read_surfaces(node_ptr) bind(C, name='read_surfaces')
      use ISO_C_BINDING
      implicit none
      type(C_PTR) :: node_ptr
    end subroutine read_surfaces
  end interface

contains

!===============================================================================
! READ_INPUT_XML calls each of the separate subroutines for reading settings,
! geometry, materials, and tallies.
!===============================================================================

  subroutine read_input_xml()

    type(VectorReal), allocatable :: nuc_temps(:) ! List of T to read for each nuclide
    type(VectorReal), allocatable :: sab_temps(:) ! List of T to read for each S(a,b)
    real(8), allocatable    :: material_temps(:)

    call read_settings_xml()
    call read_cross_sections_xml()
    call read_materials_xml(material_temps)
    call read_geometry_xml()

    ! Set up neighbor lists, convert user IDs -> indices, assign temperatures
    call finalize_geometry(material_temps, nuc_temps, sab_temps)

    if (run_mode /= MODE_PLOTTING) then
      call time_read_xs % start()
      if (run_CE) then
        ! Read continuous-energy cross sections
        call read_ce_cross_sections(nuc_temps, sab_temps)
      else
        ! Create material macroscopic data for MGXS
        call read_mgxs()
        call create_macro_xs()
      end if
      call time_read_xs % stop()
    end if

    call read_tallies_xml()

    ! Initialize distribcell_filters
    call prepare_distribcell()

    if (cmfd_run) call configure_cmfd()

    if (run_mode == MODE_PLOTTING) then
      ! Read plots.xml if it exists
      call read_plots_xml()
      if (master .and. verbosity >= 5) call print_plot()

    else
      ! Normalize atom/weight percents
      call normalize_ao()

      ! Write summary information
      if (master .and. output_summary) call write_summary()

      ! Warn if overlap checking is on
      if (master .and. check_overlaps) &
           call warning("Cell overlap checking is ON.")
    end if

  end subroutine read_input_xml

  subroutine finalize_geometry(material_temps, nuc_temps, sab_temps)
    real(8), intent(in) :: material_temps(:)
    type(VectorReal),            allocatable, intent(out) :: nuc_temps(:)
    type(VectorReal),  optional, allocatable, intent(out) :: sab_temps(:)

    ! Perform some final operations to set up the geometry
    call adjust_indices()
    call count_instance(universes(root_universe))

    ! After reading input and basic geometry setup is complete, build lists of
    ! neighboring cells for efficient tracking
    call neighbor_lists()

    ! Assign temperatures to cells that don't have temperatures already assigned
    call assign_temperatures(material_temps)

    ! Determine desired txemperatures for each nuclide and S(a,b) table
    call get_temperatures(nuc_temps, sab_temps)

    ! Check to make sure there are not too many nested coordinate levels in the
    ! geometry since the coordinate list is statically allocated for performance
    ! reasons
    if (maximum_levels(universes(root_universe)) > MAX_COORD) then
      call fatal_error("Too many nested coordinate levels in the geometry. &
           &Try increasing the maximum number of coordinate levels by &
           &providing the CMake -Dmaxcoord= option.")
    end if

  end subroutine finalize_geometry

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
    integer(C_INT32_T) :: i_start, i_end
    integer(C_INT64_T) :: seed
    integer(C_INT) :: err
    integer, allocatable :: temp_int_array(:)
    integer :: n_tracks
    logical :: file_exists
    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: filename
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_mode
    type(XMLNode) :: node_cutoff
    type(XMLNode) :: node_entropy
    type(XMLNode) :: node_ufs
    type(XMLNode) :: node_sp
    type(XMLNode) :: node_output
    type(XMLNode) :: node_res_scat
    type(XMLNode) :: node_trigger
    type(XMLNode) :: node_vol
    type(XMLNode) :: node_tab_leg
    type(XMLNode), allocatable :: node_mesh_list(:)
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
      if (n_batches <= n_inactive) then
        call fatal_error("Number of active batches must be greater than zero.")
      elseif (n_inactive < 0) then
        call fatal_error("Number of inactive batches must be non-negative.")
      elseif (n_particles <= 0) then
        call fatal_error("Number of particles must be greater than zero.")
      end if
    end if

    ! Copy random number seed if specified
    if (check_for_node(root, "seed")) then
      call get_node_value(root, "seed", seed)
      call openmc_set_seed(seed)
    end if

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
      if (master) call warning("Ignoring number of threads.")
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

    ! Check if we want to write out source
    if (check_for_node(root, "write_initial_source")) then
      call get_node_value(root, "write_initial_source", write_initial_source)
    end if

    ! Read each source
    do i = 1, n
      call external_source(i) % from_xml(node_source_list(i), path_source)
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
        call fatal_error("Number of integers specified in 'track' is not &
             &divisible by 3.  Please provide 3 integers per particle to be &
             &tracked.")
      end if

      ! Allocate space and get list of tracks
      allocate(temp_int_array(n_tracks))
      call get_node_array(root, "track", temp_int_array)

      ! Reshape into track_identifiers
      allocate(track_identifiers(3, n_tracks/3))
      track_identifiers = reshape(temp_int_array, [3, n_tracks/3])
    end if

    ! Read meshes
    call get_node_list(root, "mesh", node_mesh_list)

    ! Check for user meshes and allocate
    n = size(node_mesh_list)
    if (n > 0) then
      err = openmc_extend_meshes(n, i_start, i_end)
    end if

    do i = 1, n
      associate (m => meshes(i_start + i - 1))
        ! Instantiate mesh from XML node
        call m % from_xml(node_mesh_list(i))

        ! Add mesh to dictionary
        call mesh_dict % set(m % id, i_start + i - 1)
      end associate
    end do

    ! Shannon Entropy mesh
    if (check_for_node(root, "entropy_mesh")) then
      call get_node_value(root, "entropy_mesh", temp_int)
      if (mesh_dict % has(temp_int)) then
        index_entropy_mesh = mesh_dict % get(temp_int)
      else
        call fatal_error("Mesh " // to_str(temp_int) // " specified for &
             &Shannon entropy does not exist.")
      end if
    elseif (check_for_node(root, "entropy")) then
      call warning("Specifying a Shannon entropy mesh via the <entropy> element &
           &is deprecated. Please create a mesh using <mesh> and then reference &
           &it by specifying its ID in an <entropy_mesh> element.")

      ! Get pointer to entropy node
      node_entropy = root % child("entropy")

      err = openmc_extend_meshes(1, index_entropy_mesh)

      associate (m => meshes(index_entropy_mesh))
        ! Assign ID
        m % id = 10000

        call m % from_xml(node_entropy)
      end associate
    end if

    if (index_entropy_mesh > 0) then
      associate(m => meshes(index_entropy_mesh))
        if (.not. allocated(m % dimension)) then
          ! If the user did not specify how many mesh cells are to be used in
          ! each direction, we automatically determine an appropriate number of
          ! cells
          m % n_dimension = 3
          allocate(m % dimension(3))
          m % dimension = ceiling((n_particles/20)**(ONE/THREE))

          ! Calculate width
          m % width = (m % upper_right - m % lower_left) / m % dimension
        end if

        ! Allocate space for storing number of fission sites in each mesh cell
        allocate(entropy_p(1, product(m % dimension)))
      end associate

      ! Turn on Shannon entropy calculation
      entropy_on = .true.
    end if

    ! Uniform fission source weighting mesh
    if (check_for_node(root, "ufs_mesh")) then
      call get_node_value(root, "ufs_mesh", temp_int)
      if (mesh_dict % has(temp_int)) then
        index_ufs_mesh = mesh_dict % get(temp_int)
      else
        call fatal_error("Mesh " // to_str(temp_int) // " specified for &
             &uniform fission site method does not exist.")
      end if
    elseif (check_for_node(root, "uniform_fs")) then
      call warning("Specifying a UFS mesh via the <uniform_fs> element &
           &is deprecated. Please create a mesh using <mesh> and then reference &
           &it by specifying its ID in a <ufs_mesh> element.")

      ! Get pointer to ufs node
      node_ufs = root % child("uniform_fs")

      err = openmc_extend_meshes(1, index_ufs_mesh)

      ! Allocate mesh object and coordinates on mesh
      associate (m => meshes(index_ufs_mesh))
        ! Assign ID
        m % id = 10001

        call m % from_xml(node_ufs)
      end associate
    end if

    if (index_ufs_mesh > 0) then
      ! Allocate array to store source fraction for UFS
      allocate(source_frac(1, product(meshes(index_ufs_mesh) % dimension)))

      ! Turn on uniform fission source weighting
      ufs = .true.
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
    if (check_for_node(root, "temperature_range")) then
      call get_node_array(root, "temperature_range", temperature_range)
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

      ! Preallocate space for keff and entropy by generation
      call k_generation % reserve(n_max_batches*gen_per_batch)
      call entropy % reserve(n_max_batches*gen_per_batch)

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

  end subroutine get_run_parameters

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i, j, k, m, i_x, i_a, input_index
    integer :: n, n_mats, n_x, n_y, n_z, n_rings, n_rlats, n_hlats
    integer :: id
    integer :: univ_id
    integer :: n_cells_in_univ
    integer, allocatable :: temp_int_array(:)
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    character(:), allocatable :: region_spec
    type(Cell),     pointer :: c
    class(Lattice), pointer :: lat
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_cell
    type(XMLNode) :: node_lat
    type(XMLNode), allocatable :: node_cell_list(:)
    type(XMLNode), allocatable :: node_rlat_list(:)
    type(XMLNode), allocatable :: node_hlat_list(:)
    type(VectorInt) :: tokens
    type(VectorInt) :: rpn
    type(VectorInt) :: fill_univ_ids ! List of fill universe IDs
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains

    ! Display output message
    call write_message("Reading geometry XML file...", 5)

    ! Check if geometry.xml exists
    filename = trim(path_input) // "geometry.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Geometry XML file '" // trim(filename) // "' does not &
           &exist!")
    end if

    ! Parse geometry.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! ==========================================================================
    ! READ SURFACES FROM GEOMETRY.XML

    ! This variable is used to check whether at least one boundary condition was
    ! applied to a surface
    boundary_exists = .false.

    call read_surfaces(root % ptr)

    ! Allocate surfaces array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      surfaces(i) % ptr = surface_pointer_c(i - 1);

      if (surfaces(i) % bc() /= BC_TRANSMIT) boundary_exists = .true.

      ! Add surface to dictionary
      call surface_dict % set(surfaces(i) % id(), i)
    end do

    ! Check to make sure a boundary condition was applied to at least one
    ! surface
    if (run_mode /= MODE_PLOTTING) then
      if (.not. boundary_exists) then
        call fatal_error("No boundary conditions were applied to any surfaces!")
      end if
    end if

    ! ==========================================================================
    ! READ CELLS FROM GEOMETRY.XML

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
      if (cell_dict % has(c % id)) then
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

        ! Convert user IDs to surface indices
        do j = 1, tokens % size()
          id = tokens % data(j)
          if (id < OP_UNION) then
            if (surface_dict % has(abs(id))) then
              k = surface_dict % get(abs(id))
              tokens % data(j) = sign(k, id)
            end if
          end if
        end do

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
      call cell_dict % set(c % id, i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      univ_id = c % universe
      if (.not. cells_in_univ_dict % has(univ_id)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % set(univ_id, n_universes)
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get(univ_id)
      end if
      call cells_in_univ_dict % set(univ_id, n_cells_in_univ)

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
      if (lattice_dict % has(lat % id)) then
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
      call lattice_dict % set(lat % id, i)

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
      if (lattice_dict % has(lat % id)) then
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
      call lattice_dict % set(lat % id, n_rlats + i)

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
        n_cells_in_univ = cells_in_univ_dict % get(u % id)
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
      j = universe_dict % get(cells(i) % universe)

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

!===============================================================================
! READ_MATERIAL_XML reads data from a materials.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_cross_sections_xml()
    integer :: i, j
    logical                 :: file_exists
    character(MAX_FILE_LEN) :: env_variable
    character(MAX_LINE_LEN) :: filename
    type(XMLDocument)       :: doc
    type(XMLNode)           :: root

    ! Check if materials.xml exists
    filename = trim(path_input) // "materials.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Material XML file '" // trim(filename) // "' does not &
           &exist!")
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
      call read_ce_cross_sections_xml()
    else
      call read_mg_cross_sections_header()
    end if

    ! Creating dictionary that maps the name of the material to the entry
    do i = 1, size(libraries)
      do j = 1, size(libraries(i) % materials)
        call library_dict % set(to_lower(libraries(i) % materials(j)), i)
      end do
    end do

    ! Check that 0K nuclides are listed in the cross_sections.xml file
    if (allocated(res_scat_nuclides)) then
      do i = 1, size(res_scat_nuclides)
        if (.not. library_dict % has(to_lower(res_scat_nuclides(i)))) then
          call fatal_error("Could not find resonant scatterer " &
               // trim(res_scat_nuclides(i)) // " in cross_sections.xml file!")
        end if
      end do
    end if

  end subroutine read_cross_sections_xml

  subroutine read_materials_xml(material_temps)
    real(8), allocatable, intent(out) :: material_temps(:)

    integer :: i              ! loop index for materials
    integer :: j              ! loop index for nuclides
    integer :: k              ! loop index
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
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    real(8)                 :: val          ! value entered for density
    real(8)                 :: temp_dble    ! temporary double prec. real
    logical                 :: sum_density  ! density is sum of nuclide densities
    type(VectorChar)        :: names        ! temporary list of nuclide names
    type(VectorChar)        :: list_iso_lab ! temporary list of isotropic lab scatterers
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

    ! Display output message
    call write_message("Reading materials XML file...", 5)

    ! Check if materials.xml exists
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
      if (material_dict % has(mat % id)) then
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

      ! =======================================================================
      ! READ AND PARSE <isotropic> element

      if (check_for_node(node_mat, "isotropic")) then
        n = node_word_count(node_mat, "isotropic")
        allocate(sarray(n))
        call get_node_array(node_mat, "isotropic", sarray)
        do j = 1, n
          call list_iso_lab % push_back(sarray(j))
        end do
        deallocate(sarray)
      end if

      ! ========================================================================
      ! COPY NUCLIDES TO ARRAYS IN MATERIAL

      ! allocate arrays in Material object
      n = names % size()
      mat % n_nuclides = n
      allocate(mat % names(n))
      allocate(mat % nuclide(n))
      allocate(mat % atom_density(n))

      ALL_NUCLIDES: do j = 1, mat % n_nuclides
        ! Check that this nuclide is listed in the cross_sections.xml file
        name = trim(names % data(j))
        if (.not. library_dict % has(to_lower(name))) then
          call fatal_error("Could not find nuclide " // trim(name) &
               // " in cross_sections data file!")
        end if
        i_library = library_dict % get(to_lower(name))

        if (run_CE) then
          ! Check to make sure cross-section is continuous energy neutron table
          if (libraries(i_library) % type /= LIBRARY_NEUTRON) then
            call fatal_error("Cross-section table " // trim(name) &
                 // " is not a continuous-energy neutron table.")
          end if
        end if

        ! If this nuclide hasn't been encountered yet, we need to add its name
        ! and alias to the nuclide_dict
        if (.not. nuclide_dict % has(to_lower(name))) then
          index_nuclide    = index_nuclide + 1
          mat % nuclide(j) = index_nuclide

          call nuclide_dict % set(to_lower(name), index_nuclide)
        else
          mat % nuclide(j) = nuclide_dict % get(to_lower(name))
        end if

        ! Copy name and atom/weight percent
        mat % names(j) = name
        mat % atom_density(j) = densities % data(j)

      end do ALL_NUCLIDES

      if (run_CE) then
        ! By default, isotropic-in-lab is not used
        if (list_iso_lab % size() > 0) then
          mat % has_isotropic_nuclides = .true.
          allocate(mat % p0(n))
          mat % p0(:) = .false.

          ! Apply isotropic-in-lab treatment to specified nuclides
          do j = 1, list_iso_lab % size()
            do k = 1, n
              if (names % data(k) == list_iso_lab % data(j)) then
                mat % p0(k) = .true.
              end if
            end do
          end do
        end if
      end if

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
          allocate(mat % sab_fracs(n_sab))

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

            ! Read the fraction of nuclei affected by this S(a,b) table
            if (check_for_node(node_sab, "fraction")) then
              call get_node_value(node_sab, "fraction", mat % sab_fracs(j))
            else
              mat % sab_fracs(j) = ONE
            end if

            ! Check that this nuclide is listed in the cross_sections.xml file
            if (.not. library_dict % has(to_lower(name))) then
              call fatal_error("Could not find S(a,b) table " // trim(name) &
                   // " in cross_sections.xml file!")
            end if

            ! Find index in xs_listing and set the name and alias according to the
            ! listing
            i_library = library_dict % get(to_lower(name))

            if (run_CE) then
              ! Check to make sure cross-section is continuous energy neutron table
              if (libraries(i_library) % type /= LIBRARY_THERMAL) then
                call fatal_error("Cross-section table " // trim(name) &
                     // " is not a S(a,b) table.")
              end if
            end if

            ! If this S(a,b) table hasn't been encountered yet, we need to add its
            ! name and alias to the sab_dict
            if (.not. sab_dict % has(to_lower(name))) then
              index_sab = index_sab + 1
              mat % i_sab_tables(j) = index_sab
              call sab_dict % set(to_lower(name), index_sab)
            else
              mat % i_sab_tables(j) = sab_dict % get(to_lower(name))
            end if
          end do
        end if
      end if

      ! Add material to dictionary
      call material_dict % set(mat % id, i)
    end do

    ! Set total number of nuclides and S(a,b) tables
    n_nuclides = index_nuclide
    n_sab_tables = index_sab

    ! Close materials XML file
    call doc % clear()

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
    integer :: filter_id     ! user-specified identifier for filter
    integer :: i_filt        ! index in filters array
    integer :: i_filter_mesh ! index of mesh filter
    integer :: i_elem        ! index of entry in dictionary
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filter      ! number of filters
    integer :: n_new         ! number of new scores to add based on Yn/Pn tally
    integer :: n_scores      ! number of tot scores after adjusting for Yn/Pn tally
    integer :: n_bins        ! total new bins for this score
    integer :: n_user_trig   ! number of user-specified tally triggers
    integer :: trig_ind      ! index of triggers array for each tally
    integer :: user_trig_ind ! index of user-specified triggers for each tally
    integer :: i_start, i_end
    integer :: i_filt_start, i_filt_end
    integer(C_INT) :: err
    real(8) :: threshold     ! trigger convergence threshold
    integer :: n_order       ! moment order requested
    integer :: n_order_pos   ! oosition of Scattering order in score name string
    integer :: MT            ! user-specified MT for score
    integer :: imomstr       ! Index of MOMENT_STRS & MOMENT_N_STRS
    logical :: file_exists   ! does tallies.xml file exist?
    integer, allocatable :: temp_filter(:) ! temporary filter indices
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    type(DictCharInt) :: trigger_scores
    type(TallyFilterContainer), pointer :: f
    type(RegularMesh), pointer :: m
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_tal
    type(XMLNode) :: node_filt
    type(XMLNode) :: node_trigger
    type(XMLNode), allocatable :: node_mesh_list(:)
    type(XMLNode), allocatable :: node_tal_list(:)
    type(XMLNode), allocatable :: node_filt_list(:)
    type(XMLNode), allocatable :: node_trigger_list(:)
    type(XMLNode), allocatable :: node_deriv_list(:)
    type(DictEntryCI) :: elem

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

    ! Display output message
    call write_message("Reading tallies XML file...", 5)

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

    ! Check for <assume_separate> setting
    if (check_for_node(root, "assume_separate")) then
      call get_node_value(root, "assume_separate", assume_separate)
    end if

    ! ==========================================================================
    ! READ MESH DATA

    ! Check for user meshes and allocate
    n = size(node_mesh_list)
    if (n > 0) then
      err = openmc_extend_meshes(n, i_start, i_end)
    end if

    do i = 1, n
      m => meshes(i_start + i - 1)

      ! Instantiate mesh from XML node
      call m % from_xml(node_mesh_list(i))

      ! Add mesh to dictionary
      call mesh_dict % set(m % id, i_start + i - 1)
    end do

    ! We only need the mesh info for plotting
    if (run_mode == MODE_PLOTTING) then
      call doc % clear()
      return
    end if

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
      call tally_derivs(i) % from_xml(node_deriv_list(i))

      ! Update tally derivative dictionary
      call tally_deriv_dict % set(tally_derivs(i) % id, i)
    end do

    ! ==========================================================================
    ! READ FILTER DATA

    ! Check for user filters and allocate
    n = size(node_filt_list)
    if (n > 0) then
      err = openmc_extend_filters(n, i_start, i_end)
    end if

    READ_FILTERS: do i = 1, n
      f => filters(i_start + i - 1)

      ! Get pointer to filter xml node
      node_filt = node_filt_list(i)

      ! Copy filter id
      if (check_for_node(node_filt, "id")) then
        call get_node_value(node_filt, "id", filter_id)
      else
        call fatal_error("Must specify id for filter in tally XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (filter_dict % has(filter_id)) then
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
      case ("mesh", "universe", "material", "cell", "distribcell", &
            "cellborn", "cellfrom", "surface", "delayedgroup")
        if (.not. check_for_node(node_filt, "bins")) then
          call fatal_error("Bins not set in filter " // trim(to_str(filter_id)))
        end if
      end select

      ! Allocate according to the filter type
      err = openmc_filter_set_type(i_start + i - 1, to_c_string(temp_str))

      ! Read filter data from XML
      call f % obj % from_xml(node_filt)

      ! Set filter id
      err = openmc_filter_set_id(i_start + i - 1, filter_id)

      ! Initialize filter
      call f % obj % initialize()
    end do READ_FILTERS

    ! ==========================================================================
    ! READ TALLY DATA

    ! Check for user tallies
    n = size(node_tal_list)
    if (n == 0) then
      if (master) call warning("No tallies present in tallies.xml file!")
    end if

    ! Allocate user tallies
    if (n > 0 .and. run_mode /= MODE_PLOTTING) then
      err = openmc_extend_tallies(n, i_start, i_end)
    end if

    READ_TALLIES: do i = 1, n
      ! Allocate tally
      err = openmc_tally_set_type(i_start + i - 1, &
           C_CHAR_'generic' // C_NULL_CHAR)

      ! Get pointer to tally
      associate (t => tallies(i_start + i - 1) % obj)

      ! Get pointer to tally xml node
      node_tal = node_tal_list(i)

      ! Copy tally id
      if (check_for_node(node_tal, "id")) then
        call get_node_value(node_tal, "id", t % id)
      else
        call fatal_error("Must specify id for tally in tally XML file.")
      end if

      ! Check to make sure 'id' hasn't been used
      if (tally_dict % has(t % id)) then
        call fatal_error("Two or more tallies use the same unique ID: " &
             // to_str(t % id))
      end if

      ! Copy tally name
      if (check_for_node(node_tal, "name")) &
           call get_node_value(node_tal, "name", t % name)

      ! =======================================================================
      ! READ DATA FOR FILTERS

      ! Check if user is using old XML format and throw an error if so
      if (check_for_node(node_tal, "filter")) then
        call fatal_error("Tally filters must be specified independently of &
             &tallies in a <filter> element. The <tally> element itself should &
             &have a list of filters that apply, e.g., <filters>1 2</filters> &
             &where 1 and 2 are the IDs of filters specified outside of &
             &<tally>.")
      end if

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

        do j = 1, n_filter
          ! Get pointer to filter
          if (filter_dict % has(temp_filter(j))) then
            i_filt = filter_dict % get(temp_filter(j))
            f => filters(i_filt)
          else
            call fatal_error("Could not find filter " &
                 // trim(to_str(temp_filter(j))) // " specified on tally " &
                 // trim(to_str(t % id)))
          end if

          ! Store the index of the filter
          temp_filter(j) = i_filt
        end do

        ! Set the filters
        err = openmc_tally_set_filters(i_start + i - 1, n_filter, temp_filter)
      end if
      deallocate(temp_filter)

      ! =======================================================================
      ! READ DATA FOR NUCLIDES

      if (check_for_node(node_tal, "nuclides")) then

        ! Allocate a temporary string array for nuclides and copy values over
        allocate(sarray(node_word_count(node_tal, "nuclides")))
        call get_node_array(node_tal, "nuclides", sarray)

        if (trim(sarray(1)) == 'all') then
          ! Handle special case <nuclides>all</nuclides>
          allocate(t % nuclide_bins(n_nuclides + 1))

          ! Set bins to 1, 2, 3, ..., n_nuclides, -1
          t % nuclide_bins(1:n_nuclides) = &
               (/ (j, j=1, n_nuclides) /)
          t % nuclide_bins(n_nuclides + 1) = -1

          ! Set number of nuclide bins
          t % n_nuclide_bins = n_nuclides + 1

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
            if (.not. nuclide_dict % has(word)) then
              call fatal_error("Could not find the nuclide " &
                   // trim(word) // " specified in tally " &
                   // trim(to_str(t % id)) // " in any material.")
            end if

            ! Set bin to index in nuclides array
            t % nuclide_bins(j) = nuclide_dict % get(word)
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
          if (trigger_on) call trigger_scores % set(trim(score_name), j)

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
            t % depletion_rx = .true.

          case ('n3n', '(n,3n)')
            t % score_bins(j) = N_3N
            t % depletion_rx = .true.

          case ('n4n', '(n,4n)')
            t % score_bins(j) = N_4N
            t % depletion_rx = .true.

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

            ! Check which type of current is desired: mesh currents or
            ! surface currents
            if (t % find_filter(FILTER_SURFACE) > 0 .or. &
                 &t % find_filter(FILTER_CELL) > 0 .or. &
                 &t % find_filter(FILTER_CELLFROM) > 0) then

              ! Check to make sure that mesh currents are not desired as well
              if (t % find_filter(FILTER_MESH) > 0) then
                call fatal_error("Cannot tally other mesh currents &
                     &in the same tally as surface currents")
              end if

              t % type = TALLY_SURFACE
              t % score_bins(j) = SCORE_CURRENT

            else if (t % find_filter(FILTER_MESH) > 0) then
              t % score_bins(j) = SCORE_CURRENT
              t % type = TALLY_MESH_CURRENT

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
                m => meshes(filt % mesh)
              end select

              ! Extend the filters array so we can add a surface filter and
              ! mesh filter
              err = openmc_extend_filters(2, i_filt_start, i_filt_end)

              ! Duplicate the mesh filter since other tallies might use this
              ! filter and we need to change the dimension
              filters(i_filt_start) = filters(i_filter_mesh)

              ! We need to increase the dimension by one since we also need
              ! currents coming into and out of the boundary mesh cells.
              filters(i_filt_start) % obj % n_bins = product(m % dimension + 1)

              ! Set ID
              call openmc_get_filter_next_id(filter_id)
              err = openmc_filter_set_id(i_filt_start, filter_id)


              ! Add surface filter
              allocate(SurfaceFilter :: filters(i_filt_end) % obj)
              select type (filt => filters(i_filt_end) % obj)
              type is (SurfaceFilter)
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

                ! Set ID
                call openmc_get_filter_next_id(filter_id)
                err = openmc_filter_set_id(i_filt_end, filter_id)
              end select

              ! Copy filter indices to resized array
              n_filter = size(t % filter)
              allocate(temp_filter(n_filter + 1))
              temp_filter(1:size(t % filter)) = t % filter
              n_filter = n_filter + 1

              ! Set mesh and surface filters
              temp_filter(t % find_filter(FILTER_MESH)) = i_filt_start
              temp_filter(n_filter) = i_filt_end

              ! Set filters
              err = openmc_tally_set_filters(i_start + i - 1, n_filter, temp_filter)
              deallocate(temp_filter)
            end if

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
            t % depletion_rx = .true.
          case ('(n,p)')
            t % score_bins(j) = N_P
            t % depletion_rx = .true.
          case ('(n,d)')
            t % score_bins(j) = N_D
          case ('(n,t)')
            t % score_bins(j) = N_T
          case ('(n,3He)')
            t % score_bins(j) = N_3HE
          case ('(n,a)')
            t % score_bins(j) = N_A
            t % depletion_rx = .true.
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
              t % n_triggers = t % n_triggers + trigger_scores % size()
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

              ! Loop over all tally scores
              i_elem = 0
              do
                ! Move to next score
                call trigger_scores % next_entry(elem, i_elem)
                if (i_elem == 0) exit

                score_name = trim(elem % key)

                ! Store the score name and index in the trigger
                t % triggers(trig_ind) % score_name = score_name
                t % triggers(trig_ind) % score_index = elem % value

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

                ! Increment the overall trigger index
                trig_ind = trig_ind + 1
              end do

            ! Scores other than the "all" placeholder
            else

              ! Store the score name and index
              t % triggers(trig_ind) % score_name = trim(score_name)
              t % triggers(trig_ind) % score_index = &
                   trigger_scores % get(trim(score_name))

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
      call tally_dict % set(t % id, i)

      end associate
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

    ! Display output message
    call write_message("Reading plot XML file...", 5)

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
      if (plot_dict % has(pl % id)) then
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
          if (master) call warning("Color specifications ignored in voxel &
               &plot " // trim(to_str(pl % id)))
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
            call fatal_error("Must specify id for color specification in &
                 &plot " // trim(to_str(pl % id)))
          end if

          ! Add RGB
          if (pl % color_by == PLOT_COLOR_CELLS) then

            if (cell_dict % has(col_id)) then
              col_id = cell_dict % get(col_id)
              call get_node_array(node_col, "rgb", pl % colors(col_id) % rgb)
            else
              call fatal_error("Could not find cell " // trim(to_str(col_id)) &
                   // " specified in plot " // trim(to_str(pl % id)))
            end if

          else if (pl % color_by == PLOT_COLOR_MATS) then

            if (material_dict % has(col_id)) then
              col_id = material_dict % get(col_id)
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

              if (index_ufs_mesh < 0) then
                call fatal_error("No UFS mesh for meshlines on plot " &
                     // trim(to_str(pl % id)))
              end if

              pl % meshlines_mesh => meshes(index_ufs_mesh)

            case ('cmfd')

              if (.not. cmfd_run) then
                call fatal_error("Need CMFD run to plot CMFD mesh for &
                     &meshlines on plot " // trim(to_str(pl % id)))
              end if

              pl % meshlines_mesh => cmfd_mesh

            case ('entropy')

              if (index_entropy_mesh < 0) then
                call fatal_error("No entropy mesh for meshlines on plot " &
                     // trim(to_str(pl % id)))
              end if

              pl % meshlines_mesh => meshes(index_entropy_mesh)

            case ('tally')

              ! Ensure that there is a mesh id if the type is tally
              if (check_for_node(node_meshlines, "id")) then
                call get_node_value(node_meshlines, "id", meshid)
              else
                call fatal_error("Must specify a mesh id for meshlines tally &
                     &mesh specification in plot " // trim(to_str(pl % id)))
              end if

              ! Check if the specified tally mesh exists
              if (mesh_dict % has(meshid)) then
                pl % meshlines_mesh => meshes(mesh_dict % get(meshid))
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

                if (cell_dict % has(col_id)) then
                  iarray(j) = cell_dict % get(col_id)
                else
                  call fatal_error("Could not find cell " &
                       // trim(to_str(col_id)) // " specified in the mask in &
                       &plot " // trim(to_str(pl % id)))
                end if

              else if (pl % color_by == PLOT_COLOR_MATS) then

                if (material_dict % has(col_id)) then
                  iarray(j) = material_dict % get(col_id)
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
      call plot_dict % set(pl % id, i)

    end do READ_PLOTS

    ! Close plots XML file
    call doc % clear()

  end subroutine read_plots_xml

!===============================================================================
! READ_*_CROSS_SECTIONS_XML reads information from a cross_sections.xml file. This
! file contains a listing of the CE and MG cross sections that may be used.
!===============================================================================

  subroutine read_ce_cross_sections_xml()
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

    call write_message("Reading cross sections XML file...", 5)

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

  subroutine read_mg_cross_sections_header()
    integer :: i           ! loop index
    integer :: n_libraries
    logical :: file_exists ! does mgxs.h5 exist?
    integer(HID_T) :: file_id
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
    rev_energy_bins = energy_bins
    energy_bins = energy_bins(num_energy_groups + 1:1:-1)

    ! Get the midpoint of the energy groups
    allocate(energy_bin_avg(num_energy_groups))
    do i = 1, num_energy_groups
      energy_bin_avg(i) = HALF * (energy_bins(i) + energy_bins(i + 1))
    end do

    ! Get the minimum and maximum energies
    energy_min_neutron = energy_bins(num_energy_groups + 1)
    energy_max_neutron = energy_bins(1)

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

  subroutine read_ce_cross_sections(nuc_temps, sab_temps)
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

    allocate(nuclides(n_nuclides))
    allocate(sab_tables(n_sab_tables))

    ! Read cross sections
    do i = 1, size(materials)
      do j = 1, size(materials(i) % names)
        name = materials(i) % names(j)

        if (.not. already_read % contains(name)) then
          i_library = library_dict % get(to_lower(name))
          i_nuclide = nuclide_dict % get(to_lower(name))

          call write_message('Reading ' // trim(name) // ' from ' // &
               trim(libraries(i_library) % path), 6)

          ! Open file and make sure version is sufficient
          file_id = file_open(libraries(i_library) % path, 'r')
          call check_data_version(file_id)

          ! Read nuclide data from HDF5
          group_id = open_group(file_id, name)
          call nuclides(i_nuclide) % from_hdf5(group_id, nuc_temps(i_nuclide), &
               temperature_method, temperature_tolerance, temperature_range, &
               master, i_nuclide)
          call close_group(group_id)
          call file_close(file_id)

          ! Assign resonant scattering data
          if (res_scat_on) call nuclides(i_nuclide) % assign_0K_elastic_scattering()

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

    ! Set up logarithmic grid for nuclides
    do i = 1, size(nuclides)
      call nuclides(i) % init_grid(energy_min_neutron, &
           energy_max_neutron, n_log_bins)
    end do
    log_spacing = log(energy_max_neutron/energy_min_neutron) / n_log_bins

    do i = 1, size(materials)
      ! Skip materials with no S(a,b) tables
      if (.not. allocated(materials(i) % sab_names)) cycle

      do j = 1, size(materials(i) % sab_names)
        ! Get name of S(a,b) table
        name = materials(i) % sab_names(j)

        if (.not. already_read % contains(name)) then
          i_library = library_dict % get(to_lower(name))
          i_sab  = sab_dict % get(to_lower(name))

          call write_message('Reading ' // trim(name) // ' from ' // &
               trim(libraries(i_library) % path), 6)

          ! Open file and make sure version matches
          file_id = file_open(libraries(i_library) % path, 'r')
          call check_data_version(file_id)

          ! Read S(a,b) data from HDF5
          group_id = open_group(file_id, name)
          call sab_tables(i_sab) % from_hdf5(group_id, sab_temps(i_sab), &
               temperature_method, temperature_tolerance, temperature_range)
          call close_group(group_id)
          call file_close(file_id)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do

      ! Associate S(a,b) tables with specific nuclides
      call materials(i) % assign_sab_tables()
    end do

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
        i_material = cells(i) % material(j)
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

      ! Display message
      call write_message("Loading Multipole XS table: " // filename, 6)

      allocate(nuc % multipole)

      ! Call the read routine
      call multipole_read(filename, nuc % multipole, i_table)
      nuc % mp_present = .true.

    end associate

  end subroutine read_multipole_data

!===============================================================================
! ADJUST_INDICES changes the values for 'surfaces' for each cell and the
! material index assigned to each to the indices in the surfaces and material
! array rather than the unique IDs assigned to each surface and material. Also
! assigns boundary conditions to surfaces based on those read into the bc_dict
! dictionary
!===============================================================================

  subroutine adjust_indices()

    integer :: i                      ! index for various purposes
    integer :: j                      ! index for various purposes
    integer :: k                      ! loop index for lattices
    integer :: m                      ! loop index for lattices
    integer :: lid                    ! lattice IDs
    integer :: id                     ! user-specified id
    class(Lattice),    pointer :: lat => null()

    do i = 1, n_cells
      ! =======================================================================
      ! ADJUST UNIVERSE INDEX FOR EACH CELL
      associate (c => cells(i))

      id = c % universe
      if (universe_dict % has(id)) then
        c % universe = universe_dict % get(id)
      else
        call fatal_error("Could not find universe " // trim(to_str(id)) &
             &// " specified on cell " // trim(to_str(c % id)))
      end if

      ! =======================================================================
      ! ADJUST MATERIAL/FILL POINTERS FOR EACH CELL

      if (c % material(1) == NONE) then
        id = c % fill
        if (universe_dict % has(id)) then
          c % type = FILL_UNIVERSE
          c % fill = universe_dict % get(id)
        elseif (lattice_dict % has(id)) then
          lid = lattice_dict % get(id)
          c % type = FILL_LATTICE
          c % fill = lid
        else
          call fatal_error("Specified fill " // trim(to_str(id)) // " on cell "&
               // trim(to_str(c % id)) // " is neither a universe nor a &
               &lattice.")
        end if
      else
        do j = 1, size(c % material)
          id = c % material(j)
          if (id == MATERIAL_VOID) then
            c % type = FILL_MATERIAL
          else if (material_dict % has(id)) then
            c % type = FILL_MATERIAL
            c % material(j) = material_dict % get(id)
          else
            call fatal_error("Could not find material " // trim(to_str(id)) &
                 // " specified on cell " // trim(to_str(c % id)))
          end if
        end do
      end if
      end associate
    end do

    ! ==========================================================================
    ! ADJUST UNIVERSE INDICES FOR EACH LATTICE

    do i = 1, n_lattices
      lat => lattices(i) % obj
      select type (lat)

      type is (RectLattice)
        do m = 1, lat % n_cells(3)
          do k = 1, lat % n_cells(2)
            do j = 1, lat % n_cells(1)
              id = lat % universes(j,k,m)
              if (universe_dict % has(id)) then
                lat % universes(j,k,m) = universe_dict % get(id)
              else
                call fatal_error("Invalid universe number " &
                     &// trim(to_str(id)) // " specified on lattice " &
                     &// trim(to_str(lat % id)))
              end if
            end do
          end do
        end do

      type is (HexLattice)
        do m = 1, lat % n_axial
          do k = 1, 2*lat % n_rings - 1
            do j = 1, 2*lat % n_rings - 1
              if (j + k < lat % n_rings + 1) then
                cycle
              else if (j + k > 3*lat % n_rings - 1) then
                cycle
              end if
              id = lat % universes(j, k, m)
              if (universe_dict % has(id)) then
                lat % universes(j, k, m) = universe_dict % get(id)
              else
                call fatal_error("Invalid universe number " &
                     &// trim(to_str(id)) // " specified on lattice " &
                     &// trim(to_str(lat % id)))
              end if
            end do
          end do
        end do

      end select

      if (lat % outer /= NO_OUTER_UNIVERSE) then
        if (universe_dict % has(lat % outer)) then
          lat % outer = universe_dict % get(lat % outer)
        else
          call fatal_error("Invalid universe number " &
               &// trim(to_str(lat % outer)) &
               &// " specified on lattice " // trim(to_str(lat % id)))
        end if
      end if

    end do

  end subroutine adjust_indices

!===============================================================================
! PREPARE_DISTRIBCELL initializes any distribcell filters present and sets the
! offsets for distribcells
!===============================================================================

  subroutine prepare_distribcell()

    integer :: i, j                ! Tally, filter loop counters
    logical :: distribcell_active  ! Does simulation use distribcell?
    integer, allocatable :: univ_list(:)              ! Target offsets
    integer, allocatable :: counts(:,:)               ! Target count
    logical, allocatable :: found(:,:)                ! Target found

    ! Assume distribcell is not needed until proven otherwise.
    distribcell_active = .false.

    ! We need distribcell if any tallies have distribcell filters.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % obj % filter)
        select type(filt => filters(tallies(i) % obj % filter(j)) % obj)
        type is (DistribcellFilter)
          distribcell_active = .true.
        end select
      end do
    end do

    ! We also need distribcell if any distributed materials or distributed
    ! temperatues are present.
    if (.not. distribcell_active) then
      do i = 1, n_cells
        if (size(cells(i) % material) > 1 .or. size(cells(i) % sqrtkT) > 1) then
          distribcell_active = .true.
          exit
        end if
      end do
    end if

    ! If distribcell isn't used in this simulation then no more work left to do.
    if (.not. distribcell_active) return

    ! Make sure the number of materials and temperatures matches the number of
    ! cell instances.
    do i = 1, n_cells
      associate (c => cells(i))
        if (size(c % material) > 1) then
          if (size(c % material) /= c % instances) then
            call fatal_error("Cell " // trim(to_str(c % id)) // " was &
                 &specified with " // trim(to_str(size(c % material))) &
                 // " materials but has " // trim(to_str(c % instances)) &
                 // " distributed instances. The number of materials must &
                 &equal one or the number of instances.")
          end if
        end if
        if (size(c % sqrtkT) > 1) then
          if (size(c % sqrtkT) /= c % instances) then
            call fatal_error("Cell " // trim(to_str(c % id)) // " was &
                 &specified with " // trim(to_str(size(c % sqrtkT))) &
                 // " temperatures but has " // trim(to_str(c % instances)) &
                 // " distributed instances. The number of temperatures must &
                 &equal one or the number of instances.")
          end if
        end if
      end associate
    end do

    ! Allocate offset maps at each level in the geometry
    call allocate_offsets(univ_list, counts, found)

    ! Calculate offsets for each target distribcell
    do i = 1, n_maps
      do j = 1, n_universes
        call calc_offsets(univ_list(i), i, universes(j), counts, found)
      end do
    end do

  end subroutine prepare_distribcell

!===============================================================================
! ALLOCATE_OFFSETS determines the number of maps needed and allocates required
! memory for distribcell offset tables
!===============================================================================

  recursive subroutine allocate_offsets(univ_list, counts, found)

    integer, intent(out), allocatable     :: univ_list(:) ! Target offsets
    integer, intent(out), allocatable     :: counts(:,:)  ! Target count
    logical, intent(out), allocatable     :: found(:,:)   ! Target found

    integer      :: i, j, k   ! Loop counters
    type(SetInt) :: cell_list ! distribells to track

    ! Begin gathering list of cells in distribcell tallies
    n_maps = 0

    ! List all cells referenced in distribcell filters.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % obj % filter)
        select type(filt => filters(tallies(i) % obj % filter(j)) % obj)
        type is (DistribcellFilter)
          call cell_list % add(filt % cell)
        end select
      end do
    end do

    ! List all cells with multiple (distributed) materials or temperatures.
    do i = 1, n_cells
      if (size(cells(i) % material) > 1 .or. size(cells(i) % sqrtkT) > 1) then
        call cell_list % add(i)
      end if
    end do

    ! Compute the number of unique universes containing these distribcells
    ! to determine the number of offset tables to allocate
    do i = 1, n_universes
      do j = 1, size(universes(i) % cells)
        if (cell_list % contains(universes(i) % cells(j))) then
          n_maps = n_maps + 1
        end if
      end do
    end do

    ! Allocate the list of offset tables for each unique universe
    allocate(univ_list(n_maps))

    ! Allocate list to accumulate target distribcell counts in each universe
    allocate(counts(n_universes, n_maps))
    counts(:,:) = 0

    ! Allocate list to track if target distribcells are found in each universe
    allocate(found(n_universes, n_maps))
    found(:,:) = .false.


    ! Search through universes for distributed cells and assign each one a
    ! unique distribcell array index.
    k = 1
    do i = 1, n_universes
      do j = 1, size(universes(i) % cells)
        if (cell_list % contains(universes(i) % cells(j))) then
          cells(universes(i) % cells(j)) % distribcell_index = k
          univ_list(k) = universes(i) % id
          k = k + 1
        end if
      end do
    end do

    ! Allocate the offset tables for lattices
    do i = 1, n_lattices
      associate(lat => lattices(i) % obj)
        select type(lat)

        type is (RectLattice)
          allocate(lat % offset(n_maps, lat % n_cells(1), lat % n_cells(2), &
                   lat % n_cells(3)))
        type is (HexLattice)
          allocate(lat % offset(n_maps, 2 * lat % n_rings - 1, &
               2 * lat % n_rings - 1, lat % n_axial))
        end select

        lat % offset(:, :, :, :) = 0
      end associate
    end do

    ! Allocate offset table for fill cells
    do i = 1, n_cells
      if (cells(i) % type /= FILL_MATERIAL) then
        allocate(cells(i) % offset(n_maps))
      end if
    end do

    ! Free up memory
    call cell_list % clear()

  end subroutine allocate_offsets

end module input_xml
