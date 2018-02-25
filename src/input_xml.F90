module input_xml

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: find
  use cmfd_input,       only: configure_cmfd
  use cmfd_header,      only: index_cmfd_mesh
  use constants
  use dict_header,      only: DictIntInt, DictCharInt, DictEntryCI
  use endf,             only: reaction_name
  use error,            only: fatal_error, warning, write_message, openmc_err_msg
  use geometry,         only: neighbor_lists
  use geometry_header
#ifdef CAD
  use cad_header
#endif
  use hdf5_interface
  use list_header,      only: ListChar, ListInt, ListReal
  use material_header
  use mesh_header
  use message_passing
  use mgxs_data,        only: create_macro_xs, read_mgxs
  use mgxs_interface
  use nuclide_header
  use output,           only: title, header, print_plot
  use photon_header
  use plot_header
  use random_lcg,       only: prn
  use surface_header
  use set_header,       only: SetChar
  use settings
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with, split_string, &
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
    subroutine adjust_indices() bind(C)
    end subroutine adjust_indices

    subroutine assign_temperatures() bind(C)
    end subroutine assign_temperatures

    subroutine count_cell_instances(univ_indx) bind(C)
      import C_INT32_T
      integer(C_INT32_T), intent(in), value :: univ_indx
    end subroutine count_cell_instances

    subroutine prepare_distribcell_c(cell_list, n) &
         bind(C, name="prepare_distribcell")
      import C_INT32_T, C_INT
      integer(C_INT),     intent(in), value :: n
      integer(C_INT32_T), intent(in)        :: cell_list(n)
    end subroutine prepare_distribcell_c

    subroutine read_surfaces(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine read_surfaces

    subroutine read_cells(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine read_cells

    subroutine read_lattices(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine read_lattices

    subroutine read_settings_xml() bind(C)
    end subroutine read_settings_xml

    subroutine read_materials(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine read_materials

    function find_root_universe() bind(C) result(root)
      import C_INT32_T
      integer(C_INT32_T) :: root
    end function find_root_universe

    function maximum_levels(univ) bind(C) result(n)
      import C_INT32_T, C_INT
      integer(C_INT32_T), intent(in), value :: univ
      integer(C_INT)                        :: n
    end function maximum_levels

    subroutine set_particle_energy_bounds(particle, E_min, E_max) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: particle
      real(C_DOUBLE), value :: E_min
      real(C_DOUBLE), value :: E_max
    end subroutine
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
    call count_cell_instances(root_universe)

    ! After reading input and basic geometry setup is complete, build lists of
    ! neighboring cells for efficient tracking
    call neighbor_lists()

    ! Assign temperatures to cells that don't have temperatures already assigned
    call assign_temperatures()

    ! Determine desired temperatures for each nuclide and S(a,b) table
    call get_temperatures(nuc_temps, sab_temps)

    ! Check to make sure there are not too many nested coordinate levels in the
    ! geometry since the coordinate list is statically allocated for performance
    ! reasons
    if (maximum_levels(root_universe) > MAX_COORD) then
      call fatal_error("Too many nested coordinate levels in the geometry. &
           &Try increasing the maximum number of coordinate levels by &
           &providing the CMake -Dmaxcoord= option.")
    end if

  end subroutine finalize_geometry

!===============================================================================
! READ_SETTINGS_XML reads data from a settings.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_settings_xml_f(root_ptr) bind(C)
    type(C_PTR), value :: root_ptr

    integer :: i
    integer :: n
    integer, allocatable :: temp_int_array(:)
    integer :: n_tracks
    type(XMLNode) :: root
    type(XMLNode) :: node_sp
    type(XMLNode) :: node_res_scat
    type(XMLNode) :: node_vol
    type(XMLNode), allocatable :: node_vol_list(:)

    ! Get proper XMLNode type given pointer
    root % ptr = root_ptr

    ! Check for use of CAD geometry
    if (check_for_node(root, "dagmc")) then
#ifdef CAD
       call get_node_value_bool(root, "dagmc", dagmc)
#else
       if (dagmc) then
          call fatal_error("CAD mode unsupported for this build of OpenMC")
       end if
#endif
    end if

    if (run_mode == MODE_EIGENVALUE) then
      ! Preallocate space for keff and entropy by generation
      call k_generation % reserve(n_max_batches*gen_per_batch)
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
    else
      ! If no <source_point> tag was present, by default we keep source bank in
      ! statepoint file and write it out at statepoints intervals
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

    ! Resonance scattering parameters
    if (check_for_node(root, "resonance_scattering")) then
      node_res_scat = root % child("resonance_scattering")

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

  end subroutine read_settings_xml_f

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml()

    integer :: i
    integer :: n, n_rlats, n_hlats
    integer :: univ_id
    integer :: n_cells_in_univ
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    type(Cell),     pointer :: c
    class(Lattice), pointer :: lat
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_cell
    type(XMLNode) :: node_lat
    type(XMLNode), allocatable :: node_cell_list(:)
    type(XMLNode), allocatable :: node_rlat_list(:)
    type(XMLNode), allocatable :: node_hlat_list(:)
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains
#ifdef CAD
    if (dagmc) then
       call write_message("Reading CAD geometry...", 5)
       call load_cad_geometry()
       call allocate_surfaces()
       call allocate_cells()
       return
    end if
#endif
    
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

    call allocate_surfaces()
    ! Allocate surfaces array
    do i = 1, n_surfaces
      surfaces(i) % ptr = surface_pointer(i - 1);

      if (surfaces(i) % bc() /= BC_TRANSMIT) boundary_exists = .true.
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

    call read_cells(root % ptr)

    call allocate_cells()

    ! Get pointer to list of XML <cell>
    call get_node_list(root, "cell", node_cell_list)

    ! Get number of <cell> tags
    n_cells = size(node_cell_list)

    ! Check for no cells
    if (n_cells == 0) then
      call fatal_error("No cells found in geometry.xml!")
    end if


    n_universes = 0
    do i = 1, n_cells
      c => cells(i)

      c % ptr = cell_pointer(i - 1)

      ! Get pointer to i-th cell node
      node_cell = node_cell_list(i)

      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has(c % id())) then
        call fatal_error("Two or more cells use the same unique ID: " &
             // to_str(c % id()))
      end if

      ! Rotation matrix
      if (check_for_node(node_cell, "rotation")) then
        ! Rotations can only be applied to cells that are being filled with
        ! another universe
        if (c % fill() == C_NONE) then
          call fatal_error("Cannot apply a rotation to cell " // trim(to_str(&
               &c % id())) // " because it is not filled with another universe")
        end if

        ! Read number of rotation parameters
        n = node_word_count(node_cell, "rotation")
        if (n /= 3) then
          call fatal_error("Incorrect number of rotation parameters on cell " &
               // to_str(c % id()))
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

      ! Add cell to dictionary
      call cell_dict % set(c % id(), i)

      ! For cells, we also need to check if there's a new universe --
      ! also for every cell add 1 to the count of cells for the
      ! specified universe
      univ_id = c % universe()
      if (.not. cells_in_univ_dict % has(univ_id)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call universe_dict % set(univ_id, n_universes - 1)
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get(univ_id)
      end if
      call cells_in_univ_dict % set(univ_id, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    call read_lattices(root % ptr)

    ! Get pointer to list of XML <lattice>
    call get_node_list(root, "lattice", node_rlat_list)
    call get_node_list(root, "hex_lattice", node_hlat_list)

    ! Allocate lattices array
    n_rlats = size(node_rlat_list)
    n_hlats = size(node_hlat_list)
    allocate(lattices(n_rlats + n_hlats))

    RECT_LATTICES: do i = 1, n_rlats
      lat => lattices(i)
      lat % ptr = lattice_pointer(i - 1)

      ! Get pointer to i-th lattice
      node_lat = node_rlat_list(i)

      ! Add lattice to dictionary
      call lattice_dict % set(lat % id(), i)

    end do RECT_LATTICES

    HEX_LATTICES: do i = 1, n_hlats
      lat => lattices(n_rlats + i)
      lat % ptr = lattice_pointer(n_rlats + i - 1)

      ! Get pointer to i-th lattice
      node_lat = node_hlat_list(i)

      ! Add lattice to dictionary
      call lattice_dict % set(lat % id(), n_rlats + i)

    end do HEX_LATTICES

    ! ==========================================================================
    ! SETUP UNIVERSES

    ! Allocate universes, universe cell arrays, and assign base universe
    root_universe = find_root_universe()

    ! Clear dictionary
    call cells_in_univ_dict%clear()

    ! Close geometry XML file
    call doc % clear()

  end subroutine read_geometry_xml

  subroutine allocate_surfaces()
    integer :: i

    ! Allocate surfaces array
    allocate(surfaces(n_surfaces))

    do i = 1, n_surfaces
      surfaces(i) % ptr = surface_pointer_c(i - 1);
      ! Add surface to dictionary
      call surface_dict % set(surfaces(i) % id(), i)
    end do

    end subroutine allocate_surfaces

  subroutine allocate_cells()
    integer :: i
    type(Cell), pointer :: c
    
    ! Allocate cells array
    allocate(cells(n_cells))
    
    do i = 1, n_cells
       c => cells(i)
       
       c % ptr = cell_pointer_c(i - 1)

       ! Check to make sure 'id' hasn't been used
       if (cell_dict % has(c % id())) then
          call fatal_error("Two or more cells use the same unique ID: " &
               // to_str(c % id()))
       end if
       
       ! Add cell to dictionary
       call cell_dict % set(c % id(), i)

    end do
  end subroutine allocate_cells

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
                 & user's guide at http://openmc.readthedocs.io for &
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
               &guide at http://openmc.readthedocs.io for information on &
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
    integer :: index_element  ! index in elements
    integer :: index_sab      ! index in sab_tables
    logical :: file_exists    ! does materials.xml exist?
    character(20)           :: name         ! name of nuclide, e.g. U235
    character(3)            :: element      ! name of element, e.g. Zr
    character(MAX_WORD_LEN) :: units        ! units on density
    character(MAX_LINE_LEN) :: filename     ! absolute path to materials.xml
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

    call read_materials(root % ptr)

    ! Get pointer to list of XML <material>
    call get_node_list(root, "material", node_mat_list)

    ! Allocate materials array
    n_materials = size(node_mat_list)
    allocate(materials(n_materials))
    allocate(material_temps(n_materials))

    ! Initialize count for number of nuclides/S(a,b) tables
    index_nuclide = 0
    index_element = 0
    index_sab = 0

    do i = 1, n_materials
      mat => materials(i)

      mat % ptr = material_pointer(i - 1)

      ! Get pointer to i-th material node
      node_mat = node_mat_list(i)

      ! Check if material is depletable
      if (check_for_node(node_mat, "depletable")) then
        call get_node_value(node_mat, "depletable", mat % depletable)
      end if

      ! Copy material name
      if (check_for_node(node_mat, "name")) then
        call get_node_value(node_mat, "name", mat % name)
      end if

      ! Get material default temperature
      if (check_for_node(node_mat, "temperature")) then
        call get_node_value(node_mat, "temperature", material_temps(i))
      else
        material_temps(i) = -1.0
      end if

      ! Get pointer to density element
      if (check_for_node(node_mat, "density")) then
        node_dens = node_mat % child("density")
      else
        call fatal_error("Must specify density element in material " &
             // trim(to_str(mat % id())))
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
               // trim(to_str(mat % id())) // ".")
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
               // "' specified on material " // trim(to_str(mat % id())))
        end select
      end if

      ! Issue error if elements are provided
      call get_node_list(node_mat, "element", node_ele_list)

      if (size(node_ele_list) > 0) then
        call fatal_error("Unable to add an element to material " &
             // trim(to_str(mat % id())) // " since the element option has &
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
             &material " // trim(to_str(mat % id())))
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
             // trim(to_str(mat % id())))
      else if (size(node_macro_list) == 1) then

        node_nuc = node_macro_list(1)

        ! Check for empty name on nuclide
        if (.not. check_for_node(node_nuc, "name")) then
          call fatal_error("No name specified on macroscopic data in material " &
               // trim(to_str(mat % id())))
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
                 // trim(to_str(mat % id())))
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
      allocate(mat % element(n))
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

        ! If the corresponding element hasn't been encountered yet and photon
        ! transport will be used, we need to add its symbol to the element_dict
        if (photon_transport) then
          element = name(1:scan(name, '0123456789') - 1)

          ! Make sure photon cross section data is available
          if (.not. library_dict % has(to_lower(element))) then
            call fatal_error("Could not find element " // trim(element) &
                 // " in cross_sections data file!")
          end if

          if (.not. element_dict % has(element)) then
            index_element = index_element + 1
            mat % element(j) = index_element

            call element_dict % set(element, index_element)
          else
            mat % element(j) = element_dict % get(element)
          end if
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
             // to_str(mat % id()))
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
      call material_dict % set(mat % id(), i)
    end do

    ! Set total number of nuclides and S(a,b) tables
    n_nuclides = index_nuclide
    n_elements = index_element
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
    integer :: l             ! loop over bins
    integer :: filter_id     ! user-specified identifier for filter
    integer :: tally_id      ! user-specified identifier for filter
    integer :: i_filt        ! index in filters array
    integer :: i_elem        ! index of entry in dictionary
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filter      ! number of filters
    integer :: n_scores      ! number of scores
    integer :: n_user_trig   ! number of user-specified tally triggers
    integer :: trig_ind      ! index of triggers array for each tally
    integer :: user_trig_ind ! index of user-specified triggers for each tally
    integer :: i_start, i_end
    integer(C_INT) :: err
    real(8) :: threshold     ! trigger convergence threshold
    integer :: MT            ! user-specified MT for score
    logical :: file_exists   ! does tallies.xml file exist?
    integer, allocatable :: temp_filter(:) ! temporary filter indices
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: word
    character(MAX_WORD_LEN) :: score_name
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: sarray(:)
    type(DictCharInt) :: trigger_scores
    type(TallyFilterContainer), pointer :: f
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_tal
    type(XMLNode) :: node_filt
    type(XMLNode) :: node_trigger
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
    call read_meshes(root % ptr)

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

      ! Make sure bins have been set
      select case(temp_str)
      case ("energy", "energyout", "mu", "polar", "azimuthal")
        if (.not. check_for_node(node_filt, "bins")) then
          call fatal_error("Bins not set in filter " // trim(to_str(filter_id)))
        end if
      case ("mesh", "meshsurface", "universe", "material", "cell", "distribcell", &
            "cellborn", "cellfrom", "surface", "delayedgroup")
        if (.not. check_for_node(node_filt, "bins")) then
          call fatal_error("Bins not set in filter " // trim(to_str(filter_id)))
        end if
      end select

      ! Allocate according to the filter type
      err = openmc_filter_set_type(i_start + i - 1, to_c_string(temp_str))
      if (err /= 0) call fatal_error(to_f_string(openmc_err_msg))

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
      err = openmc_tally_allocate(i_start + i - 1, &
           C_CHAR_'generic' // C_NULL_CHAR)

      ! Get pointer to tally
      associate (t => tallies(i_start + i - 1) % obj)

      ! Get pointer to tally xml node
      node_tal = node_tal_list(i)

      ! Copy and set tally id
      if (check_for_node(node_tal, "id")) then
        call get_node_value(node_tal, "id", tally_id)
        err = openmc_tally_set_id(i_start + i - 1, tally_id)
        if (err /= 0) call fatal_error(to_f_string(openmc_err_msg))
      else
        call fatal_error("Must specify id for tally in tally XML file.")
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

        ! Append the score to the list of possible trigger scores
        do j = 1, n_words
          sarray(j) = to_lower(sarray(j))
          score_name = trim(sarray(j))

          if (trigger_on) call trigger_scores % set(trim(score_name), j)

        end do
        n_scores = n_words

        ! Allocate score storage accordingly
        allocate(t % score_bins(n_scores))

        ! Check the validity of the scores and their filters
        do j = 1, n_scores
          score_name = sarray(j)

          ! Check if delayed group filter is used with any score besides
          ! delayed-nu-fission or decay-rate
          if ((score_name /= 'delayed-nu-fission' .and. &
               score_name /= 'decay-rate') .and. &
               t % find_filter(FILTER_DELAYEDGROUP) > 0) then
            call fatal_error("Cannot tally " // trim(score_name) // " with a &
                 &delayedgroup filter.")
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

          case ('total', '(n,total)')
            t % score_bins(j) = SCORE_TOTAL
            if (t % find_filter(FILTER_ENERGYOUT) > 0) then
              call fatal_error("Cannot tally total reaction rate with an &
                   &outgoing energy filter.")
            end if

          case ('scatter')
            t % score_bins(j) = SCORE_SCATTER
            if (t % find_filter(FILTER_ENERGYOUT) > 0 .or. &
                 t % find_filter(FILTER_LEGENDRE) > 0) then
              ! Set tally estimator to analog
              t % estimator = ESTIMATOR_ANALOG
            end if

          case ('nu-scatter')
            t % score_bins(j) = SCORE_NU_SCATTER

            ! Set tally estimator to analog for CE mode
            ! (MG mode has all data available without a collision being
            ! necessary)
            if (run_CE) then
              t % estimator = ESTIMATOR_ANALOG
            else
              if (t % find_filter(FILTER_ENERGYOUT) > 0 .or. &
                   t % find_filter(FILTER_LEGENDRE) > 0) then
                ! Set tally estimator to analog
                t % estimator = ESTIMATOR_ANALOG
              end if
            end if

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

              ! Check to make sure that mesh surface currents are not desired as well
              if (t % find_filter(FILTER_MESHSURFACE) > 0) then
                call fatal_error("Cannot tally mesh surface currents &
                     &in the same tally as normal surface currents")
              end if

              t % type = TALLY_SURFACE
              t % score_bins(j) = SCORE_CURRENT

            else if (t % find_filter(FILTER_MESHSURFACE) > 0) then
              t % score_bins(j) = SCORE_CURRENT
              t % type = TALLY_MESH_SURFACE

              ! Check to make sure that current is the only desired response
              ! for this tally
              if (n_words > 1) then
                call fatal_error("Cannot tally other scores in the &
                     &same tally as surface currents")
              end if
            else
                call fatal_error("Cannot tally currents without surface &
                     &type filters")
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
            ! First look for deprecated scores
            if (starts_with(trim(score_name), 'scatter-') .or. &
                 starts_with(trim(score_name), 'nu-scatter-') .or. &
                 starts_with(trim(score_name), 'total-y') .or. &
                 starts_with(trim(score_name), 'flux-y')) then
              call fatal_error(trim(score_name) // " is no longer available.")
            end if

            ! Assume that user has specified an MT number
            MT = int(str_to_int(score_name))

            if (MT /= ERROR_INT) then
              ! Specified score was an integer
              if (MT > 1) then
                t % score_bins(j) = MT
              else
                call fatal_error("Invalid MT on <scores>: " // trim(score_name))
              end if

            else
              ! Specified score was not an integer
              call fatal_error("Unknown scoring function: " // trim(score_name))
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

        ! Deallocate temporary string array of scores
        deallocate(sarray)

        ! Check that no duplicate scores exist
        do j = 1, n_scores - 1
          do k = j + 1, n_scores
            if (t % score_bins(j) == t % score_bins(k)) then
              call fatal_error("Duplicate score of type '" // trim(&
                   reaction_name(t % score_bins(j))) // "' found in tally " &
                   // trim(to_str(t % id)))
            end if
          end do
        end do

        ! Check if tally is compatible with particle type
        if (photon_transport) then
          if (t % find_filter(FILTER_PARTICLE) == 0) then
            do j = 1, n_scores
              select case (t % score_bins(j))
              case (SCORE_INVERSE_VELOCITY)
                call fatal_error("Particle filter must be used with photon &
                     &transport on and inverse velocity score")
              case (SCORE_FLUX, SCORE_TOTAL, SCORE_SCATTER, SCORE_NU_SCATTER, &
                   SCORE_ABSORPTION, SCORE_FISSION, SCORE_NU_FISSION, &
                   SCORE_CURRENT, SCORE_EVENTS, SCORE_DELAYED_NU_FISSION, &
                   SCORE_PROMPT_NU_FISSION, SCORE_DECAY_RATE)
                call warning("Particle filter is not used with photon transport&
                     & on and " // trim(to_str(t % score_bins(j))) // " score")
              end select
            end do
          else
            select type(filt => filters(t % find_filter(FILTER_PARTICLE)) % obj)
            type is (ParticleFilter)
              do l = 1, filt % n_bins
                if (filt % particles(l) == ELECTRON .or. filt % particles(l) == POSITRON) then
                  t % estimator = ESTIMATOR_ANALOG
                end if
              end do
            end select
          end if
        else
          if (t % find_filter(FILTER_PARTICLE) > 0) then
            select type(filt => filters(t % find_filter(FILTER_PARTICLE)) % obj)
            type is (ParticleFilter)
              do l = 1, filt % n_bins
                if (filt % particles(l) /= NEUTRON) then
                  call warning("Particle filter other than NEUTRON used with &
                       &photon transport turned off. All tallies for particle &
                       &type " // trim(to_str(filt % particles(l))) // " will have no scores")
                end if
              end do
            end select
          end if
        end if
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
    integer(C_INT) :: err, idx
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

              pl % index_meshlines_mesh = index_ufs_mesh

            case ('cmfd')

              if (.not. cmfd_run) then
                call fatal_error("Need CMFD run to plot CMFD mesh for &
                     &meshlines on plot " // trim(to_str(pl % id)))
              end if

              pl % index_meshlines_mesh = index_cmfd_mesh

            case ('entropy')

              if (index_entropy_mesh < 0) then
                call fatal_error("No entropy mesh for meshlines on plot " &
                     // trim(to_str(pl % id)))
              end if

              pl % index_meshlines_mesh = index_entropy_mesh

            case ('tally')

              ! Ensure that there is a mesh id if the type is tally
              if (check_for_node(node_meshlines, "id")) then
                call get_node_value(node_meshlines, "id", meshid)
              else
                call fatal_error("Must specify a mesh id for meshlines tally &
                     &mesh specification in plot " // trim(to_str(pl % id)))
              end if

              ! Check if the specified tally mesh exists
              err = openmc_get_mesh_index(meshid, idx)
              if (err /= 0) then
                call fatal_error("Could not find mesh " &
                     // trim(to_str(meshid)) // " specified in meshlines for &
                     &plot " // trim(to_str(pl % id)))
              end if
              pl % index_meshlines_mesh = idx

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
        case ('photon')
          libraries(i) % type = LIBRARY_PHOTON
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
    energy_min(NEUTRON) = energy_bins(num_energy_groups + 1)
    energy_max(NEUTRON) = energy_bins(1)
    call set_particle_energy_bounds(NEUTRON, energy_min(NEUTRON), &
         energy_max(NEUTRON))

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
            awr = get_awr_c(mat % nuclide(j))
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
              awr = get_awr_c(mat % nuclide(j))
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
    integer :: i_element
    integer :: i_sab
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    logical :: mp_found     ! if windowed multipole libraries were found
    character(MAX_WORD_LEN) :: name
    character(3) :: element
    type(SetChar) :: already_read
    type(SetChar) :: element_already_read

    allocate(nuclides(n_nuclides))
    allocate(elements(n_elements))
    allocate(sab_tables(n_sab_tables))
    if (photon_transport .and. electron_treatment == ELECTRON_TTB) then
      allocate(ttb(n_materials))
    end if

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
            energy_min(NEUTRON) = max(energy_min(NEUTRON), &
                 nuclides(i_nuclide) % grid(1) % energy(1))
            energy_max(NEUTRON) = min(energy_max(NEUTRON), nuclides(i_nuclide) % &
                 grid(1) % energy(size(nuclides(i_nuclide) % grid(1) % energy)))
            call set_particle_energy_bounds(NEUTRON, energy_min(NEUTRON), &
                 energy_max(NEUTRON))
          end if

          ! Add name and alias to dictionary
          call already_read % add(name)

          ! Check if elemental data has been read, if needed
          element = name(1:scan(name, '0123456789') - 1)
          if (photon_transport) then
            if (.not. element_already_read % contains(element)) then
              ! Read photon interaction data from HDF5 photon library
              i_library = library_dict % get(to_lower(element))
              i_element = element_dict % get(element)
              call write_message('Reading ' // trim(element) // ' from ' // &
                   trim(libraries(i_library) % path), 6)

              ! Open file and make sure version is sufficient
              file_id = file_open(libraries(i_library) % path, 'r')
              call check_data_version(file_id)

              ! Read element data from HDF5
              group_id = open_group(file_id, element)
              call elements(i_element) % from_hdf5(group_id)
              call close_group(group_id)
              call file_close(file_id)

              ! Determine if minimum/maximum energy for this element is
              ! greater/less than the previous
              if (size(elements(i_element) % energy) >= 1) then
                energy_min(PHOTON) = max(energy_min(PHOTON), &
                     exp(elements(i_element) % energy(1)))
                energy_max(PHOTON) = min(energy_max(PHOTON), &
                     exp(elements(i_element) % energy(size(elements(i_element) &
                     % energy))))
                call set_particle_energy_bounds(PHOTON, energy_min(PHOTON), &
                     energy_max(PHOTON))
              end if

              ! Add element to set
              call element_already_read % add(element)
            end if
          end if

          ! Read multipole file into the appropriate entry on the nuclides array
          if (temperature_multipole) call read_multipole_data(i_nuclide)
        end if

        ! Check if material is fissionable
        if (nuclides(materials(i) % nuclide(j)) % fissionable) then
          materials(i) % fissionable = .true.
        end if
      end do

      ! Generate material bremsstrahlung data for electrons and positrons
      if (photon_transport .and. electron_treatment == ELECTRON_TTB) then
        call bremsstrahlung_init(ttb(i) % electron, i, ELECTRON)
        call bremsstrahlung_init(ttb(i) % positron, i, POSITRON)
      end if
    end do

    if (photon_transport .and. electron_treatment == ELECTRON_TTB) then
      ! Deallocate element bremsstrahlung DCS and stopping power data since
      ! only the material bremsstrahlung data is needed
      do i = 1, size(elements)
        if (allocated(elements(i) % stopping_power_collision)) &
             deallocate(elements(i) % stopping_power_collision)
        if (allocated(elements(i) % stopping_power_radiative)) &
             deallocate(elements(i) % stopping_power_radiative)
        if (allocated(elements(i) % dcs)) deallocate(elements(i) % dcs)
        if (allocated(ttb_k_grid)) deallocate(ttb_k_grid)
      end do

      ! Determine if minimum/maximum energy for bremsstrahlung is greater/less
      ! than the current minimum/maximum
      if (size(ttb_e_grid) >= 1) then
        energy_min(PHOTON) = max(energy_min(PHOTON), ttb_e_grid(1))
        energy_max(PHOTON) = min(energy_max(PHOTON), ttb_e_grid(size(ttb_e_grid)))
        call set_particle_energy_bounds(PHOTON, energy_min(PHOTON), &
             energy_max(PHOTON))
      end if

      ! Take logarithm of energies since they are log-log interpolated
      ttb_e_grid = log(ttb_e_grid)
    end if

    ! Set up logarithmic grid for nuclides
    do i = 1, size(nuclides)
      call nuclides(i) % init_grid(energy_min(NEUTRON), &
           energy_max(NEUTRON), n_log_bins)
    end do
    log_spacing = log(energy_max(NEUTRON)/energy_min(NEUTRON)) / n_log_bins

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
      ! If a nuclide is present in a material that's not used in the model, its
      ! grid has not been allocated
      if (size(nuclides(i) % grid) > 0) then
        if (nuclides(i) % grid(1) % energy(size(nuclides(i) % grid(1) % energy)) &
             == energy_max(NEUTRON)) then
          call write_message("Maximum neutron transport energy: " // &
               trim(to_str(energy_max(NEUTRON))) // " eV for " // &
               trim(adjustl(nuclides(i) % name)), 7)
          exit
        end if
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
      call nuc % multipole % from_hdf5(filename)
      nuc % mp_present = .true.

    end associate

  end subroutine read_multipole_data

!===============================================================================
! PREPARE_DISTRIBCELL initializes any distribcell filters present and sets the
! offsets for distribcells
!===============================================================================

  subroutine prepare_distribcell()

    integer :: i, j
    type(SetInt)  :: cell_list  ! distribcells to track
    integer(C_INT32_T), allocatable :: cell_list_c(:)

    ! Find all cells listed in a distribcell filter.
    do i = 1, n_tallies
      do j = 1, size(tallies(i) % obj % filter)
        select type(filt => filters(tallies(i) % obj % filter(j)) % obj)
        type is (DistribcellFilter)
          call cell_list % add(filt % cell)
        end select
      end do
    end do

    allocate(cell_list_c(cell_list % size()))
    do i = 1, cell_list % size()
      cell_list_c(i) = cell_list % get_item(i) - 1
    end do
    call prepare_distribcell_c(cell_list_c, cell_list % size())

  end subroutine prepare_distribcell

end module input_xml
