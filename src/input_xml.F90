module input_xml

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: find
  use constants
  use dict_header,      only: DictIntInt, DictCharInt, DictEntryCI
  use endf,             only: reaction_name
  use error,            only: fatal_error, warning, write_message, openmc_err_msg
  use geometry_header
#ifdef DAGMC
  use dagmc_header
#endif
  use hdf5_interface
  use material_header
  use mesh_header
  use message_passing
  use mgxs_interface
  use nuclide_header
  use photon_header
  use random_lcg,       only: prn
  use surface_header
  use settings
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, &
                              starts_with, ends_with, to_c_string
  use tally
  use tally_header,     only: openmc_extend_tallies
  use tally_derivative_header
  use tally_filter_header
  use tally_filter
  use volume_header
  use xml_interface

  implicit none
  save

  interface
    subroutine count_cell_instances(univ_indx) bind(C)
      import C_INT32_T
      integer(C_INT32_T), intent(in), value :: univ_indx
    end subroutine count_cell_instances

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

    subroutine read_plots(node_ptr) bind(C)
      import C_PTR
      type(C_PTR) :: node_ptr
    end subroutine read_plots

    subroutine set_particle_energy_bounds(particle, E_min, E_max) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: particle
      real(C_DOUBLE), value :: E_min
      real(C_DOUBLE), value :: E_max
    end subroutine
  end interface

contains

!===============================================================================
! READ_SETTINGS_XML reads data from a settings.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_settings_xml_f(root_ptr) bind(C)
    type(C_PTR), value :: root_ptr

    integer :: i
    integer :: n
    type(XMLNode) :: root
    type(XMLNode) :: node_vol
    type(XMLNode), allocatable :: node_vol_list(:)

    ! Get proper XMLNode type given pointer
    root % ptr = root_ptr

    call get_node_list(root, "volume_calc", node_vol_list)
    n = size(node_vol_list)
    allocate(volume_calcs(n))
    do i = 1, n
      node_vol = node_vol_list(i)
      call volume_calcs(i) % from_xml(node_vol)
    end do

  end subroutine read_settings_xml_f


#ifdef DAGMC

!===============================================================================
! READ_GEOMETRY_DAGMC reads data from a DAGMC .h5m file, checking
! for material properties and surface boundary conditions
! some universe information is spoofed for now
!===============================================================================

  subroutine read_geometry_dagmc()

    integer :: i
    integer :: univ_id
    integer :: n_cells_in_univ
    logical :: file_exists
    character(MAX_LINE_LEN) :: filename
    type(Cell),     pointer :: c
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains

    ! Check if dagmc.h5m exists
    filename = trim(path_input) // "dagmc.h5m"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      call fatal_error("Geometry DAGMC file '" // trim(filename) // "' does not &
           &exist!")
    end if

    call write_message("Reading DAGMC geometry...", 5)
    call load_dagmc_geometry()
    call allocate_surfaces()
    call allocate_cells()

    ! setup universe data structs
    do i = 1, n_cells
      c => cells(i)
      ! additional metadata spoofing
      univ_id = c % universe()

      if (.not. cells_in_univ_dict % has(univ_id)) then
        n_universes = n_universes + 1
        n_cells_in_univ = 1
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get(univ_id)
      end if
      call cells_in_univ_dict % set(univ_id, n_cells_in_univ)
    end do

    root_universe = find_root_universe()

  end subroutine read_geometry_dagmc

#endif

!===============================================================================
! READ_GEOMETRY_XML reads data from a geometry.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_geometry_xml() bind(C)

    integer :: i, n
    integer :: univ_id
    integer :: n_cells_in_univ
    real(8) :: phi, theta, psi
    logical :: file_exists
    logical :: boundary_exists
    character(MAX_LINE_LEN) :: filename
    type(Cell),     pointer :: c
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_cell
    type(XMLNode), allocatable :: node_cell_list(:)
    type(VectorInt) :: univ_ids      ! List of all universe IDs
    type(DictIntInt) :: cells_in_univ_dict ! Used to count how many cells each
                                           ! universe contains
#ifdef DAGMC
    if (dagmc) then
      call read_geometry_dagmc()
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

    ! Allocate surfaces array
    allocate(surfaces(surfaces_size()))
    do i = 1, size(surfaces)
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
        call univ_ids % push_back(univ_id)
      else
        n_cells_in_univ = 1 + cells_in_univ_dict % get(univ_id)
      end if
      call cells_in_univ_dict % set(univ_id, n_cells_in_univ)

    end do

    ! ==========================================================================
    ! READ LATTICES FROM GEOMETRY.XML

    call read_lattices(root % ptr)

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
    allocate(surfaces(surfaces_size()))

    do i = 1, size(surfaces)
      surfaces(i) % ptr = surface_pointer(i - 1);
    end do

  end subroutine allocate_surfaces

  subroutine allocate_cells()
    integer :: i
    type(Cell), pointer :: c

    ! Allocate cells array
    allocate(cells(n_cells))

    do i = 1, n_cells
      c => cells(i)
      c % ptr = cell_pointer(i - 1)
      ! Check to make sure 'id' hasn't been used
      if (cell_dict % has(c % id())) then
        call fatal_error("Two or more cells use the same unique ID: " &
               // to_str(c % id()))
      end if
      ! Add cell to dictionary
      call cell_dict % set(c % id(), i)
    end do
  end subroutine allocate_cells

  subroutine read_materials_xml() bind(C)
    logical :: file_exists    ! does materials.xml exist?
    character(MAX_LINE_LEN) :: filename     ! absolute path to materials.xml
    type(XMLDocument) :: doc
    type(XMLNode) :: root

    interface
      function nuclides_size() bind(C) result(n)
        import C_INT
        integer(C_INT) :: n
      end function

      function elements_size() bind(C) result(n)
        import C_INT
        integer(C_INT) :: n
      end function
    end interface

    ! Display output message
    call write_message("Reading materials XML file...", 5)

    doc % ptr = C_NULL_PTR

#ifdef DAGMC
    if (dagmc) then
      doc % ptr = read_uwuw_materials()
    end if
#endif

    if (.not. c_associated(doc % ptr)) then
      ! Check if materials.xml exists
      filename = trim(path_input) // "materials.xml"
      inquire(FILE=filename, EXIST=file_exists)
      if (.not. file_exists) then
        call fatal_error("Material XML file '" // trim(filename) // "' does not &
             &exist!")
    end if

    ! Parse materials.xml file
    call doc % load_file(filename)

    end if

    root = doc % document_element()
    call read_materials(root % ptr)

    ! Set total number of nuclides and elements
    n_nuclides = nuclides_size()
    n_elements = elements_size()
    allocate(nuclides(n_nuclides))

    ! Close materials XML file
    call doc % clear()

  end subroutine read_materials_xml

!===============================================================================
! READ_TALLIES_XML reads data from a tallies.xml file and parses it, checking
! for errors and placing properly-formatted data in the right data structures
!===============================================================================

  subroutine read_tallies_xml() bind(C)

    integer :: i             ! loop over user-specified tallies
    integer :: j             ! loop over words
    integer :: l             ! loop over bins
    integer :: filter_id     ! user-specified identifier for filter
    integer :: tally_id      ! user-specified identifier for filter
    integer :: deriv_id
    integer :: i_filt        ! index in filters array
    integer :: n             ! size of arrays in mesh specification
    integer :: n_words       ! number of words read
    integer :: n_filter      ! number of filters
    integer :: i_start, i_end
    integer(C_INT) :: err
    logical :: file_exists   ! does tallies.xml file exist?
    integer, allocatable :: temp_filter(:) ! temporary filter indices
    logical :: has_energyout
    integer :: particle_filter_index
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: temp_str
    type(TallyFilterContainer), pointer :: f
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_tal
    type(XMLNode) :: node_filt
    type(XMLNode), allocatable :: node_tal_list(:)
    type(XMLNode), allocatable :: node_filt_list(:)
    type(TallyDerivative), pointer :: deriv

    interface
      subroutine tally_init_from_xml(tally_ptr, xml_node) bind(C)
        import C_PTR
        type(C_PTR), value :: tally_ptr
        type(C_PTR) :: xml_node
      end subroutine

      subroutine tally_set_scores(tally_ptr, xml_node) bind(C)
        import C_PTR
        type(C_PTR), value :: tally_ptr
        type(C_PTR) :: xml_node
      end subroutine

      subroutine tally_set_nuclides(tally_ptr, xml_node) bind(C)
        import C_PTR
        type(C_PTR), value :: tally_ptr
        type(C_PTR) :: xml_node
      end subroutine

      subroutine tally_init_triggers(tally_ptr, i_tally, xml_node) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: tally_ptr
        integer(C_INT), value :: i_tally
        type(C_PTR) :: xml_node
      end subroutine

      subroutine read_tally_derivatives(node_ptr) bind(C)
        import C_PTR
        type(C_PTR) :: node_ptr
      end subroutine
    end interface

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

    call read_tally_derivatives(root % ptr)

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

      call tally_init_from_xml(t % ptr, node_tal % ptr)

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
                 // trim(to_str(t % id())))
          end if

          ! Store the index of the filter
          temp_filter(j) = i_filt - 1
        end do

        ! Set the filters
        err = openmc_tally_set_filters(i_start + i - 1, n_filter, temp_filter)
      end if
      deallocate(temp_filter)

      ! Check for the presence of certain filter types
      has_energyout = (t % energyout_filter() > 0)
      particle_filter_index = 0
      do j = 1, t % n_filters()
        select type (filt => filters(t % filter(j) + 1) % obj)
        type is (ParticleFilter)
          particle_filter_index = j
        end select
      end do

      ! Change the tally estimator if a filter demands it
      do j = 1, t % n_filters()
        select type (filt => filters(t % filter(j) + 1) % obj)
        type is (EnergyoutFilter)
          call t % set_estimator(ESTIMATOR_ANALOG)
        type is (LegendreFilter)
          call t % set_estimator(ESTIMATOR_ANALOG)
        type is (SphericalHarmonicsFilter)
          if (filt % cosine() == COSINE_SCATTER) then
            call t % set_estimator(ESTIMATOR_ANALOG)
          end if
        type is (SpatialLegendreFilter)
          call t % set_estimator(ESTIMATOR_COLLISION)
        type is (ZernikeFilter)
          call t % set_estimator(ESTIMATOR_COLLISION)
        type is (ZernikeRadialFilter)
          call t % set_estimator(ESTIMATOR_COLLISION)
        end select
      end do

      ! =======================================================================
      ! READ DATA FOR NUCLIDES

      call tally_set_nuclides(t % ptr, node_tal % ptr)

      ! =======================================================================
      ! READ DATA FOR SCORES

      call tally_set_scores(t % ptr, node_tal % ptr)

      if (check_for_node(node_tal, "scores")) then
        n_words = node_word_count(node_tal, "scores")

        ! Check if tally is compatible with particle type
        if (photon_transport) then
          if (particle_filter_index == 0) then
            do j = 1, t % n_score_bins()
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
            select type(filt => filters(particle_filter_index) % obj)
            type is (ParticleFilter)
              do l = 1, filt % n_bins
                if (filt % particles(l) == ELECTRON .or. filt % particles(l) == POSITRON) then
                  call t % set_estimator(ESTIMATOR_ANALOG)
                end if
              end do
            end select
          end if
        else
          if (particle_filter_index > 0) then
            select type(filt => filters(particle_filter_index) % obj)
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
             // trim(to_str(t % id())) // ".")
      end if

      ! Check for a tally derivative.
      if (check_for_node(node_tal, "derivative")) then
        call get_node_value(node_tal, "derivative", deriv_id)

        ! Find the derivative with the given id, and store it's index.
        do j = 0, n_tally_derivs() - 1
          deriv => tally_deriv_c(j)
          if (deriv % id == deriv_id) then
            call t % set_deriv(j)
            ! Only analog or collision estimators are supported for differential
            ! tallies.
            if (t % estimator() == ESTIMATOR_TRACKLENGTH) then
              call t % set_estimator(ESTIMATOR_COLLISION)
            end if
            ! We found the derivative we were looking for; exit the do loop.
            exit
          end if
          if (j == n_tally_derivs()) then
            call fatal_error("Could not find derivative " &
                 // trim(to_str(deriv_id)) // " specified on tally " &
                 // trim(to_str(t % id())))
          end if
        end do

        deriv => tally_deriv_c(t % deriv())
        if (deriv % variable == DIFF_NUCLIDE_DENSITY &
             .or. deriv % variable == DIFF_TEMPERATURE) then
          do j = 1, t % n_nuclide_bins()
            if (has_energyout .and. t % nuclide_bins(j) == -1) then
              call fatal_error("Error on tally " // trim(to_str(t % id())) &
                   // ": Cannot use a 'nuclide_density' or 'temperature' &
                   &derivative on a tally with an outgoing energy filter and &
                   &'total' nuclide rate. Instead, tally each nuclide in the &
                   &material individually.")
              ! Note that diff tallies with these characteristics would work
              ! correctly if no tally events occur in the perturbed material
              ! (e.g. pertrubing moderator but only tallying fuel), but this
              ! case would be hard to check for by only reading inputs.
            end if
          end do
        end if
      end if

      ! If settings.xml trigger is turned on, create tally triggers
      if (trigger_on) then
        !TODO: off-by-one
        call tally_init_triggers(t % ptr, i_start + i - 1 - 1, node_tal % ptr)
      end if

      ! =======================================================================
      ! SET TALLY ESTIMATOR

      ! Check if user specified estimator
      if (check_for_node(node_tal, "estimator")) then
        temp_str = ''
        call get_node_value(node_tal, "estimator", temp_str)
        select case(trim(temp_str))
        case ('analog')
          call t % set_estimator(ESTIMATOR_ANALOG)

        case ('tracklength', 'track-length', 'pathlength', 'path-length')
          ! If the estimator was set to an analog estimator, this means the
          ! tally needs post-collision information
          if (t % estimator() == ESTIMATOR_ANALOG) then
            call fatal_error("Cannot use track-length estimator for tally " &
                 // to_str(t % id()))
          end if

          ! Set estimator to track-length estimator
          call t % set_estimator(ESTIMATOR_TRACKLENGTH)

        case ('collision')
          ! If the estimator was set to an analog estimator, this means the
          ! tally needs post-collision information
          if (t % estimator() == ESTIMATOR_ANALOG) then
            call fatal_error("Cannot use collision estimator for tally " &
                 // to_str(t % id()))
          end if

          ! Set estimator to collision estimator
          call t % set_estimator(ESTIMATOR_COLLISION)

        case default
          call fatal_error("Invalid estimator '" // trim(temp_str) &
               // "' on tally " // to_str(t % id()))
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

  subroutine read_plots_xml() bind(C)

    logical :: file_exists              ! does plots.xml file exist?
    character(MAX_LINE_LEN) :: filename ! absolute path to plots.xml
    type(XMLDocument) :: doc
    type(XMLNode) :: root

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

    call read_plots(root % ptr)

    ! Close plots XML file
    call doc % clear()

  end subroutine read_plots_xml

  subroutine read_mg_cross_sections_header() bind(C)
    integer :: i           ! loop index
    logical :: file_exists ! does mgxs.h5 exist?
    integer(HID_T) :: file_id
    character(kind=C_CHAR), pointer :: string(:)

    interface
      subroutine read_mg_cross_sections_header_c(file_id) bind(C)
        import HID_T
        integer(HID_T), value :: file_id
      end subroutine

      function path_cross_sections_c() result(ptr) bind(C)
        import C_PTR
        type(C_PTR) :: ptr
      end function
    end interface

    call c_f_pointer(path_cross_sections_c(), string, [255])
    path_cross_sections = to_f_string(string)

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

    ! Set up energy bins on C++ side
    call read_mg_cross_sections_header_c(file_id)

    ! Get the minimum and maximum energies
    call set_particle_energy_bounds(NEUTRON, &
         energy_bins(num_energy_groups + 1), energy_bins(1))

    ! Close MGXS HDF5 file
    call file_close(file_id)

  end subroutine read_mg_cross_sections_header

end module input_xml
