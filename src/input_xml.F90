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
  use output,           only: title, header
  use photon_header
  use random_lcg,       only: prn
  use surface_header
  use set_header,       only: SetChar
  use settings
  use stl_vector,       only: VectorInt, VectorReal, VectorChar
  use string,           only: to_lower, to_str, str_to_int, str_to_real, &
                              starts_with, ends_with, split_string, &
                              zero_padded, to_c_string
  use tally
  use tally_header,     only: openmc_extend_tallies
  use tally_derivative_header
  use tally_filter_header
  use tally_filter
  use trigger_header
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
!$omp parallel
!$omp critical (ReadTallyDeriv)
      call tally_derivs(i) % from_xml(node_deriv_list(i))
!$omp end critical (ReadTallyDeriv)
!$omp end parallel

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
            word = sarray(j)

            ! Search through nuclides
            k = nuclide_map_get(to_c_string(word))
            if (k == -1) then
              call fatal_error("Could not find the nuclide " &
                   // trim(word) // " specified in tally " &
                   // trim(to_str(t % id)) // " in any material.")
            end if

            ! Set bin to index in nuclides array
            t % nuclide_bins(j) = k
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
    energy_min(NEUTRON) = energy_bins(num_energy_groups + 1)
    energy_max(NEUTRON) = energy_bins(1)
    call set_particle_energy_bounds(NEUTRON, energy_min(NEUTRON), &
         energy_max(NEUTRON))

    ! Close MGXS HDF5 file
    call file_close(file_id)

  end subroutine read_mg_cross_sections_header

end module input_xml
