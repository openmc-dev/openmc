module cmfd_input

  use global

#ifdef PETSC
  use petscsys
#endif

  implicit none
  private
  public :: configure_cmfd 

contains

!===============================================================================
! CONFIGURE_CMFD initializes PETSc and CMFD parameters
!===============================================================================

  subroutine configure_cmfd()

    use cmfd_header,  only: allocate_cmfd

#ifdef PETSC
    integer :: new_comm ! new mpi communicator
#endif
    integer :: color    ! color group of processor

    ! Read in cmfd input file
    call read_cmfd_xml()

    ! Assign color
    if (master) then
      color = 1
    else
      color = 2
    end if

    ! Split up procs
#ifdef PETSC
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, 0, new_comm, mpi_err)
#endif

    ! assign to PETSc
#ifdef PETSC
    PETSC_COMM_WORLD = new_comm

    ! Initialize PETSc on all procs
    call PetscInitialize(PETSC_NULL_CHARACTER, mpi_err)
#endif

    ! Initialize timers
    call time_cmfd % reset()
    call time_cmfdbuild % reset()
    call time_cmfdsolve % reset()

    ! Allocate cmfd object
    call allocate_cmfd(cmfd, n_batches)

  end subroutine configure_cmfd

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_cmfd_xml()

    use error,   only: fatal_error, warning
    use global
    use output,  only: write_message
    use string,  only: lower_case
    use xml_interface
    use, intrinsic :: ISO_FORTRAN_ENV

    integer :: ng
    integer, allocatable :: iarray(:)
    logical :: file_exists ! does cmfd.xml exist?
    logical :: found
    character(MAX_LINE_LEN) :: filename
    character(MAX_LINE_LEN) :: temp_str
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_mesh => null()

    ! Read cmfd input file
    filename = trim(path_input) // "cmfd.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! CMFD is optional unless it is in on from settings
      if (cmfd_on) then
        message = "No CMFD XML file, '" // trim(filename) // "' does not exist!"
        call fatal_error()
      end if
      return
    else

      ! Tell user
      message = "Reading CMFD XML file..."
      call write_message(5)

    end if

    ! Parse cmfd.xml file
    call open_xmldoc(doc, filename)

    ! Get pointer to mesh XML node
    call get_node_ptr(doc, "mesh", node_mesh, found = found)

    ! Check if mesh is there
    if (.not.found) then
      message = "No CMFD mesh specified in CMFD XML file."
      call fatal_error()
    end if

    ! Set spatial dimensions in cmfd object
    call get_node_array(node_mesh, "dimension", cmfd % indices(1:3))

    ! Get number of energy groups
    if (check_for_node(node_mesh, "energy")) then
      ng = get_arraysize_double(node_mesh, "energy")
      if(.not.allocated(cmfd%egrid)) allocate(cmfd%egrid(ng))
      call get_node_array(node_mesh, "energy", cmfd%egrid)
      cmfd % indices(4) = ng - 1 ! sets energy group dimension
    else
      if(.not.allocated(cmfd % egrid)) allocate(cmfd % egrid(2))
      cmfd % egrid = (/0.0_8,20.0_8/)
      cmfd % indices(4) = 1 ! one energy group
    end if
    
    ! Set global albedo
    if (check_for_node(node_mesh, "albedo")) then
      call get_node_array(node_mesh, "albedo", cmfd % albedo)
    else
      cmfd % albedo = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
    end if

    ! Get acceleration map
    if (check_for_node(node_mesh, "map")) then
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2), &
           cmfd % indices(3)))
      if (get_arraysize_integer(node_mesh, "map") /= &
          product(cmfd % indices(1:3))) then
        message = 'FATAL==>CMFD coremap not to correct dimensions'
        call fatal_error() 
      end if
      allocate(iarray(get_arraysize_integer(node_mesh, "map")))
      call get_node_array(node_mesh, "map", iarray)
      cmfd % coremap = reshape(iarray,(cmfd % indices(1:3)))
      cmfd_coremap = .true.
      deallocate(iarray)
    end if

    ! Check for normalization constant
    if (check_for_node(doc, "norm")) then
      call get_node_value(doc, "norm", cmfd % norm)
    end if

    ! Set feedback logical
    if (check_for_node(doc, "feedback")) then
      call get_node_value(doc, "feedback", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_feedback = .true.
    end if

    ! Set downscatter logical
    if (check_for_node(doc, "downscatter")) then
      call get_node_value(doc, "downscatter", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_downscatter = .true.
    end if

    ! Set the solver type
    if (check_for_node(doc, "solver")) &
      call get_node_value(doc, "solver", cmfd_solver_type)

    ! Set monitoring
    if (check_for_node(doc, "snes_monitor")) then
      call get_node_value(doc, "snes_monitor", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_snes_monitor = .true.
    end if
    if (check_for_node(doc, "ksp_monitor")) then
      call get_node_value(doc, "ksp_monitor", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_ksp_monitor = .true.
    end if
    if (check_for_node(doc, "power_monitor")) then
      call get_node_value(doc, "power_monitor", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_power_monitor = .true.
    end if

    ! Output logicals
    if (check_for_node(doc, "write_matrices")) then
      call get_node_value(doc, "write_matices", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_write_matrices = .true.
    end if

    ! Run an adjoint calc
    if (check_for_node(doc, "run_adjoint")) then
      call get_node_value(doc, "run_adjoint", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
        cmfd_run_adjoint = .true.
    end if

    ! Batch to begin cmfd
    if (check_for_node(doc, "begin")) &
      call get_node_value(doc, "begin", cmfd_begin) 

    ! Tally during inactive batches
    if (check_for_node(doc, "inactive")) then
      call get_node_value(doc, "inactive", temp_str)
      call lower_case(temp_str)
      if (trim(temp_str) == 'false' .or. trim(temp_str) == '0') &
        cmfd_tally_on = .false.
    end if

    ! Inactive batch flush window
    if (check_for_node(doc, "inactive_flush")) &
      call get_node_value(doc, "inactive_flush", cmfd_inact_flush(1))
    if (check_for_node(doc, "num_flushes")) &
      call get_node_value(doc, "num_flushes", cmfd_inact_flush(2))

    ! Last flush before active batches
    if (check_for_node(doc, "active_flush")) &
      call get_node_value(doc, "active_flush", cmfd_act_flush)

    ! Get display
    if (check_for_node(doc, "display")) &
      call get_node_value(doc, "display", cmfd_display)
    if (trim(cmfd_display) == 'dominance' .and. &
        trim(cmfd_solver_type) /= 'power') then
      message = 'Dominance Ratio only aviable with power iteration solver'
      call warning()
      cmfd_display = ''
    end if

    ! Create tally objects
    call create_cmfd_tally(doc)

    ! Close CMFD XML file
    call close_xmldoc(doc)

  end subroutine read_cmfd_xml

!===============================================================================
! CREATE_CMFD_TALLY creates the tally object for OpenMC to process for CMFD
! accleration.
! There are 3 tally types:
!   1: Only an energy in filter-> flux,total,p1 scatter
!   2: Energy in and energy out filter-> nu-scatter,nu-fission
!   3: Surface current
!===============================================================================

  subroutine create_cmfd_tally(doc)

    use constants,        only: MAX_LINE_LEN
    use error,            only: fatal_error, warning
    use mesh_header,      only: StructuredMesh
    use string
    use tally,            only: setup_active_cmfdtallies
    use tally_header,     only: TallyObject, TallyFilter
    use tally_initialize, only: add_tallies
    use xml_interface

    type(Node), pointer :: doc ! pointer to XML doc info

    character(MAX_LINE_LEN) :: temp_str ! temp string
    integer :: i           ! loop counter
    integer :: n           ! size of arrays in mesh specification
    integer :: ng          ! number of energy groups (default 1)
    integer :: n_filters   ! number of filters
    integer :: i_filter_mesh ! index for mesh filter
    integer :: iarray3(3) ! temp integer array
    real(8) :: rarray3(3) ! temp double array
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter) :: filters(N_FILTER_TYPES) ! temporary filters
    type(Node), pointer :: node_mesh => null()

    ! Set global variables if they are 0 (this can happen if there is no tally
    ! file)
    if (n_meshes == 0) n_meshes = n_user_meshes + n_cmfd_meshes

    ! Allocate mesh
    if (.not. allocated(meshes)) allocate(meshes(n_meshes))
    m => meshes(n_user_meshes+1)

    ! Set mesh id
    m % id = n_user_meshes + 1 

    ! Set mesh type to rectangular
    m % type = LATTICE_RECT

    ! Get pointer to mesh XML node
    call get_node_ptr(doc, "mesh", node_mesh)

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

      ! Set upper right coordinate and width
      m % upper_right = rarray3(1:n)
      m % width = (m % upper_right - m % lower_left) / real(m % dimension, 8)
    end if

    ! Set volume fraction
    m % volume_frac = ONE/real(product(m % dimension),8)

    ! Add mesh to dictionary
    call mesh_dict % add_key(m % id, n_user_meshes + 1)

    ! Allocate tallies
    call add_tallies("cmfd", n_cmfd_tallies)

    ! Begin loop around tallies
    do i = 1, n_cmfd_tallies

      ! Point t to tally variable
      t => cmfd_tallies(i)

      ! Set reset property
      if (check_for_node(doc, "reset")) then
        call get_node_value(doc, "reset", temp_str)
        call lower_case(temp_str)
        if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') &
          t % reset = .true.
      end if

      ! Set up mesh filter
      n_filters = 1
      filters(n_filters) % type = FILTER_MESH
      filters(n_filters) % n_bins = product(m % dimension)
      allocate(filters(n_filters) % int_bins(1))
      filters(n_filters) % int_bins(1) = n_user_meshes + 1
      t % find_filter(FILTER_MESH) = n_filters

      ! Read and set incoming energy mesh filter
      if (check_for_node(node_mesh, "energy")) then
        n_filters = n_filters + 1
        filters(n_filters) % type = FILTER_ENERGYIN
        ng = get_arraysize_double(node_mesh, "energy") 
        filters(n_filters) % n_bins = ng - 1
        allocate(filters(n_filters) % real_bins(ng))
        call get_node_array(node_mesh, "energy", &
             filters(n_filters) % real_bins)
        t % find_filter(FILTER_ENERGYIN) = n_filters
      end if

      ! Set number of nucilde bins
      allocate(t % nuclide_bins(1))
      t % nuclide_bins(1) = -1
      t % n_nuclide_bins = 1

      ! Record tally id which is equivalent to loop number
      t % id = i_cmfd_tallies + i

      if (i == 1) then

        ! Set label
        t % label = "CMFD flux, total, scatter-1"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Set tally type to volume
        t % type = TALLY_VOLUME

        ! Allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! Allocate scoring bins 
        allocate(t % score_bins(3))
        t % n_score_bins = 3
        t % n_user_score_bins = 3

        ! Allocate scattering order data
        allocate(t % scatt_order(3))
        t % scatt_order = 0
        
        ! Set macro_bins
        t % score_bins(1)  = SCORE_FLUX
        t % score_bins(2)  = SCORE_TOTAL
        t % score_bins(3)  = SCORE_SCATTER_N
        t % scatt_order(3) = 1

      else if (i == 2) then

        ! Set label
        t % label = "CMFD neutron production"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Set tally type to volume
        t % type = TALLY_VOLUME

        ! read and set outgoing energy mesh filter
        if (check_for_node(node_mesh, "energy")) then
          n_filters = n_filters + 1
          filters(n_filters) % type = FILTER_ENERGYOUT
          ng = get_arraysize_double(node_mesh, "energy")
          filters(n_filters) % n_bins = ng - 1
          allocate(filters(n_filters) % real_bins(ng))
          call get_node_array(node_mesh, "energy", &
               filters(n_filters) % real_bins)
          t % find_filter(FILTER_ENERGYOUT) = n_filters
        end if

        ! Allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! deallocate filters bins array
        if (check_for_node(node_mesh, "energy")) &
             deallocate(filters(n_filters) % real_bins)

        ! Allocate macro reactions
        allocate(t % score_bins(2))
        t % n_score_bins = 2
        t % n_user_score_bins = 2

        ! Allocate scattering order data
        allocate(t % scatt_order(2))
        t % scatt_order = 0

        ! Set macro_bins
        t % score_bins(1) = SCORE_NU_SCATTER
        t % score_bins(2) = SCORE_NU_FISSION

      else if (i == 3) then

        ! Set label
        t % label = "CMFD surface currents"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Add extra filter for surface
        n_filters = n_filters + 1
        filters(n_filters) % type = FILTER_SURFACE
        filters(n_filters) % n_bins = 2 * m % n_dimension
        allocate(filters(n_filters) % int_bins(2 * m % n_dimension))
        if (m % n_dimension == 2) then
          filters(n_filters) % int_bins = (/ IN_RIGHT, OUT_RIGHT, IN_FRONT, &
               OUT_FRONT /)
        elseif (m % n_dimension == 3) then
          filters(n_filters) % int_bins = (/ IN_RIGHT, OUT_RIGHT, IN_FRONT, &
               OUT_FRONT, IN_TOP, OUT_TOP /)
        end if
        t % find_filter(FILTER_SURFACE) = n_filters

        ! Allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! Deallocate filters bins array
        deallocate(filters(n_filters) % int_bins)

        ! Allocate macro reactions
        allocate(t % score_bins(1))
        t % n_score_bins = 1
        t % n_user_score_bins = 1

        ! Allocate scattering order data
        allocate(t % scatt_order(1))
        t % scatt_order = 0

        ! Set macro bins
        t % score_bins(1) = SCORE_CURRENT
        t % type = TALLY_SURFACE_CURRENT

        ! We need to increase the dimension by one since we also need
        ! currents coming into and out of the boundary mesh cells.
        i_filter_mesh = t % find_filter(FILTER_MESH)
        t % filters(i_filter_mesh) % n_bins = product(m % dimension + 1)

      end if

      ! Deallocate filter bins
      deallocate(filters(1) % int_bins)
      if (check_for_node(node_mesh, "energy")) &
        deallocate(filters(2) % real_bins)

    end do

    ! Put cmfd tallies into active tally array and turn tallies on
!$omp parallel
    call setup_active_cmfdtallies()
!$omp end parallel
    tallies_on = .true.

  end subroutine create_cmfd_tally

end module cmfd_input
