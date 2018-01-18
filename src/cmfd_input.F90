module cmfd_input

  use, intrinsic :: ISO_C_BINDING

  use cmfd_header
  use mesh_header, only: mesh_dict
  use mgxs_header, only: energy_bins
  use tally
  use tally_header
  use timer_header

  implicit none
  private
  public :: configure_cmfd

contains

!===============================================================================
! CONFIGURE_CMFD initializes CMFD parameters
!===============================================================================

  subroutine configure_cmfd()

    ! Read in cmfd input file
    call read_cmfd_xml()

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

    use constants, only: ZERO, ONE
    use error,     only: fatal_error, warning, write_message
    use string,    only: to_lower
    use xml_interface
    use, intrinsic :: ISO_FORTRAN_ENV

    integer :: i, g
    integer :: ng
    integer :: n_params
    integer, allocatable :: iarray(:)
    integer, allocatable :: int_array(:)
    logical :: file_exists ! does cmfd.xml exist?
    logical :: found
    character(MAX_LINE_LEN) :: filename
    real(8) :: gs_tol(2)
    type(XMLDocument) :: doc
    type(XMLNode) :: root
    type(XMLNode) :: node_mesh

    ! Read cmfd input file
    filename = trim(path_input) // "cmfd.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      ! CMFD is optional unless it is in on from settings
      if (cmfd_run) then
        call fatal_error("No CMFD XML file, '" // trim(filename) // "' does not&
             & exist!")
      end if
      return
    else

      ! Tell user
      call write_message("Reading CMFD XML file...", 5)

    end if

    ! Parse cmfd.xml file
    call doc % load_file(filename)
    root = doc % document_element()

    ! Get pointer to mesh XML node
    node_mesh = root % child("mesh")

    ! Check if mesh is there
    if (.not. node_mesh % associated()) then
      call fatal_error("No CMFD mesh specified in CMFD XML file.")
    end if

    ! Set spatial dimensions in cmfd object
    call get_node_array(node_mesh, "dimension", cmfd % indices(1:3))

    ! Get number of energy groups
    if (check_for_node(node_mesh, "energy")) then
      ng = node_word_count(node_mesh, "energy")
      if(.not. allocated(cmfd%egrid)) allocate(cmfd%egrid(ng))
      call get_node_array(node_mesh, "energy", cmfd%egrid)
      cmfd % indices(4) = ng - 1 ! sets energy group dimension
      ! If using MG mode, check to see if these egrid points at least match
      ! the MG Data breakpoints
      if (.not. run_CE) then
        do i = 1, ng
          found = .false.
          do g = 1, num_energy_groups + 1
            if (cmfd % egrid(i) == energy_bins(g)) then
              found = .true.
              exit
            end if
          end do
          if (.not. found) then
            call fatal_error("CMFD energy mesh boundaries must align with&
                             & boundaries of multi-group data!")
          end if
        end do
      end if
    else
      if(.not.allocated(cmfd % egrid)) allocate(cmfd % egrid(2))
      cmfd % egrid = [ ZERO, energy_max_neutron ]
      cmfd % indices(4) = 1 ! one energy group
    end if

    ! Set global albedo
    if (check_for_node(node_mesh, "albedo")) then
      call get_node_array(node_mesh, "albedo", cmfd % albedo)
    else
      cmfd % albedo = [ ONE, ONE, ONE, ONE, ONE, ONE ]
    end if

    ! Get acceleration map
    if (check_for_node(node_mesh, "map")) then
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2), &
           cmfd % indices(3)))
      if (node_word_count(node_mesh, "map") /= &
           product(cmfd % indices(1:3))) then
        call fatal_error('CMFD coremap not to correct dimensions')
      end if
      allocate(iarray(node_word_count(node_mesh, "map")))
      call get_node_array(node_mesh, "map", iarray)
      cmfd % coremap = reshape(iarray,(cmfd % indices(1:3)))
      cmfd_coremap = .true.
      deallocate(iarray)
    end if

    ! Check for normalization constant
    if (check_for_node(root, "norm")) then
      call get_node_value(root, "norm", cmfd % norm)
    end if

    ! Set feedback logical
    if (check_for_node(root, "feedback")) then
      call get_node_value(root, "feedback", cmfd_feedback)
    end if

    ! Set downscatter logical
    if (check_for_node(root, "downscatter")) then
      call get_node_value(root, "downscatter", cmfd_downscatter)
    end if

    ! Reset dhat parameters
    if (check_for_node(root, "dhat_reset")) then
      call get_node_value(root, "dhat_reset", dhat_reset)
    end if

    ! Set monitoring
    if (check_for_node(root, "power_monitor")) then
      call get_node_value(root, "power_monitor", cmfd_power_monitor)
    end if

    ! Output logicals
    if (check_for_node(root, "write_matrices")) then
      call get_node_value(root, "write_matrices", cmfd_write_matrices)
    end if

    ! Run an adjoint calc
    if (check_for_node(root, "run_adjoint")) then
      call get_node_value(root, "run_adjoint", cmfd_run_adjoint)
    end if

    ! Batch to begin cmfd
    if (check_for_node(root, "begin")) &
         call get_node_value(root, "begin", cmfd_begin)

    ! Check for cmfd tally resets
    if (check_for_node(root, "tally_reset")) then
      n_cmfd_resets = node_word_count(root, "tally_reset")
    else
      n_cmfd_resets = 0
    end if
    if (n_cmfd_resets > 0) then
      allocate(int_array(n_cmfd_resets))
      call get_node_array(root, "tally_reset", int_array)
      do i = 1, n_cmfd_resets
        call cmfd_reset % add(int_array(i))
      end do
      deallocate(int_array)
    end if

    ! Get display
    if (check_for_node(root, "display")) &
         call get_node_value(root, "display", cmfd_display)

    ! Read in spectral radius estimate and tolerances
    if (check_for_node(root, "spectral")) &
         call get_node_value(root, "spectral", cmfd_spectral)
    if (check_for_node(root, "shift")) &
         call get_node_value(root, "shift", cmfd_shift)
    if (check_for_node(root, "ktol")) &
         call get_node_value(root, "ktol", cmfd_ktol)
    if (check_for_node(root, "stol")) &
         call get_node_value(root, "stol", cmfd_stol)
    if (check_for_node(root, "gauss_seidel_tolerance")) then
      n_params = node_word_count(root, "gauss_seidel_tolerance")
      if (n_params /= 2) then
        call fatal_error('Gauss Seidel tolerance is not 2 parameters &
                   &(absolute, relative).')
      end if
      call get_node_array(root, "gauss_seidel_tolerance", gs_tol)
      cmfd_atoli = gs_tol(1)
      cmfd_rtoli = gs_tol(2)
    end if

    ! Create tally objects
    call create_cmfd_tally(root)

    ! Close CMFD XML file
    call doc % clear()

  end subroutine read_cmfd_xml

!===============================================================================
! CREATE_CMFD_TALLY creates the tally object for OpenMC to process for CMFD
! accleration.
! There are 3 tally types:
!   1: Only an energy in filter-> flux,total,p1 scatter
!   2: Energy in and energy out filter-> nu-scatter,nu-fission
!   3: Mesh current
!===============================================================================

  subroutine create_cmfd_tally(root)

    use constants,        only: MAX_LINE_LEN
    use error,            only: fatal_error, warning
    use mesh_header,      only: RegularMesh, openmc_extend_meshes
    use string
    use tally,            only: openmc_tally_set_type
    use tally_header,     only: openmc_extend_tallies
    use tally_filter_header
    use tally_filter
    use xml_interface

    type(XMLNode), intent(in) :: root ! XML root element

    logical :: energy_filters
    integer :: i           ! loop counter
    integer :: n           ! size of arrays in mesh specification
    integer(C_INT32_T) :: ng  ! number of energy groups (default 1)
    integer :: n_filter    ! number of filters
    integer :: i_start, i_end
    integer :: i_filt_start, i_filt_end
    integer(C_INT32_T), allocatable :: filter_indices(:)
    integer(C_INT) :: err
    integer :: i_filt     ! index in filters array
    integer :: filt_id
    integer :: iarray3(3) ! temp integer array
    real(8) :: rarray3(3) ! temp double array
    real(C_DOUBLE), allocatable :: energies(:)
    type(RegularMesh), pointer :: m
    type(XMLNode) :: node_mesh

    err = openmc_extend_meshes(1, i_start)

    ! Allocate mesh
    cmfd_mesh => meshes(i_start)
    m => meshes(i_start)

    ! Set mesh id
    m % id = i_start

    ! Set mesh type to rectangular
    m % type = LATTICE_RECT

    ! Get pointer to mesh XML node
    node_mesh = root % child("mesh")

    ! Determine number of dimensions for mesh
    n = node_word_count(node_mesh, "dimension")
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
      call fatal_error("All entries on the <dimension> element for a tally mesh&
           & must be positive.")
    end if

    ! Read dimensions in each direction
    m % dimension = iarray3(1:n)

    ! Read mesh lower-left corner location
    if (m % n_dimension /= node_word_count(node_mesh, "lower_left")) then
      call fatal_error("Number of entries on <lower_left> must be the same as &
           &the number of entries on <dimension>.")
    end if
    call get_node_array(node_mesh, "lower_left", m % lower_left)

    ! Make sure both upper-right or width were specified
    if (check_for_node(node_mesh, "upper_right") .and. &
         check_for_node(node_mesh, "width")) then
      call fatal_error("Cannot specify both <upper_right> and <width> on a &
           &tally mesh.")
    end if

    ! Make sure either upper-right or width was specified
    if (.not.check_for_node(node_mesh, "upper_right") .and. &
         .not.check_for_node(node_mesh, "width")) then
      call fatal_error("Must specify either <upper_right> and <width> on a &
           &tally mesh.")
    end if

    if (check_for_node(node_mesh, "width")) then
      ! Check to ensure width has same dimensions
      if (node_word_count(node_mesh, "width") /= &
           node_word_count(node_mesh, "lower_left")) then
        call fatal_error("Number of entries on <width> must be the same as the &
             &number of entries on <lower_left>.")
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
        call fatal_error("Number of entries on <upper_right> must be the same &
             &as the number of entries on <lower_left>.")
      end if

      ! Check that upper-right is above lower-left
      call get_node_array(node_mesh, "upper_right", rarray3(1:n))
      if (any(rarray3(1:n) < m % lower_left)) then
        call fatal_error("The <upper_right> coordinates must be greater than &
             &the <lower_left> coordinates on a tally mesh.")
      end if

      ! Set upper right coordinate and width
      m % upper_right = rarray3(1:n)
      m % width = (m % upper_right - m % lower_left) / real(m % dimension, 8)
    end if

    ! Set volume fraction
    m % volume_frac = ONE/real(product(m % dimension),8)

    ! Add mesh to dictionary
    call mesh_dict % set(m % id, i_start)

    ! Determine number of filters
    energy_filters = check_for_node(node_mesh, "energy")
    n = merge(5, 3, energy_filters)

    ! Extend filters array so we can add CMFD filters
    err = openmc_extend_filters(n, i_filt_start, i_filt_end)

    ! Set up mesh filter
    i_filt = i_filt_start
    err = openmc_filter_set_type(i_filt, C_CHAR_'mesh' // C_NULL_CHAR)
    call openmc_get_filter_next_id(filt_id)
    err = openmc_filter_set_id(i_filt, filt_id)
    err = openmc_mesh_filter_set_mesh(i_filt, i_start)

    if (energy_filters) then
      ! Read and set incoming energy mesh filter
      i_filt = i_filt + 1
      err = openmc_filter_set_type(i_filt, C_CHAR_'energy' // C_NULL_CHAR)
      call openmc_get_filter_next_id(filt_id)
      err = openmc_filter_set_id(i_filt, filt_id)

      ! Get energies and set bins
      ng = node_word_count(node_mesh, "energy")
      allocate(energies(ng))
      call get_node_array(node_mesh, "energy", energies)
      err = openmc_energy_filter_set_bins(i_filt, ng, energies)

      ! Read and set outgoing energy mesh filter
      i_filt = i_filt + 1
      err = openmc_filter_set_type(i_filt, C_CHAR_'energyout' // C_NULL_CHAR)
      call openmc_get_filter_next_id(filt_id)
      err = openmc_filter_set_id(i_filt, filt_id)
      err = openmc_energy_filter_set_bins(i_filt, ng, energies)
    end if

    ! Duplicate the mesh filter for the mesh current tally since other
    ! tallies use this filter and we need to change the dimension
    i_filt = i_filt + 1
    err = openmc_filter_set_type(i_filt, C_CHAR_'mesh' // C_NULL_CHAR)
    call openmc_get_filter_next_id(filt_id)
    err = openmc_filter_set_id(i_filt, filt_id)
    err = openmc_mesh_filter_set_mesh(i_filt, i_start)

    ! We need to increase the dimension by one since we also need
    ! currents coming into and out of the boundary mesh cells.
    filters(i_filt) % obj % n_bins = product(m % dimension + 1)

    ! Set up surface filter
    i_filt = i_filt + 1
    allocate(SurfaceFilter :: filters(i_filt) % obj)
    select type(filt => filters(i_filt) % obj)
    type is(SurfaceFilter)
      filt % id = i_filt
      filt % n_bins = 4 * m % n_dimension
      allocate(filt % surfaces(4 * m % n_dimension))
      if (m % n_dimension == 2) then
        filt % surfaces = (/ OUT_LEFT, IN_LEFT, IN_RIGHT, OUT_RIGHT, &
             OUT_BACK, IN_BACK, IN_FRONT, OUT_FRONT /)
      elseif (m % n_dimension == 3) then
        filt % surfaces = (/ OUT_LEFT, IN_LEFT, IN_RIGHT, OUT_RIGHT, &
             OUT_BACK, IN_BACK, IN_FRONT, OUT_FRONT, &
             OUT_BOTTOM, IN_BOTTOM, IN_TOP, OUT_TOP /)
      end if
      filt % current = .true.
      ! Add filter to dictionary
      call filter_dict % set(filt % id, i_filt)
    end select

    ! Initialize filters
    do i = i_filt_start, i_filt_end
      select type (filt => filters(i) % obj)
      type is (SurfaceFilter)
        ! Don't do anything
      class default
        call filt % initialize()
      end select
    end do

    ! Allocate tallies
    err = openmc_extend_tallies(3, i_start, i_end)
    cmfd_tallies => tallies(i_start:i_end)

    ! Begin loop around tallies
    do i = 1, size(cmfd_tallies)
      ! Allocate tally
      err = openmc_tally_set_type(i_start + i - 1, C_CHAR_'generic' // C_NULL_CHAR)

      ! Point t to tally variable
      associate (t => cmfd_tallies(i) % obj)

      ! Set reset property
      if (check_for_node(root, "reset")) then
        call get_node_value(root, "reset", t % reset)
      end if

      ! Set the incoming energy mesh filter index in the tally find_filter
      ! array
      n_filter = 1
      if (energy_filters) then
        n_filter = n_filter + 1
      end if

      ! Set number of nucilde bins
      allocate(t % nuclide_bins(1))
      t % nuclide_bins(1) = -1
      t % n_nuclide_bins = 1

      ! Record tally id which is equivalent to loop number
      t % id = i_start + i - 1

      if (i == 1) then

        ! Set name
        t % name = "CMFD flux, total, scatter-1"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Set tally type to volume
        t % type = TALLY_VOLUME

        ! Allocate and set filters
        allocate(filter_indices(n_filter))
        filter_indices(1) = i_filt_start
        if (energy_filters) then
          filter_indices(2) = i_filt_start + 1
        end if
        err = openmc_tally_set_filters(i_start + i - 1, n_filter, filter_indices)
        deallocate(filter_indices)

        ! Allocate scoring bins
        allocate(t % score_bins(3))
        t % n_score_bins = 3
        t % n_user_score_bins = 3

        ! Allocate scattering order data
        allocate(t % moment_order(3))
        t % moment_order = 0

        ! Set macro_bins
        t % score_bins(1)  = SCORE_FLUX
        t % score_bins(2)  = SCORE_TOTAL
        t % score_bins(3)  = SCORE_SCATTER_N
        t % moment_order(3) = 1

      else if (i == 2) then

        ! Set name
        t % name = "CMFD neutron production"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Set tally type to volume
        t % type = TALLY_VOLUME

        ! Set the incoming energy mesh filter index in the tally find_filter
        ! array
        if (energy_filters) then
          n_filter = n_filter + 1
        end if

        ! Allocate and set indices in filters array
        allocate(filter_indices(n_filter))
        filter_indices(1) = i_filt_start
        if (energy_filters) then
          filter_indices(2) = i_filt_start + 1
          filter_indices(3) = i_filt_start + 2
        end if
        err = openmc_tally_set_filters(i_start + i - 1, n_filter, filter_indices)
        deallocate(filter_indices)

        ! Allocate macro reactions
        allocate(t % score_bins(2))
        t % n_score_bins = 2
        t % n_user_score_bins = 2

        ! Allocate scattering order data
        allocate(t % moment_order(2))
        t % moment_order = 0

        ! Set macro_bins
        t % score_bins(1) = SCORE_NU_SCATTER
        t % score_bins(2) = SCORE_NU_FISSION

      else if (i == 3) then

        ! Set name
        t % name = "CMFD surface currents"

        ! Set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! Set the surface filter index in the tally find_filter array
        n_filter = n_filter + 1

        ! Allocate and set filters
        allocate(filter_indices(n_filter))
        filter_indices(1) = i_filt_end - 1
        filter_indices(n_filter) = i_filt_end
        if (energy_filters) then
          filter_indices(2) = i_filt_start + 1
        end if
        err = openmc_tally_set_filters(i_start + i - 1, n_filter, filter_indices)
        deallocate(filter_indices)

        ! Allocate macro reactions
        allocate(t % score_bins(1))
        t % n_score_bins = 1
        t % n_user_score_bins = 1

        ! Allocate scattering order data
        allocate(t % moment_order(1))
        t % moment_order = 0

        ! Set macro bins
        t % score_bins(1) = SCORE_CURRENT
        t % type = TALLY_MESH_CURRENT
      end if

      ! Make CMFD tallies active from the start
      t % active = .true.

      end associate
    end do

  end subroutine create_cmfd_tally

end module cmfd_input
