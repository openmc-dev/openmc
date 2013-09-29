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

    integer :: new_comm ! new mpi communicator
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
# ifdef PETSC 
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, 0, new_comm, mpi_err)
# endif

    ! assign to PETSc
# ifdef PETSC
    PETSC_COMM_WORLD = new_comm

    ! Initialize PETSc on all procs
    call PetscInitialize(PETSC_NULL_CHARACTER, mpi_err)
# endif

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
    use output,  only: write_message
    use string,  only: lower_case
    use xml_data_cmfd_t
    use, intrinsic :: ISO_FORTRAN_ENV

    integer :: ng                       ! number of energy groups
    logical :: file_exists              ! does cmfd.xml exist?
    character(MAX_LINE_LEN) :: filename ! name of input file

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
    call read_xml_file_cmfd_t(filename)

    ! Set spatial dimensions in cmfd object (structed Cartesian mesh)
    cmfd % indices(1:3) = mesh_ % dimension(1:3)

    ! Get number of energy groups or set to 1 group default
    if (associated(mesh_ % energy)) then
      ng = size(mesh_ % energy)
      if(.not.allocated(cmfd % egrid)) allocate(cmfd % egrid(ng))
      cmfd % egrid = mesh_ % energy 
      cmfd % indices(4) = ng - 1 ! sets energy group dimension
    else
      if(.not.allocated(cmfd % egrid)) allocate(cmfd % egrid(2))
      cmfd % egrid = (/0.0_8,20.0_8/)
      cmfd % indices(4) = 1 ! one energy group
    end if
    
    ! Set global albedo, these can be overwritten by coremap
    if (associated(mesh_ % albedo)) then
      cmfd % albedo = mesh_ % albedo
    else
      cmfd % albedo = (/1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
    end if

    ! Get core map overlay for a subset mesh for CMFD 
    if (associated(mesh_ % map)) then

      ! Allocate a core map with appropriate dimensions
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2), &
           cmfd % indices(3)))

      ! Check and make sure it is of correct size
      if (size(mesh_ % map) /= product(cmfd % indices(1:3))) then
        message = 'CMFD coremap not to correct dimensions'
        call fatal_error() 
      end if

      ! Reshape core map vector into (x,y,z) array
      cmfd % coremap = reshape(mesh_ % map,(cmfd % indices(1:3)))

      ! Indicate to cmfd that a core map overlay is active
      cmfd_coremap = .true.

      ! Write to stdout that a core map overlay is active
      message = "Core Map Overlay Activated"
      call write_message(10)

    end if

    ! Get normalization constant for source (default is 1.0 from XML) 
    cmfd % norm = norm_

    ! Set feedback logical
    call lower_case(feedback_)
    if (feedback_ == 'true' .or. feedback_ == '1') cmfd_feedback = .true.

    ! Set downscatter logical
    call lower_case(downscatter_)
    if (downscatter_ == 'true' .or. downscatter_ == '1') &
         cmfd_downscatter = .true.

    ! Set the solver type (default power from XML)
    cmfd_solver_type = solver_(1:10)

    ! Set convergence monitoring 
    call lower_case(snes_monitor_)
    call lower_case(ksp_monitor_)
    call lower_case(power_monitor_)
    if (snes_monitor_ == 'true' .or. snes_monitor_ == '1') &
         cmfd_snes_monitor = .true.
    if (ksp_monitor_ == 'true' .or. ksp_monitor_ == '1') &
         cmfd_ksp_monitor = .true.
    if (power_monitor_ == 'true' .or. power_monitor_ == '1') &
         cmfd_power_monitor = .true.

    ! Output logicals
    call lower_case(write_matrices_)
    if (write_matrices_ == 'true' .or. write_matrices_ == '1') &
         cmfd_write_matrices = .true.

    ! Run an adjoint calc
    call lower_case(run_adjoint_)
    if (run_adjoint_ == 'true' .or. run_adjoint_ == '1') &
         cmfd_run_adjoint = .true.

    ! Batch to begin cmfd (default is 1 from XML)
    cmfd_begin = begin_

    ! Tally during inactive batches (by default we will always tally from 1)
    call lower_case(inactive_)
    if (inactive_ == 'false' .or. inactive_ == '0') cmfd_tally_on = .false.

    ! Inactive batch flush window
    cmfd_inact_flush(1) = inactive_flush_ ! the interval of batches
    cmfd_inact_flush(2) = num_flushes_    ! number of times to do this

    ! Last flush before active batches
    cmfd_act_flush = active_flush_

    ! Get display
    cmfd_display = display_
    if (trim(cmfd_display) == 'dominance' .and. &
        trim(cmfd_solver_type) /= 'power') then
      message = 'Dominance Ratio only aviable with power iteration solver'
      call warning()
      cmfd_display = ''
    end if

    ! Create tally objects
    call create_cmfd_tally()

  end subroutine read_cmfd_xml

!===============================================================================
! CREATE_CMFD_TALLY creates the tally object for OpenMC to process for CMFD
! accleration.
! There are 3 tally types:
!   1: Only an energy in filter-> flux,total,p1 scatter
!   2: Energy in and energy out filter-> nu-scatter,nu-fission
!   3: Surface current
!===============================================================================

  subroutine create_cmfd_tally()

    use error,            only: fatal_error, warning
    use mesh_header,      only: StructuredMesh
    use string
    use tally,            only: setup_active_cmfdtallies
    use tally_header,     only: TallyObject, TallyFilter
    use tally_initialize, only: add_tallies
    use xml_data_cmfd_t

    integer :: i           ! loop counter
    integer :: n           ! size of arrays in mesh specification
    integer :: ng          ! number of energy groups (default 1)
    integer :: n_filters   ! number of filters
    integer :: i_filter_mesh ! index for mesh filter
    character(MAX_LINE_LEN) :: filename ! name of cmfd file
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter) :: filters(N_FILTER_TYPES) ! temporary filters

    ! Parse cmfd.xml file
     filename = trim(path_input) // "cmfd.xml"
     call read_xml_file_cmfd_t(filename)

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

    ! Determine number of dimensions for mesh
    n = size(mesh_ % dimension)
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
    if (any(mesh_ % dimension <= 0)) then
       message = "All entries on the <dimension> element for a tally mesh &
            &must be positive."
       call fatal_error()
    end if

    ! Read dimensions in each direction
    m % dimension = mesh_ % dimension

    ! Read mesh lower-left corner location
    if (m % n_dimension /= size(mesh_ % lower_left)) then
       message = "Number of entries on <lower_left> must be the same as &
            &the number of entries on <dimension>."
       call fatal_error()
    end if
    m % lower_left = mesh_ % lower_left

    ! Make sure either upper-right or width was specified
    if (associated(mesh_ % upper_right) .and. &
         associated(mesh_ % width)) then
       message = "Cannot specify both <upper_right> and <width> on a &
             &tally mesh."
       call fatal_error()
    end if

    ! Make sure either upper-right or width was specified
    if (.not. associated(mesh_ % upper_right) .and. &
         .not. associated(mesh_ % width)) then
       message = "Must specify either <upper_right> and <width> on a &
            &tally mesh."
       call fatal_error()
    end if

    if (associated(mesh_ % width)) then
       ! Check to ensure width has same dimensions
       if (size(mesh_ % width) /= size(mesh_ % lower_left)) then
          message = "Number of entries on <width> must be the same as the &
               &number of entries on <lower_left>."
          call fatal_error()
       end if

       ! Check for negative widths
       if (any(mesh_ % width < ZERO)) then
          message = "Cannot have a negative <width> on a tally mesh."
          call fatal_error()
       end if

       ! Set width and upper right coordinate
       m % width = mesh_ % width
       m % upper_right = m % lower_left + m % dimension * m % width

    elseif (associated(mesh_ % upper_right)) then
       ! Check to ensure width has same dimensions
       if (size(mesh_ % upper_right) /= size(mesh_ % lower_left)) then
          message = "Number of entries on <upper_right> must be the same as &
               &the number of entries on <lower_left>."
          call fatal_error()
       end if

       ! Check that upper-right is above lower-left
       if (any(mesh_ % upper_right < mesh_ % lower_left)) then
          message = "The <upper_right> coordinates must be greater than the &
               &<lower_left> coordinates on a tally mesh."
          call fatal_error()
       end if

       ! Set width and upper right coordinate
       m % upper_right = mesh_ % upper_right
       m % width = (m % upper_right - m % lower_left) / m % dimension
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
      call lower_case(reset_)
      if (reset_ == 'true' .or. reset_ == '1') t % reset = .true.

      ! Set up mesh filter
      n_filters = 1
      filters(n_filters) % type = FILTER_MESH
      filters(n_filters) % n_bins = product(m % dimension)
      allocate(filters(n_filters) % int_bins(1))
      filters(n_filters) % int_bins(1) = n_user_meshes + 1
      t % find_filter(FILTER_MESH) = n_filters

      ! Read and set incoming energy mesh filter
      if (associated(mesh_ % energy)) then
        n_filters = n_filters + 1
        filters(n_filters) % type = FILTER_ENERGYIN
        ng = size(mesh_ % energy)
        filters(n_filters) % n_bins = ng - 1
        allocate(filters(n_filters) % real_bins(ng))
        filters(n_filters) % real_bins = mesh_ % energy
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

        ! Read and set outgoing energy mesh filter
        if (associated(mesh_ % energy)) then
          n_filters = n_filters + 1
          filters(n_filters) % type = FILTER_ENERGYOUT
          ng = size(mesh_ % energy)
          filters(n_filters) % n_bins = ng - 1
          allocate(filters(n_filters) % real_bins(ng))
          filters(n_filters) % real_bins = mesh_ % energy
          t % find_filter(FILTER_ENERGYOUT) = n_filters
        end if

        ! Allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! Deallocate filters bins array
        if (associated(mesh_ % energy)) &
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
      if (associated(mesh_ % energy)) deallocate(filters(2) % real_bins)

    end do

    ! Put cmfd tallies into active tally array and turn tallies on
    call setup_active_cmfdtallies()
    tallies_on = .true.

  end subroutine create_cmfd_tally

end module cmfd_input
