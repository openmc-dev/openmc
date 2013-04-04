module cmfd_input

  implicit none
  private
  public :: configure_cmfd 

contains

!===============================================================================
! CONFIGURE_CMFD
!===============================================================================

  subroutine configure_cmfd()

# ifdef PETSC
    use cmfd_message_passing,   only: petsc_init_mpi
# endif

    ! read in cmfd input file
    call read_cmfd_xml()

    ! initialize petsc on mpi
# ifdef PETSC
    call petsc_init_mpi()
# endif

  end subroutine configure_cmfd

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_cmfd_xml()
    
    use error,   only: fatal_error
    use global
    use output,  only: write_message
    use string,  only: lower_case
    use xml_data_cmfd_t
    use, intrinsic :: ISO_FORTRAN_ENV

    integer :: ng
    logical :: file_exists ! does cmfd.xml exist?
    character(MAX_LINE_LEN) :: filename

    ! read cmfd infput file
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

      ! tell user
      message = "Reading CMFD XML file..."
      call write_message(5)

    end if

    ! parse cmfd.xml file
    call read_xml_file_cmfd_t(filename)

    ! set spatial dimensions in cmfd object
    cmfd % indices(1:3) = mesh_ % dimension(1:3) ! sets spatial dimensions

    ! get number of energy groups
    if (associated(mesh_ % energy)) then
      ng = size(mesh_ % energy)
      if(.not.allocated(cmfd%egrid)) allocate(cmfd%egrid(ng))
      cmfd%egrid = mesh_ % energy 
      cmfd % indices(4) = ng - 1 ! sets energy group dimension
    else
      if(.not.allocated(cmfd%egrid)) allocate(cmfd%egrid(2))
      cmfd%egrid = (/0.0_8,20.0_8/)
      cmfd % indices(4) = 1 ! one energy group
    end if

    ! set global albedo
    cmfd % albedo = mesh_ % albedo

    ! get acceleration map
    if (associated(mesh_ % map)) then
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2), &
           cmfd % indices(3)))
      if (size(mesh_ % map) /= product(cmfd % indices(1:3))) then
        message = 'FATAL==>CMFD coremap not to correct dimensions'
        call fatal_error() 
      end if
      cmfd % coremap = reshape(mesh_ % map,(cmfd % indices(1:3)))
      cmfd_coremap = .true.
   end if

    ! check for core map activation by printing note
    if (cmfd_coremap .and. master) then
      message = "Core Map Overlay Activated"
      call write_message()
    end if

    ! check for normalization constant
    cmfd % norm = norm_

    ! set feedback logical
    call lower_case(feedback_)
    if (feedback_ == 'true' .or. feedback_ == '1') cmfd_feedback = .true.

    ! set balance logical
    ! call lower_case(balance_)
    ! if (balance_ == 'true' .or. balance == '1') cmfd_balance = .true.

    ! set downscatter logical
    ! call lower_case(downscatter_)
    ! if (downscatter_ == 'true' .or. downscatter == '1') &
    !      cmfd_downscatter = downscatter_

    ! set 2 group fix
    call lower_case(run_2grp_)
    if (run_2grp_ == 'true' .or. run_2grp_ == '1') cmfd_run_2grp = .true.

    ! set the solver type
    cmfd_solver_type = solver_(1:10)

    ! set monitoring 
    call lower_case(snes_monitor_)
    call lower_case(ksp_monitor_)
    call lower_case(power_monitor_)
    if (snes_monitor_ == 'true' .or. snes_monitor_ == '1') &
         cmfd_snes_monitor = .true.
    if (ksp_monitor_ == 'true' .or. ksp_monitor_ == '1') &
         cmfd_ksp_monitor = .true.
    if (power_monitor_ == 'true' .or. power_monitor_ == '1') &
         cmfd_power_monitor = .true.

    ! output logicals
    call lower_case(write_balance_)
    call lower_case(write_matrices_)
    ! call lower_case(write_hdf5_)
    if (write_balance_ == 'true' .or. write_balance_ == '1') &
         cmfd_write_balance = .true.
    if (write_matrices_ == 'true' .or. write_matrices_ == '1') &
         cmfd_write_matrices = .true.
    ! if (write_hdf5_ == 'true' .or. write_hdf5_ == '1') &
    !      cmfd_write_hdf5 = .true.

    ! run an adjoint calc
    call lower_case(run_adjoint_)
    if (run_adjoint_ == 'true' .or. run_adjoint_ == '1') &
         cmfd_run_adjoint = .true.

    ! batch to begin cmfd
    cmfd_begin = begin_

    ! tally during inactive batches
    call lower_case(inactive_)
    if (inactive_ == 'false' .or. inactive_ == '0') cmfd_tally_on = .false.

    ! inactive batch flush window
    cmfd_inact_flush(1) = inactive_flush_
    cmfd_inact_flush(2) = num_flushes_

    ! last flush before active batches
    cmfd_act_flush = active_flush_

    ! tolerance on keff
    cmfd_keff_tol = keff_tol_
    
    ! create tally objects
    call create_cmfd_tally()

    ! set number of CMFD processors and report to user
    n_procs_cmfd = n_cmfd_procs_ 
    if (master) write(OUTPUT_UNIT,'(A,1X,I0,1X,A)') "CMFD Running on", &
       n_procs_cmfd," processors."

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
    use global
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
    character(MAX_LINE_LEN) :: filename
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter) :: filters(N_FILTER_TYPES) ! temporary filters

    ! parse cmfd.xml file
     filename = trim(path_input) // "cmfd.xml"
     call read_xml_file_cmfd_t(filename)

    ! set global variables if they are 0 (this can happen if there is no tally
    ! file)
    if (n_meshes == 0) n_meshes = n_user_meshes + n_cmfd_meshes

    ! allocate mesh
    if (.not. allocated(meshes)) allocate(meshes(n_meshes))
    m => meshes(n_user_meshes+1)

    ! set mesh id
    m % id = n_user_meshes + 1 

    ! set mesh type to rectangular
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

    ! allocate tallies
    call add_tallies("cmfd", n_cmfd_tallies)

    ! begin loop around tallies
    do i = 1, n_cmfd_tallies

      ! point t to tally variable
      t => cmfd_tallies(i)

      ! set reset property
      call lower_case(reset_)
      if (reset_ == 'true' .or. reset_ == '1') t % reset = .true.

      ! set up mesh filter
      n_filters = 1
      filters(n_filters) % type = FILTER_MESH
      filters(n_filters) % n_bins = product(m % dimension)
      allocate(filters(n_filters) % int_bins(1))
      filters(n_filters) % int_bins(1) = n_user_meshes + 1
      t % find_filter(FILTER_MESH) = n_filters

      ! read and set incoming energy mesh filter
      if (associated(mesh_ % energy)) then
        n_filters = n_filters + 1
        filters(n_filters) % type = FILTER_ENERGYIN
        ng = size(mesh_ % energy)
        filters(n_filters) % n_bins = ng - 1
        allocate(filters(n_filters) % real_bins(ng))
        filters(n_filters) % real_bins = mesh_ % energy
        t % find_filter(FILTER_ENERGYIN) = n_filters
      end if

      ! set number of nucilde bins
      allocate(t % nuclide_bins(1))
      t % nuclide_bins(1) = -1
      t % n_nuclide_bins = 1

      ! record tally id which is equivalent to loop number
      t % id = i_cmfd_tallies + i

      if (i == 1) then

        ! set label
        t % label = "CMFD flux, total, scatter-1, diffusion"

        ! set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! set tally type to volume
        t % type = TALLY_VOLUME

        ! allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! allocate scoring bins 
        allocate(t % score_bins(4))
        t % n_score_bins = 4
        t % n_user_score_bins = 4

        ! allocate scattering order data
        allocate(t % scatt_order(4))
        t % scatt_order = 0
        
        ! set macro_bins
        t % score_bins(1)  = SCORE_FLUX
        t % score_bins(2)  = SCORE_TOTAL
        t % score_bins(3)  = SCORE_SCATTER_N
        t % scatt_order(3) = 1
        t % score_bins(4)  = SCORE_DIFFUSION

      else if (i == 2) then

        ! set label
        t % label = "CMFD neutron production"

        ! set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! set tally type to volume
        t % type = TALLY_VOLUME

        ! read and set outgoing energy mesh filter
        if (associated(mesh_ % energy)) then
          n_filters = n_filters + 1
          filters(n_filters) % type = FILTER_ENERGYOUT
          ng = size(mesh_ % energy)
          filters(n_filters) % n_bins = ng - 1
          allocate(filters(n_filters) % real_bins(ng))
          filters(n_filters) % real_bins = mesh_ % energy
          t % find_filter(FILTER_ENERGYOUT) = n_filters
        end if

        ! allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! deallocate filters bins array
        if (associated(mesh_ % energy)) &
             deallocate(filters(n_filters) % real_bins)

        ! allocate macro reactions
        allocate(t % score_bins(2))
        t % n_score_bins = 2
        t % n_user_score_bins = 2

        ! allocate scattering order data
        allocate(t % scatt_order(2))
        t % scatt_order = 0

        ! set macro_bins
        t % score_bins(1) = SCORE_NU_SCATTER
        t % score_bins(2) = SCORE_NU_FISSION

      else if (i == 3) then

        ! set label
        t % label = "CMFD surface currents"

        ! set tally estimator to analog
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

        ! allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! deallocate filters bins array
        deallocate(filters(n_filters) % int_bins)

        ! allocate macro reactions
        allocate(t % score_bins(1))
        t % n_score_bins = 1
        t % n_user_score_bins = 1

        ! allocate scattering order data
        allocate(t % scatt_order(1))
        t % scatt_order = 0

        ! set macro bins
        t % score_bins(1) = SCORE_CURRENT
        t % type = TALLY_SURFACE_CURRENT

        ! we need to increase the dimension by one since we also need
        ! currents coming into and out of the boundary mesh cells.
        i_filter_mesh = t % find_filter(FILTER_MESH)
        t % filters(i_filter_mesh) % n_bins = product(m % dimension + 1)

      end if

      ! deallocate filter bins
      deallocate(filters(1) % int_bins)
      if (associated(mesh_ % energy)) deallocate(filters(2) % real_bins)

    end do

    call setup_active_cmfdtallies()
    tallies_on = .true.

  end subroutine create_cmfd_tally

end module cmfd_input
