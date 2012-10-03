module cmfd_input

  implicit none
  private
  public :: configure_cmfd 

contains

!===============================================================================
! CONFIGURE_CMFD
!===============================================================================

  subroutine configure_cmfd()

    use cmfd_message_passing,   only: petsc_init_mpi
    use global,  only: cmfd, cmfd_write_hdf5, master

#ifdef HDF5
    use hdf5
    integer(Fortran_Integer) :: hdf5_err
#endif

    ! read in cmfd input file
    call read_cmfd_xml()

    ! write out summary to standard out
!   call write_cmfd_input_summary(cmfd_tally_size)

    ! initialize petsc on mpi
    call petsc_init_mpi()

    ! Create a new file using default properties.
# ifdef HDF5
    if (cmfd_write_hdf5 .and. master)                              &
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F,cmfd%file_id, hdf5_err)
# endif

  end subroutine configure_cmfd

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_cmfd_xml()
    
    use error
    use global
    use output
    use string
    use xml_data_cmfd_t

    integer :: j           ! iteration counter
    integer :: ng          ! number of energy groups
    integer :: n_words     ! number of words read
    logical :: file_exists ! does cmfd.xml exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)

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
      cmfd % indices(4) = ng  ! sets energy group dimension
    else
      if(.not.allocated(cmfd%egrid)) allocate(cmfd%egrid(2))
      cmfd%egrid = (/0.0_8,20.0_8/)
      cmfd % indices(4) = 1 ! one energy group
    end if

    ! set global albedo
    cmfd % albedo = mesh_ % albedo

    ! get acceleration map
    if (associated(mesh_ % map)) then
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2),            &
     &         cmfd % indices(3)))
      if (size(mesh_ % map) /= product(cmfd % indices(1:3))) then
        write(*,*) 'FATAL==>CMFD coremap not to correct dimensions'
        stop
      end if
      cmfd % coremap = reshape(mesh_ % map,(cmfd % indices(1:3)))
      cmfd_coremap = .TRUE.
   end if

    ! check for core map activation by printing note
    if (cmfd_coremap .and. master) write(*,*)"Core Map Overlay Activated"

    ! check for normalization constant
    cmfd % norm = norm_

    ! set feedback logical
    cmfd_feedback = feedback_

    ! set balance logical
    ! cmfd_balance = balance_

    ! set downscatter logical
    ! cmfd_downscatter = downscatter_

    ! set 2 group fix
    cmfd_run_2grp = run_2grp_

    ! set the solver type
    cmfd_solver_type = solver_

    ! set monitoring 
    cmfd_snes_monitor = snes_monitor_
    cmfd_ksp_monitor = ksp_monitor_
    cmfd_power_monitor = power_monitor_

    ! output logicals
    cmfd_write_balance = write_balance_
    cmfd_write_matrices = write_matrices_
    cmfd_write_hdf5 = write_hdf5_

    ! run an adjoint calc
    cmfd_run_adjoint = run_adjoint_

    ! batch to begin cmfd
    cmfd_begin = begin_

    ! tally during inactive batches
    cmfd_tally_on = inactive_

    ! inactive batch flush window
!   cmfd_inact_flush = inactive_flush_

    ! last flush before active batches
    cmfd_act_flush = active_flush_

    ! tolerance on keff
    cmfd_keff_tol = keff_tol_
    
    ! create tally objects
    call create_cmfd_tally()

    ! set number of CMFD processors and report to user
    n_procs_cmfd = n_cmfd_procs_ 
    if (master) write(*,*) "CMFD Running on",n_procs_cmfd," processors."

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

    use datatypes,     only: dict_add_key, dict_get_key
    use error,         only: fatal_error, warning
    use global
    use mesh_header,   only: StructuredMesh
    use string
    use tally_header,  only: TallyObject, TallyScore
    use xml_data_cmfd_t

    integer :: i           ! loop counter
    integer :: j           ! loop counter
    integer :: id          ! user-specified identifier
    integer :: i_mesh      ! index in mesh array
    integer :: n           ! size of arrays in mesh specification
    integer :: ng          ! number of energy groups (default 1)
    integer :: n_filters   ! number of filters
    integer :: filters(N_FILTER_TYPES) ! temp list of filters
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    ! parse cmfd.xml file
     filename = trim(path_input) // "cmfd.xml"
     call read_xml_file_cmfd_t(filename)

    ! set global variables
    n_meshes = n_user_meshes + n_cmfd_meshes
    n_tallies = n_user_tallies + n_cmfd_tallies
    n_analog_tallies = n_user_analog_tallies + n_cmfd_analog_tallies
    n_tracklength_tallies = n_user_tracklength_tallies + n_cmfd_tracklength_tallies
    n_current_tallies = n_user_current_tallies + n_cmfd_current_tallies

    ! Allocate list of pointers for tallies by type
    if (.not. allocated(analog_tallies))      allocate(analog_tallies(n_analog_tallies))
    if (.not. allocated(tracklength_tallies)) allocate(tracklength_tallies(n_tracklength_tallies))
    if (.not. allocated(current_tallies))     allocate(current_tallies(n_current_tallies))

    ! allocate mesh
    if (.not. allocated(meshes)) allocate(meshes(n_meshes))
    m => meshes(n_user_meshes+1)

    ! set mesh id
    m % id = n_user_meshes + 1 

    ! set mesh type to rectangular
    m % type = LATTICE_RECT

    ! determine number of dimensions for mesh
    n = size(mesh_ % dimension)
    if (n /= 2 .and. n /= 3) then
      message = "Mesh must be two or three dimensions."
      call fatal_error()
    end if
    m % n_dimension = n

    ! allocate attribute arrays
    allocate(m % dimension(n))
    allocate(m % lower_left(n))
    allocate(m % width(n))
    allocate(m % upper_right(n))

    ! read dimensions in each direction
    m % dimension = mesh_ % dimension

    ! read mesh lower left location
    if (m % n_dimension /= size(mesh_ % lower_left)) then
      message = "Number of entries on <lower_left> must be the same as " // &
                "the number of entries on <dimension>."
      call fatal_error()
    end if
    m % lower_left = mesh_ % lower_left

    ! read mesh widths
    if (size(mesh_ % width) /= size(mesh_ % lower_left)) then
       message = "Number of entries on <width> must be the same as " // &
                 "the number of entries on <lower_left>."
       call fatal_error()
    end if
    m % width = mesh_ % width

    ! set upper right coordinate
    m % upper_right = m % lower_left + m % dimension * m % width

    ! add mesh to dictionary
    call dict_add_key(mesh_dict, m % id, n_user_meshes + 1)

    ! allocate tallies
    if (.not. allocated(tallies)) allocate(tallies(n_tallies))

    ! begin loop around tallies
    do i = n_user_tallies+1,n_tallies

      ! set n filters to 0
      n_filters = 0
      filters = 0

      ! point t to tally variable
      t => tallies(i)

      ! set reset property
      if (reset_) t % reset = .true.

      ! allocate arrays for number of bins and stride in scores array
      allocate(t % n_filter_bins(N_FILTER_TYPES))
      allocate(t % stride(N_FILTER_TYPES))

      ! initialize number of bins and stride
      t % n_filter_bins = 0
      t % stride = 0

      ! record tally id which is equivalent to loop number
      t % id = i

      ! set mesh filter mesh id
      t % mesh = n_user_meshes + 1
      t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) +        &
     &                                 product(m % dimension)
      n_filters = n_filters + 1
      filters(n_filters) = FILTER_MESH

      ! read and set incoming energy mesh filter
      if (associated(mesh_ % energy)) then
        ng = size(mesh_ % energy)
        allocate(t % energy_in(ng))
        t % energy_in = mesh_ % energy 
        t % n_filter_bins(FILTER_ENERGYIN) = ng - 1
        n_filters = n_filters + 1
        filters(n_filters) = FILTER_ENERGYIN
      end if

      if (i == n_user_tallies+1) then

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

        ! set macro_bins
        t % score_bins(1) = SCORE_FLUX
        t % score_bins(2) = SCORE_TOTAL
        t % score_bins(3) = SCORE_SCATTER_1
        t % score_bins(4) = SCORE_DIFFUSION

        ! Increment the appropriate index and set pointer
        analog_tallies(n_user_analog_tallies + 1) = i

      else if (i == n_user_tallies + 2) then

        ! set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! set tally type to volume
        t % type = TALLY_VOLUME

        ! read and set outgoing energy mesh filter
        if (associated(mesh_ % energy)) then
          ng = size(mesh_ % energy)
          allocate(t % energy_out(ng))
          t % energy_out = mesh_ % energy 
          t % n_filter_bins(FILTER_ENERGYOUT) = ng - 1
          n_filters = n_filters + 1
          filters(n_filters) = FILTER_ENERGYOUT
        end if

        ! allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! allocate macro reactions
        allocate(t % score_bins(2))
        t % n_score_bins = 2

        ! set macro_bins
        t % score_bins(1) = SCORE_NU_SCATTER
        t % score_bins(2) = SCORE_NU_FISSION

        ! Increment the appropriate index and set pointer
        analog_tallies(n_user_analog_tallies + 2) = i

      else if (i == n_user_tallies + 3) then

        ! set tally estimator to analog
        t % estimator = ESTIMATOR_ANALOG

        ! allocate and set filters
        t % n_filters = n_filters
        allocate(t % filters(n_filters))
        t % filters = filters(1:n_filters)

        ! allocate macro reactions
        allocate(t % score_bins(1))
        t % n_score_bins = 1

        ! set macro bins
        t % score_bins(1) = SCORE_CURRENT
        t % type = TALLY_SURFACE_CURRENT

        ! since the number of bins for the mesh filter was already set
        ! assuming it was a flux tally, we need to adjust the number of
        ! bins
        t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH)        &
       &                                 - product(m % dimension)

        ! get pointer to mesh
        id = t % mesh
        i_mesh = dict_get_key(mesh_dict, id)
        m => meshes(i_mesh)

        ! we need to increase the dimension by one since we also need
        ! currents coming into and out of the boundary mesh cells.
        if (size(m % dimension) == 2) then  ! these lines need changing
          t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) +    &
         &                                 product(m % dimension + 1) * 4
        elseif (size(m % dimension) == 3) then
          t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) +    &
         &                                 product(m % dimension + 1) * 6
        end if

        ! Increment the appropriate index and set pointer
        current_tallies(n_user_current_tallies + 1) = i 

      end if

    end do

  end subroutine create_cmfd_tally

!===============================================================================
! READ_CMFD_HDF5 writes an hdf5 output file with the cmfd object for restarts
!===============================================================================

  subroutine read_cmfd_hdf5()

    use cmfd_header, only: allocate_cmfd
    use global,      only: cmfd,cmfd_coremap

#ifdef HDF5
    use global, only: hdf5_output_file,hdf5_err
    use hdf5
    use hdf5_interface, only: hdf5_open_output, hdf5_close_output

    integer(HID_T) :: dataset_id ! Dataset identifier
    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4
    integer(HSIZE_T), dimension(5) :: dim5
    integer :: nx ! number of mesh cells in x direction
    integer :: ny ! number of mesh cells in y direction
    integer :: nz ! number of mesh cells in z direction
    integer :: ng ! number of energy groups
    integer :: core_map_int 

    ! open output file
    call hdf5_open_output()

    ! read indices to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/indices",dataset_id,hdf5_err)
    dim1 = (/4/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%indices,dim1,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! get indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate cmfd object
    call allocate_cmfd(cmfd)

    ! read totalxs to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/totalxs",dataset_id,hdf5_err)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%totalxs,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read p1scattxs to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/p1scattxs",dataset_id,hdf5_err)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%p1scattxs,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read scattxs to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/scattxs",dataset_id,hdf5_err)
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%scattxs,dim5,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read scattxs to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/nfissxs",dataset_id,hdf5_err)
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%nfissxs,dim5,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read diffcof to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/diffcof",dataset_id,hdf5_err)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%diffcof,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read current to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/current",dataset_id,hdf5_err)
    dim5 = (/12,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%current,dim5,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read flux to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/flux",dataset_id,hdf5_err)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%flux,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read dtilde to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/dtilde",dataset_id,hdf5_err)
    dim5 = (/6,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dtilde,dim5,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read dhat to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/dhat",dataset_id,hdf5_err)
    dim5 = (/6,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dhat,dim5,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read albedo to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/albedo",dataset_id,hdf5_err)
    dim1 = (/6/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%albedo,dim1,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read hxyz to cmfd object
    call h5dopen_f(hdf5_output_file,"cmfd/hxyz",dataset_id,hdf5_err)
    dim4 = (/3,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%hxyz,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! read in core_map logical 
    call h5dopen_f(hdf5_output_file,"cmfd/coremap_active",dataset_id,hdf5_err)
    dim1 = (/1/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,core_map_int,dim1,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! now set coremap depending on logical
    select case(core_map_int)

      ! core map is not active
      case(0)

        ! set logical to false
        cmfd_coremap = .FALSE.

      ! core map is active
      case(1)

        ! set logical to true
        cmfd_coremap = .TRUE.

        ! allocate coremap in cmfd obj
        allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2),          &
       &         cmfd % indices(3)))


        ! read coremap to cmfd object
        call h5dopen_f(hdf5_output_file,"cmfd/coremap",dataset_id,hdf5_err)
        dim3 = (/nx,ny,nz/)
        call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%coremap,dim3,hdf5_err)
        call h5dclose_f(dataset_id,hdf5_err)

        ! read mat_dim to cmfd object
        call h5dopen_f(hdf5_output_file,"cmfd/mat_dim",dataset_id,hdf5_err)
        dim1 = (/1/)
        call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%mat_dim,dim1,hdf5_err)
        call h5dclose_f(dataset_id,hdf5_err)

      ! something is wrong
      case default

        write(*,*) 'FATAL ==> Could not detect core map in hdf5 output file'
        stop

      end select

    ! close output file
    call hdf5_close_output()

#endif

  end subroutine read_cmfd_hdf5

end module cmfd_input
