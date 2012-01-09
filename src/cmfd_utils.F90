module cmfd_utils

  use cmfd_header
  use constants
  use datatypes,     only: dict_add_key, dict_get_key
  use error,         only: fatal_error, warning
  use global
  use mesh,          only: mesh_indices_to_bin
  use mesh_header,   only: StructuredMesh
  use string
  use tally_header,  only: TallyObject, TallyScore

  implicit none

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_cmfd_xml()

    use xml_data_cmfd_t

    integer :: ng=1        ! number of energy groups (default 1)
    integer :: n_words     ! number of words read
    logical :: file_exists ! does cmfd.xml exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)

    ! read cmfd infput file
    filename = trim(path_input) // "cmfd.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      write(*,*) "Cannot perform CMFD"
      STOP
    end if

    ! parse cmfd.xml file
    call read_xml_file_cmfd_t(filename)

    ! set spatial dimensions in cmfd object
    cmfd % indices(1:3) = mesh_ % dimension(1:3) ! sets spatial dimensions

    ! get number of energy groups
    if (len_trim(mesh_ % energy) > 0) then
      call split_string(mesh_ % energy, words, n_words)
      ng = n_words - 1
    end if
    cmfd % indices(4) = ng  ! sets energy group dimension

    ! set global albedo
    cmfd % albedo = mesh_ % albedo

    ! get acceleration map
    if (associated(mesh_ % map)) then
      allocate(cmfd % coremap(cmfd % indices(1), cmfd % indices(2),            &
     &         cmfd % indices(3)))
      cmfd % coremap = reshape(mesh_ % map,(cmfd % indices(1:3)))
   end if

    ! create tally objects
    call create_cmfd_tally()

  end subroutine read_cmfd_xml

!===============================================================================
! ALLOCATE_CMFD allocates all of the space for the cmfd object based on tallies
!===============================================================================

  subroutine allocate_cmfd()

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate flux, cross sections and diffusion coefficient
    if (.not. allocated(cmfd % flux)) allocate(cmfd % flux(ng,nx,ny,nz))
    if (.not. allocated(cmfd % totalxs)) allocate(cmfd % totalxs(ng,nx,ny,nz))
    if (.not. allocated(cmfd % p1scattxs)) allocate(cmfd % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(cmfd % scattxs)) allocate(cmfd % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(cmfd % nfissxs)) allocate(cmfd % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(cmfd % diffcof)) allocate(cmfd % diffcof(ng,nx,ny,nz))

    ! allocate dtilde and dhat
    if (.not. allocated(cmfd % dtilde)) allocate(cmfd % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(cmfd % dhat)) allocate(cmfd % dhat(6,ng,nx,ny,nz))

    ! allocate dimensions for each box (here for general case)
    if (.not. allocated(cmfd % hxyz)) allocate(cmfd % hxyz(3,nx,ny,nz))

    ! allocate cmfd fission source pdf
    !allocate( cmfd % sourcepdf(ng,nx,ny,nz) )

    ! allocate surface currents
    if (.not. allocated(cmfd % current)) allocate(cmfd % current(12,ng,nx,ny,nz))

    ! allocate for coremap
    if (.not. allocated(cmfd % coremap)) allocate(cmfd % coremap(nx,ny,nz))

  end subroutine allocate_cmfd

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(g,i,j,k,ng,nx,ny)

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index
    integer :: ng              ! max energy groups
    integer :: nx              ! maximum cells in x direction
    integer :: ny              ! maximum cells in y direction

    ! local variables
    integer :: nidx            ! index in matrix

    ! check if coremap is used
    if (allocated(cmfd % coremap)) then

      ! get idx from core map
      nidx = ng*(cmfd % coremap(i,j,k)) + (ng - g)

    else

      ! compute index
      nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

!===============================================================================
! PRINT_CMFD is a test routine to check if info from tally is being accessed 
!===============================================================================

  subroutine print_cmfd()

    integer :: bins(TALLY_TYPES)       ! bin for tally_types, for filters
    integer :: ijk(3)                  ! indices for mesh cell where tally is
    integer :: score_index             ! index in tally score to get value

    real(8) :: tally_val               ! value of tally being extracted    

    type(TallyObject), pointer :: t    ! pointer for a tally object
    type(StructuredMesh), pointer :: m ! pointer for mesh object

    ! associate pointers with objects
    t => tallies(3)
    m => meshes(t % mesh)

    ! set all bins to 1
    bins = 1

    ! get mesh indices, first we will first force to 1,1,1
!   ijk = (/ 1, 1, 1 /)

    ! apply filters, here we will just try a mesh filter first
 !  bins(T_MESH) = mesh_indices_to_bin(m,ijk)

    ! calculate score index from bins
 !  score_index = sum((bins - 1) * t%stride) + 1

    ! get value from tally object
 !  tally_val = t%scores(score_index,2)%val

    ! write value to file
 !  write(7,*) "Tally value is:",tally_val

    ! Left Surface
    ijk = (/ 1-1, 1, 1 /)
    score_index = sum(t % stride(1:3) * ijk) + IN_RIGHT
    print *, "Outgiong Current from Left", t % scores(score_index,1) % val
    score_index = sum(t % stride(1:3) * ijk) + OUT_RIGHT
    print *, "Incoming Current from Left", t % scores(score_index,1) % val


  end subroutine print_cmfd

!===============================================================================
! CREATE_CMFD_TALLY creates the tally object for OpenMC to process for CMFD
! accleration.
! There are 3 tally types:
!   1: Only an energy in filter-> flux,total,p1 scatter
!   2: Energy in and energy out filter-> nu-scatter,nu-fission
!   3: Surface current
!===============================================================================

  subroutine create_cmfd_tally()

    use xml_data_cmfd_t 

    integer :: i           ! loop counter
    integer :: j           ! loop counter
    integer :: id          ! user-specified identifier
    integer :: index       ! index in mesh array
    integer :: n           ! size of arrays in mesh specification
    integer :: ng=1        ! number of energy groups (default 1)
    integer :: n_words     ! number of words read
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    ! parse cmfd.xml file
     filename = trim(path_input) // "cmfd.xml"
     call read_xml_file_cmfd_t(filename)

    ! allocate mesh
    n_meshes = 1
    allocate(meshes(n_meshes))
    m => meshes(1)

    ! set mesh id
    m % id = 1

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
    allocate(m % origin(n))
    allocate(m % width(n))
    allocate(m % upper_right(n))

    ! read dimensions in each direction
    m % dimension = mesh_ % dimension

    ! read mesh origin location
    if (m % n_dimension /= size(mesh_ % origin)) then
      message = "Number of entries on <origin> must be the same as " // &
                "the number of entries on <dimension>."
      call fatal_error()
    end if
    m % origin = mesh_ % origin

    ! read mesh widths
    if (size(mesh_ % width) /= size(mesh_ % origin)) then
       message = "Number of entries on <width> must be the same as " // &
                 "the number of entries on <origin>."
       call fatal_error()
    end if
    m % width = mesh_ % width

    ! set upper right coordinate
    m % upper_right = m % origin + m % dimension * m % width

    ! add mesh to dictionary
    call dict_add_key(mesh_dict, m % id, 1)

    ! allocate tallies
    n_tallies = 3
    allocate(tallies(n_tallies))

    ! begin loop around tallies
    do i = 1,n_tallies
      t => tallies(i)

      ! allocate arrays for number of bins and stride in scores array
      allocate(t % n_bins(TALLY_TYPES))
      allocate(t % stride(TALLY_TYPES))

      ! initialize number of bins and stride
      t % n_bins = 0
      t % stride = 0

      ! record tally id which is equivalent to loop number
      t % id = i

      ! set mesh filter mesh id = 1
      t % mesh = 1
      m => meshes(1)
      t % n_bins(T_MESH) = t % n_bins(T_MESH) + product(m % dimension)

      ! read and set incoming energy mesh filter
      if (len_trim(mesh_ % energy) > 0) then
        call split_string(mesh_ % energy,words,n_words)
        ng = n_words
        allocate(t % energy_in(n_words))
        do j = 1,n_words
          t % energy_in(j) = str_to_real(words(j))
        end do
        t % n_bins(T_ENERGYIN) = n_words - 1
      end if

      if (i == 1) then

        ! allocate macro reactions
        allocate(t % macro_bins(3))
        t % n_macro_bins = 3 

        ! set macro_bins
        t % macro_bins(1) % scalar = MACRO_FLUX
        t % macro_bins(2) % scalar = MACRO_TOTAL
        t % macro_bins(3) % scalar = MACRO_SCATTER_1

      else if (i == 2) then

        ! read and set outgoing energy mesh filter
        if (len_trim(mesh_ % energy) > 0) then
          call split_string(mesh_ % energy, words, n_words)
          allocate(t % energy_out(n_words))
          do j = 1, n_words
            t % energy_out(j) = str_to_real(words(j))
          end do
          t % n_bins(T_ENERGYOUT) = n_words - 1
        end if

        ! allocate macro reactions
        allocate(t % macro_bins(2))
        t % n_macro_bins = 2

        ! set macro_bins
        t % macro_bins(1) % scalar = MACRO_NU_SCATTER
        t % macro_bins(2) % scalar = MACRO_NU_FISSION

      else if (i == 3) then

        ! allocate macro reactions
        allocate(t % macro_bins(1))
        t % n_macro_bins = 1

        ! set macro bins
        t % macro_bins(1) % scalar = MACRO_CURRENT
        t % surface_current = .true.

        ! since the number of bins for the mesh filter was already set
        ! assuming it was a flux tally, we need to adjust the number of
        ! bins
        t % n_bins(T_MESH) = t % n_bins(T_MESH) - product(m % dimension)

        ! get pointer to mesh
        id = t % mesh
        index = dict_get_key(mesh_dict, id)
        m => meshes(index)

        ! we need to increase the dimension by one since we also need
        ! currents coming into and out of the boundary mesh cells.
        if (size(m % dimension) == 2) then
          t % n_bins(T_MESH) = t % n_bins(T_MESH) + &
       &                       product(m % dimension + 1) * 4
        elseif (size(m % dimension) == 3) then
          t % n_bins(T_MESH) = t % n_bins(T_MESH) + &
                               product(m % dimension + 1) * 6
        end if

      end if

    end do

  end subroutine create_cmfd_tally

!===============================================================================
! NEUTRON_BALANCE writes a file that contains n. bal. info for all cmfd mesh 
!===============================================================================

  subroutine neutron_balance()

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z
    integer :: g                 ! iteration counter for g
    integer :: h                 ! iteration counter for outgoing groups
    integer :: l                 ! iteration counter for leakage
    integer :: io_error          ! error for opening file unit
    real(8) :: leakage           ! leakage term in neutron balance
    real(8) :: interactions      ! total number of interactions in balance
    real(8) :: scattering        ! scattering term in neutron balance
    real(8) :: fission           ! fission term in neutron balance
    real(8) :: res               ! residual of neutron balance
    character(MAX_FILE_LEN) :: filename
    character(30)           :: label

    ! open cmfd file for output
    filename = "cmfd.out"
    open(FILE=filename, UNIT=UNIT_CMFD, STATUS='replace', ACTION='write',      &
         IOSTAT=io_error)

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! begin loop around space and energy groups
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUPG: do g = 1,ng

            ! get leakage
            leakage = 0.0
            LEAK: do l = 1,3

              leakage = leakage + ((cmfd % current(4*l,g,i,j,k) -              &
             &          cmfd % current(4*l-1,g,i,j,k))) -                      &
             &         ((cmfd % current(4*l-2,g,i,j,k) -                       &
             &          cmfd % current(4*l-3,g,i,j,k)))

            end do LEAK

            ! interactions
            interactions = cmfd % totalxs(g,i,j,k) * cmfd % flux(g,i,j,k)

            ! get scattering and fission
            scattering = 0.0
            fission = 0.0
            GROUPH: do h = 1,ng

              scattering = scattering + cmfd % scattxs(h,g,i,j,k) *            &
             &                          cmfd % flux(h,i,j,k)

              fission = fission + cmfd % nfissxs(h,g,i,j,k) *                  &
             &                    cmfd % flux(h,i,j,k)

            end do GROUPH

            ! compute residual
            res = leakage + interactions - scattering - (ONE/keff)*fission

            ! write output
            label = "MESH (" // trim(int4_to_str(i)) // ". " //                &
           &        trim(int4_to_str(j)) // ", " // trim(int4_to_str(k)) //    &
           &        ")  GROUP " // trim(int4_to_str(g))
            write(UNIT=UNIT_CMFD, FMT='(A,T35,A)') label,                      &
           &      trim(real_to_str(res))

          end do GROUPG

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! close file
    close(UNIT=UNIT_CMFD)

    STOP

  end subroutine neutron_balance

!===============================================================================
! SET_COREMAP is a routine that sets the core mapping information
!===============================================================================

  subroutine set_coremap()

    integer :: kount=1           ! counter for unique fuel assemblies
    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z

    ! extract spatial indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)

    ! count how many fuel assemblies exist
    cmfd % mat_dim = sum(cmfd % coremap - 1)

    ! begin loops over spatial indices
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check for reflector
          if (cmfd % coremap(i,j,k) == 1) then

            ! reset value to 99999
            cmfd % coremap(i,j,k) = 99999

          else

            ! must be a fuel --> give unique id number
            cmfd % coremap(i,j,k) = kount
            kount = kount + 1

          end if

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine set_coremap

!===============================================================================
! GET_REFLECTOR_ALBEDO is a function that calculates the albedo to the reflector
!===============================================================================

  function get_reflector_albedo(l,g,i,j,k)

    ! function variable
    real(8) :: get_reflector_albedo ! reflector albedo

    ! local variable
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: l                    ! iteration counter for leakages
    integer :: shift_idx            ! parameter to shift index by +1 or -1
    real(8) :: current(12)          ! partial currents for all faces of mesh cell            
    real(8) :: albedo               ! the albedo

    ! get partial currents from object
    current = cmfd%current(:,g,i,j,k)

    ! define xyz and +/- indices
    shift_idx = -2*mod(l,2) + 1          ! shift neig by -1 or +1

    ! calculate albedo
    albedo = (current(2*l-1)/current(2*l))**(shift_idx)

    ! assign to function variable
    get_reflector_albedo = albedo

  end function get_reflector_albedo

!===============================================================================
! WRITE_HDF5 writes an hdf5 output file with the cmfd object for restarts 
!===============================================================================

  subroutine write_hdf5()

    use hdf5 

    character(LEN=7), parameter :: filename = "cmfd.h5"  ! File name
    character(LEN=4), parameter :: grpname = "cmfd"       ! Group name

    integer(HID_T) :: file_id             ! File identifier
    integer(HID_T) :: group_id            ! Group identifier
    integer(HID_T) :: dataspace_id        ! Data space identifier
    integer(HID_T) :: dataset_id          ! Dataset identifier
    integer        :: error               ! Error flag

    integer(HSIZE_T), dimension(1) :: dim1 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(3) :: dim3 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(4) :: dim4 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(5) :: dim5 ! vector for hdf5 dimensions

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! initialize FORTRAN interface.
    call h5open_f(error)

    ! create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    ! create the CMFD group
    call h5gcreate_f(file_id, grpname, group_id, error)

    ! write indices from cmfd object
    dim1 = (/4/)
    call h5screate_simple_f(1,dim1,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/indices",H5T_NATIVE_INTEGER,dataspace_id,   &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%indices,dim1,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write totalxs from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/totalxs",H5T_NATIVE_DOUBLE,dataspace_id,    &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%totalxs,dim4,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write p1scattxs from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/p1scattxs",H5T_NATIVE_DOUBLE,dataspace_id,  &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%p1scattxs,dim4,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write scattxs from cmfd object
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/scattxs",H5T_NATIVE_DOUBLE,dataspace_id,    &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%scattxs,dim5,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write nfissxs from cmfd object
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/nfissxs",H5T_NATIVE_DOUBLE,dataspace_id,    &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%nfissxs,dim5,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write diffcof from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/diffcof",H5T_NATIVE_DOUBLE,dataspace_id,    &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%diffcof,dim4,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write current from cmfd object
    dim5 = (/12,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/current",H5T_NATIVE_DOUBLE,dataspace_id,    &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%current,dim5,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write flux from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/flux",H5T_NATIVE_DOUBLE,dataspace_id,       &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%flux,dim4,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write dtilde from cmfd object
    dim5 = (/6,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/dtilde",H5T_NATIVE_DOUBLE,dataspace_id,     &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dtilde,dim5,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write dhat from cmfd object
    dim5 = (/6,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/dhat",H5T_NATIVE_DOUBLE,dataspace_id,       &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dhat,dim5,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write albedo from cmfd object
    dim1 = (/6/)
    call h5screate_simple_f(1,dim1,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/albedo",H5T_NATIVE_DOUBLE,dataspace_id,     &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%albedo,dim1,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write hxyz from cmfd object
    dim4 = (/3,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/hxyz",H5T_NATIVE_DOUBLE,dataspace_id,       &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%hxyz,dim4,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write coremap from cmfd object
    dim3 = (/nx,ny,nz/)
    call h5screate_simple_f(3,dim3,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/coremap",H5T_NATIVE_INTEGER,dataspace_id,   &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%coremap,dim3,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! write mat_dim from cmfd object
    dim1 = (/1/)
    call h5screate_simple_f(1,dim1,dataspace_id,error)
    call h5dcreate_f(file_id,"cmfd/mat_dim",H5T_NATIVE_INTEGER,dataspace_id,   &
   &                 dataset_id,error)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%mat_dim,dim1,error)
    call h5sclose_f(dataspace_id,error)
    call h5dclose_f(dataset_id,error)

    ! close the CMFD group
    call h5gclose_f(group_id, error)

    ! terminate access to the file.
    call h5fclose_f(file_id, error)

    ! close FORTRAN interface.
    call h5close_f(error)

  end subroutine write_hdf5

!===============================================================================
! READ_HDF5 writes an hdf5 output file with the cmfd object for restarts 
!===============================================================================

  subroutine read_hdf5()

    use hdf5

    character(LEN=8), parameter :: filename = "cmfd.h5"  ! File name

    integer(HID_T) :: file_id             ! File identifier
    integer(HID_T) :: dataset_id          ! Dataset identifier
    integer        :: error               ! Error flag

    integer(HSIZE_T), dimension(1) :: dim1
    integer(HSIZE_T), dimension(3) :: dim3
    integer(HSIZE_T), dimension(4) :: dim4 
    integer(HSIZE_T), dimension(5) :: dim5

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups

    ! initialize the Fortran interface
    call h5open_f(error)

    ! open the restart file
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)

    ! read indices to cmfd object 
    call h5dopen_f(file_id,"cmfd/indices",dataset_id,error)
    dim1 = (/4/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%indices,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! get indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate cmfd object
    call allocate_cmfd()

    ! read totalxs to cmfd object
    call h5dopen_f(file_id,"cmfd/totalxs",dataset_id,error)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%totalxs,dim4,error)
    call h5dclose_f(dataset_id,error)

    ! read p1scattxs to cmfd object
    call h5dopen_f(file_id,"cmfd/p1scattxs",dataset_id,error)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%p1scattxs,dim4,error)
    call h5dclose_f(dataset_id,error)

    ! read scattxs to cmfd object
    call h5dopen_f(file_id,"cmfd/scattxs",dataset_id,error)
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%scattxs,dim5,error)
    call h5dclose_f(dataset_id,error)

    ! read scattxs to cmfd object
    call h5dopen_f(file_id,"cmfd/nfissxs",dataset_id,error)
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%nfissxs,dim5,error)
    call h5dclose_f(dataset_id,error)

    ! read diffcof to cmfd object
    call h5dopen_f(file_id,"cmfd/diffcof",dataset_id,error)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%diffcof,dim4,error)
    call h5dclose_f(dataset_id,error)

    ! read current to cmfd object
    call h5dopen_f(file_id,"cmfd/current",dataset_id,error)
    dim5 = (/12,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%current,dim5,error)
    call h5dclose_f(dataset_id,error)

    ! read flux to cmfd object
    call h5dopen_f(file_id,"cmfd/flux",dataset_id,error)
    dim4 = (/ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%flux,dim4,error)
    call h5dclose_f(dataset_id,error)

    ! read dtilde to cmfd object
    call h5dopen_f(file_id,"cmfd/dtilde",dataset_id,error)
    dim5 = (/6,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dtilde,dim5,error)
    call h5dclose_f(dataset_id,error)

    ! read dhat to cmfd object
    call h5dopen_f(file_id,"cmfd/dhat",dataset_id,error)
    dim5 = (/6,ng,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dhat,dim5,error)
    call h5dclose_f(dataset_id,error)

    ! read albedo to cmfd object
    call h5dopen_f(file_id,"cmfd/albedo",dataset_id,error)
    dim1 = (/6/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%albedo,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! read hxyz to cmfd object
    call h5dopen_f(file_id,"cmfd/hxyz",dataset_id,error)
    dim4 = (/3,nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%hxyz,dim4,error)
    call h5dclose_f(dataset_id,error)

    ! read coremap to cmfd object
    call h5dopen_f(file_id,"cmfd/coremap",dataset_id,error)
    dim3 = (/nx,ny,nz/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%coremap,dim3,error)
    call h5dclose_f(dataset_id,error)

    ! read mat_dim to cmfd object
    call h5dopen_f(file_id,"cmfd/mat_dim",dataset_id,error)
    dim1 = (/1/)
    call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%mat_dim,dim1,error)
    call h5dclose_f(dataset_id,error)

    ! Close the file
    call h5fclose_f(file_id,error)

    ! Close the Fortran interface
    call h5close_f(error)

  end subroutine read_hdf5

end module cmfd_utils
