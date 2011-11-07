module cmfd_utils

  use datatypes,     only: dict_add_key, dict_get_key
  use error,         only: fatal_error, warning
  use global
  use mesh,          only: mesh_indices_to_bin
  use mesh_header,   only: StructuredMesh
  use string,        only: lower_case, str_to_real, split_string
  use tally_header,  only: TallyObject, TallyScore

  implicit none

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_input()

  end subroutine read_input

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

    ! compute index
    nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

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
    t => tallies(1)
    m => meshes(t % mesh)

    ! set all bins to 1
    bins = 1

    ! get mesh indices, first we will first force to 1,1,1
    ijk = (/ 1, 1, 1 /)

    ! apply filters, here we will just try a mesh filter first
    bins(T_MESH) = mesh_indices_to_bin(m,ijk)

    ! calculate score index from bins
    score_index = sum((bins - 1) * t%stride) + 1

    ! get value from tally object
    tally_val = t%scores(score_index,1)%val

    ! write value to file
    write(7,*) "Tally value is:",tally_val

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
    integer :: n_words     ! number of words read
    logical :: file_exists ! does cmfd.xml file exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)
    type(TallyObject),    pointer :: t => null()
    type(StructuredMesh), pointer :: m => null()

    ! check if cmfd.xml exists
    filename = trim(path_input) // "cmfd.xml"
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
      write(*,*) "Cannot perform CMFD"
      STOP
    end if

    ! parse cmfd.xml file
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
      if (len_trim(energy_) > 0) then
        call split_string(energy_,words,n_words)
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
        if (len_trim(energy_) > 0) then
          call split_string(energy_, words, n_words)
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

end module cmfd_utils
