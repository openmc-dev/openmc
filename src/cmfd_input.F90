module cmfd_input

  implicit none
  private
  public :: read_cmfd_xml

contains

!===============================================================================
! READ_INPUT reads the CMFD input file and organizes it into a data structure
!===============================================================================

  subroutine read_cmfd_xml()

    use global
    use string
    use xml_data_cmfd_t

    integer :: ng=1        ! number of energy groups (default 1)
    integer :: n_words     ! number of words read
    logical :: file_exists ! does cmfd.xml exist?
    character(MAX_LINE_LEN) :: filename
    character(MAX_WORD_LEN) :: words(MAX_WORDS)

    ! read cmfd infput file
    filename = "cmfd.xml"
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

    ! check for core map activation by printing note
    if (allocated(cmfd % coremap)) print *,"Core Map Overlay Activated"

    ! create tally objects
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

end module cmfd_input
