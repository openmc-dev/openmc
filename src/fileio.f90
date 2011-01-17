module fileio

  use global
  use types,  only: Cell, Surface, CellList, SurfaceList
  use string, only: split_string_wl, lower_case, string_to_real, split_string
  use output, only: message, warning, error

  implicit none

contains

!------------------------------------------------------------------------------

  subroutine read_command_line()

    integer :: argc
    character(250) :: msg
    logical :: file_exists
    character(7) :: readable

    argc = COMMAND_ARGUMENT_COUNT()
    if ( argc > 0 ) then
       call GET_COMMAND_ARGUMENT(1,inputfile)
    else
       msg = "No input file specified!"
       call error( msg )
    end if

    ! Check if input file exists
    inquire( FILE=inputfile, EXIST=file_exists, READ=readable )
    if ( .not. file_exists ) then
       msg = "Input file '" // trim(inputfile) // "' does not exist!"
       call error( msg )
    elseif ( readable(1:3) /= 'YES' ) then
       msg = "Input file '" // trim(inputfile) // "' is not readable! &
            &Change file permissions with chmod command."
       call error( msg )
    end if

  end subroutine read_command_line

!------------------------------------------------------------------------------

  subroutine read_input(filename)

    character(*), intent(in) :: filename

    character(250) :: line
    character(32) :: words(max_words)
    integer :: in = 7
    integer :: ioError
    integer :: n, i

    type(CellList), pointer :: c_head => null()
    type(CellList), pointer :: c_tail => null()
    type(CellList), pointer :: c_current => null()
    
    type(SurfaceList), pointer :: s_head => null()
    type(SurfaceList), pointer :: s_tail => null()
    type(SurfaceList), pointer :: s_current => null()

!!$    type(MaterialList), pointer :: m_head
!!$    type(MaterialList), pointer :: m_tail
!!$    type(MaterialList), pointer :: m_current

    ncell = 0
    nsurf = 0
    nmat  = 0


    open(file=filename, unit=in, status='old', action='read')

    do
       read(unit=in, fmt='(A250)', iostat=ioError) line
       if ( ioError /= 0 ) exit
       call split_string_wl(line, words, n)

       select case ( trim(words(1)) )
       case ('cell')
          ! Allocate a new cell
          if ( .not. associated(c_head) ) then
             allocate( c_head )
             c_tail => c_head
          else
             allocate( c_tail%next )
             c_tail => c_tail%next
          end if

          ! Parse data
          call parse_cell( c_tail, words, n )
          ncell = ncell + 1
          
       case ('surface')
          ! Allocate a new surface
          if ( .not. associated(s_head) ) then
             allocate( s_head )
             s_tail => s_head
          else
             allocate( s_tail%next )
             s_tail => s_tail%next
          end if

          ! Parse data
          call parse_surface( s_tail, words, n )
          nsurf = nsurf + 1

       case ( 'material' )
          nmat = nmat + 1
       case default
          ! do nothing
       end select

    end do

    close(unit=in)

    ! Convert linked list of Cells into normal array
    allocate( cells(ncell) )
    c_current => c_head
    do i = 1, ncell
       ! Allocate space on cells(i)
       n = size(c_current%boundary_list)
       allocate( cells(i)%boundary_list(n) )

       ! Copy cell
       cells(i)%id = c_current%id
       cells(i)%material = c_current%material
       cells(i)%boundary_list = c_current%boundary_list

!!$       print *, 'Cell: ', cells(i)%id
!!$       print *, '   Material ', cells(i)%material
!!$       print *, '   Surfaces ', cells(i)%boundary_list
       
       c_head => c_head%next
       deallocate( c_current )
       c_current => c_head
    end do

    ! Convert linked list of Surfaces into normal array
    allocate( surfaces(nsurf) )
    s_current => s_head
    do i = 1, nsurf
       ! Allocate space on surfaces(i)
       n = size(s_current%coeffs)
       allocate( surfaces(i)%coeffs(n) )

       ! Copy surface
       surfaces(i)%id = s_current%id
       surfaces(i)%type = s_current%type
       surfaces(i)%coeffs = s_current%coeffs

!!$       print *, 'Surface: ', surfaces(i)%id
!!$       print *, '   Type   ', surfaces(i)%type
!!$       print *, '   Coeffs ', surfaces(i)%coeffs
       
       s_head => s_head%next
       deallocate( s_current )
       s_current => s_head
    end do

    ! Convert linked list of Materials into normal array
    

  end subroutine read_input

!------------------------------------------------------------------------------

  subroutine parse_cell( cell, words, n_words )

    type(CellList), pointer, intent(inout) :: cell
    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

    integer :: readError
    integer :: i
    character(250) :: msg
    character(32) :: word

    ! Read cell identifier
    read(words(2), fmt='(I8)', iostat=readError) cell%id
    if ( readError > 0 ) then
       msg = "Invalid cell name: " // words(2)
       call error( msg )
    end if

    ! Read cell material
    read(words(3), fmt='(I8)', iostat=readError) cell%material

    ! Read list of surfaces
    allocate( cell%boundary_list(n_words-3) )
    do i = 1, n_words-3
       word = words(i+3)
       if ( word(1:1) == '(' ) then
          cell%boundary_list(i) = OP_LEFT_PAREN
       elseif ( word(1:1) == ')' ) then
          cell%boundary_list(i) = OP_RIGHT_PAREN
       elseif ( word(1:1) == ':' ) then
          cell%boundary_list(i) = OP_UNION
       elseif ( word(1:1) == '#' ) then
          cell%boundary_list(i) = OP_DIFFERENCE
       else
          read(word, fmt='(I8)', iostat=readError) cell%boundary_list(i)
       end if
    end do

  end subroutine parse_cell

!------------------------------------------------------------------------------

  subroutine parse_surface( surface, words, n_words )

    type(SurfaceList), pointer, intent(inout) :: surface
    character(*), intent(in) :: words(max_words)
    integer, intent(in) :: n_words

    integer :: readError
    integer :: i
    integer :: coeffs_reqd
    character(250) :: msg
    character(32) :: word

    ! Read surface identifier
    read(words(2), fmt='(I8)', iostat=readError) surface%id
    if ( readError > 0 ) then
       msg = "Invalid surface name: " // words(2)
       call error( msg )
    end if

    ! Read surface type
    word = words(3)
    call lower_case(word)
    select case ( trim(word) )
       case ( 'px' )
          surface%type = SURF_PX
          coeffs_reqd  = 1
       case ( 'py' )
          surface%type = SURF_PY
          coeffs_reqd  = 1
       case ( 'pz' )
          surface%type = SURF_PZ
          coeffs_reqd  = 1
       case ( 'plane' )
          surface%type = SURF_PLANE
          coeffs_reqd  = 4
       case ( 'cylx' )
          surface%type = SURF_CYL_X
          coeffs_reqd  = 3
       case ( 'cyly' )
          surface%type = SURF_CYL_Y
          coeffs_reqd  = 3
       case ( 'cylz' )
          surface%type = SURF_CYL_Z
          coeffs_reqd  = 3
       case ( 'sph' )
          surface%type = SURF_SPHERE
          coeffs_reqd  = 4
       case ( 'boxx' )
          surface%type = SURF_BOX_X
          coeffs_reqd  = 4
       case ( 'boxy' )
          surface%type = SURF_BOX_Y
          coeffs_reqd  = 4
       case ( 'boxz' )
          surface%type = SURF_BOX_Z
          coeffs_reqd  = 4
       case ( 'box' ) 
          surface%type = SURF_BOX
          coeffs_reqd  = 6
       case ( 'gq' )
          surface%type = SURF_GQ
          coeffs_reqd  = 10
       case default
          msg = "Invalid surface type: " // words(3)
          call error( msg )
    end select

    ! Make sure there are enough coefficients for surface type
    if ( n_words-3 < coeffs_reqd ) then
       msg = "Not enough coefficients for surface: " // words(2)
       call error( msg )
    end if
    
    ! Read list of surfaces
    allocate( surface%coeffs(n_words-3) )
    do i = 1, n_words-3
       surface%coeffs(i) = string_to_real(words(i+3))
    end do

  end subroutine parse_surface

!------------------------------------------------------------------------------

end module fileio

    

    
