module cross_section

  use global
  use string, only: split_string, str_to_real
  use data_structures, only: dict_create, dict_add_key, dict_get_key
  use output, only: error, message
  use types, only: xsData

  implicit none

contains

!=====================================================================
! READ_XSDATA reads the data in a SERPENT xsdata file and builds a
! dictionary to find cross-section information later on.
!=====================================================================

  subroutine read_xsdata( path )

    character(*), intent(in) :: path

    type(xsData), pointer :: iso => null()
    character(250) :: line, msg
    character(100) :: words(max_words)
    character(100) :: filename
    integer :: i, n
    integer :: in = 7
    logical :: file_exists
    character(7) :: readable
    integer :: count
    integer :: index
    integer :: ioError

    msg = "Reading cross-section summary file..."
    call message(msg, 5)

    ! Construct filename
    n = len(path)
    if (path(n:n) == '/') then
       filename = trim(path) // 'xsdata'
    else
       filename = trim(path) // '/xsdata'
    end if

    ! Check if xsdata exists and is readable
    inquire( FILE=filename, EXIST=file_exists, READ=readable )
    if ( .not. file_exists ) then
       msg = "Cross section summary '" // trim(filename) // "' does not exist!"
       call error( msg )
    elseif ( readable(1:3) == 'NO' ) then
       msg = "Cross section summary '" // trim(filename) // "' is not readable!" &
            & // "Change file permissions with chmod command."
       call error( msg )
    end if

    ! open xsdata file
    open(file=filename, unit=in, status='old', &
         & action='read', iostat=ioError)
    if (ioError /= 0) then
       msg = "Error while opening file: " // filename
       call error(msg)
    end if

    ! determine how many lines
    count = 0
    do
       read(unit=in, fmt='(A250)', iostat=ioError) line
       if ( ioError < 0 ) then
          ! reached end of file
          exit
       elseif ( ioError > 0 ) then
          msg = "Unknown error while reading file: " // filename
          close(unit=in)
          call error( msg )
       end if
       count = count + 1
    end do
    allocate(xsdatas(count))

    ! create dictionary
    call dict_create(xsdata_dict)

    ! read actual lines
    index = 0
    rewind(in)
    do
       read(unit=in, fmt='(A250)', iostat=ioError) line
       if ( ioError < 0 ) exit
       index = index + 1
       call split_string(line, words, n)
       if ( n == 0 ) cycle ! skip blank line

       ! Check to make sure there are enough arguments
       if ( n < 9 ) then
          msg = "Not enough arguments on xsdata line: " // line
          close(unit=in)
          call error( msg )
       end if

       iso => xsdatas(index)

       ! store data
       iso%alias = words(1)
       iso%id = words(2)
       iso%type = str_to_int(words(3))
       iso%zaid = str_to_int(words(4))
       iso%isomeric = str_to_int(words(5))
       iso%awr = str_to_real(words(6))
       iso%temp = str_to_real(words(7))
       iso%binary = str_to_int(words(8))
       iso%path = words(9)

       ! create dictionary entry
       call dict_add_key(xsdata_dict, iso%alias, index)
    end do

    close(unit=in)

  end subroutine read_xsdata

!=====================================================================
! MATERIAL_TOTAL_XS calculates the total macroscopic cross section of
! each material for use in sampling path lengths
!=====================================================================

  subroutine material_total_xs()

    type(Material),      pointer :: mat => null()
    type(AceContinuous), pointer :: table => null()
    integer :: i, j, k
    integer :: index
    integer :: IE
    real(8) :: xs
    real(8) :: density
    real(8) :: density_i
    real(8) :: val
    real(8) :: r
    real(8) :: E_i
    real(8) :: E_i1
    real(8) :: sigma_i
    real(8) :: sigma_i1
    character(250) :: msg

    msg = "Creating material total cross-sections..."
    call message(msg, 4)

    ! loop over all materials
    do i = 1, n_materials

       ! allocate storage for total xs
       mat => materials(i)
       density = mat%atom_density
       allocate(mat%total_xs(n_grid))

       ! initialize total cross-section
       mat % total_xs = ZERO

       ! loop over each isotope
       do j = 1, mat%n_isotopes
          index = mat % table(j)
          table => xs_continuous(index)

          ! find atom density of isotope
          density_i = density * mat % atom_percent(j)

          ! loop over points in union energy grid
          do k = 1, n_grid
             ! find nuclide grid index
             IE      = table % grid_index(k)

             ! get nuclide grid energy and cross-section value
             E_i     = table % energy(IE)
             sigma_i = table % sigma_t(IE)

             ! interpolate as necessary
             if (E_i < e_grid(k)) then
                E_i1     = table % energy(IE + 1)
                sigma_i1 = table % sigma_t(IE + 1)

                r = (e_grid(k) - E_i)/(E_i1 - E_i)
                val = (ONE - r)*sigma_i + r*sigma_i1
             else
                val = sigma_i
             end if

             ! add contribution to total cross-section
             mat%total_xs(k) = mat%total_xs(k) + density_i * val
          end do
       end do
    end do

  end subroutine material_total_xs

end module cross_section
