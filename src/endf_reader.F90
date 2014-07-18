module endf_reader

  use error, only: fatal_error
  use global
  use output,only: write_message

  integer, allocatable      :: urr_zaids(:)  ! ZAID #'s for URR nuclides
  type ENDFFilename
    character(MAX_LINE_LEN) :: filename      ! ENDF filename for URR nuclide
  end type ENDFFilename
  type(ENDFFilename), allocatable :: urr_endf_filenames(:) ! ENDF filename

contains

!===============================================================================
! READ_ENDF reads an ENDF file
!===============================================================================

  subroutine read_endf()

    integer       :: i           ! index over nuclides w/ an otf URR treatment
    logical       :: file_exists ! does ENDF file exist?
    character(7)  :: readable    ! is ENDF file readable?
    integer       :: in = 11     ! input unit
    character(80) :: record

    do i = 1, n_otf_urr_nuclides

      inquire(file = trim(path_endf)//trim(urr_endf_filenames(i)%filename), &
        & exist = file_exists, read = readable)
      if (.not. file_exists) then
        message = 'ENDF file '//trim(urr_endf_filenames(i)%filename) &
          //' does not exist.'
        call fatal_error()
      else if (readable(1:3) == 'NO') then
        message = 'ENDF file '//trim(urr_endf_filenames(i)%filename)// &
          & ' is not readable.  Change file permissions with chmod command.'
        call fatal_error()
      end if

      ! display message
      message = "Loading ENDF file: "//trim(urr_endf_filenames(i)%filename)
      call write_message(6)

      open(unit = in, &
           file = trim(path_endf)//trim(urr_endf_filenames(i)%filename))

      read(in, 10) record
10    format(A80)
      write(*,10) record

      close(in)

    end do

  end subroutine read_endf

end module endf_reader
