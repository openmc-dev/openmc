module particle_track

  use constants
  use global
  !use output_interface  only: file_create, file_close, write_double1Darray
  use string,           only: to_str

  implicit none

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK opens a particle track output file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.  It
! should also probably write a header that identifies the file as a particle
! track and maybe adds particle identifying information (batch #, etc.).
!===============================================================================

  subroutine initialize_particle_track()
    character(MAX_FILE_LEN) :: filename

    filename = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.binary'
    open(UNIT=UNIT_TRACK, FILE=filename, STATUS='replace', &
         ACCESS='stream')
    !file_create(filename, 'serial')

  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK outputs particle position to a binary file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.  Perhaps
! it should also be made somehow more general so that it can output
! information other than just particle position.
!===============================================================================

  subroutine write_particle_track()
    write(UNIT_TRACK) p % coord0 % xyz
    !write_double1Darray(p % coord0 % xyz, 'coordinates', 3)
  end subroutine write_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK closes the particle track file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.
!===============================================================================

  subroutine finalize_particle_track()
    close(UNIT=UNIT_TRACK)
    !file_close('serial')
  end subroutine finalize_particle_track

end module particle_track
