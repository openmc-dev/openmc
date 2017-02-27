module particle_restart_write

  use bank_header,     only: Bank
  use global
  use hdf5_interface
  use particle_header, only: Particle
  use string,          only: to_str

  use hdf5

  implicit none
  private
  public ::  write_particle_restart

contains

!===============================================================================
! WRITE_PARTICLE_RESTART is the main routine that writes out the particle file
!===============================================================================

  subroutine write_particle_restart(p)
    type(Particle), intent(in) :: p

    integer(HID_T) :: file_id
    character(MAX_FILE_LEN) :: filename

    ! Dont write another restart file if in particle restart mode
    if (run_mode == MODE_PARTICLE) return

    ! Set up file name
    filename = trim(path_output) // 'particle_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(p%id)) // '.h5'

!$omp critical (WriteParticleRestart)
    ! Create file
    file_id = file_create(filename)

    associate (src => source_bank(current_work))
      ! Write filetype and version info
      call write_attribute(file_id, 'filetype', 'particle restart')
      call write_attribute(file_id, 'version', VERSION_PARTICLE_RESTART)
      call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
      call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

      ! Write data to file
      call write_dataset(file_id, 'current_batch', current_batch)
      call write_dataset(file_id, 'generations_per_batch', gen_per_batch)
      call write_dataset(file_id, 'current_generation', current_gen)
      call write_dataset(file_id, 'n_particles', n_particles)
      select case(run_mode)
      case (MODE_FIXEDSOURCE)
        call write_dataset(file_id, 'run_mode', 'fixed source')
      case (MODE_EIGENVALUE)
        call write_dataset(file_id, 'run_mode', 'eigenvalue')
      case (MODE_PARTICLE)
        call write_dataset(file_id, 'run_mode', 'particle restart')
      end select
      call write_dataset(file_id, 'id', p%id)
      call write_dataset(file_id, 'weight', src%wgt)
      call write_dataset(file_id, 'energy', src%E)
      call write_dataset(file_id, 'xyz', src%xyz)
      call write_dataset(file_id, 'uvw', src%uvw)
    end associate

    ! Close file
    call file_close(file_id)
!$omp end critical (WriteParticleRestart)

  end subroutine write_particle_restart

end module particle_restart_write
