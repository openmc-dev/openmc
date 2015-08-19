module particle_restart_write

  use bank_header,      only: Bank
  use global
  use output_interface, only: BinaryOutput
  use particle_header,  only: Particle
  use string,           only: to_str

  implicit none
  private
  public ::  write_particle_restart

  ! Binary output file
  type(BinaryOutput) :: pr

contains

!===============================================================================
! WRITE_PARTICLE_RESTART is the main routine that writes out the particle file
!===============================================================================

  subroutine write_particle_restart(p)

    type(Particle), intent(in) :: p

    character(MAX_FILE_LEN) :: filename
    type(Bank), pointer     :: src => null()

    ! Dont write another restart file if in particle restart mode
    if (run_mode == MODE_PARTICLE) return

    ! Set up file name
    filename = trim(path_output) // 'particle_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(p % id))
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

!$omp critical (WriteParticleRestart)
    ! Create file
    call pr % file_create(filename)

    ! Get information about source particle
    select case (run_mode)
    case (MODE_EIGENVALUE)
      src => source_bank(current_work)
    case (MODE_FIXEDSOURCE)
      src => source_site
    end select

    ! Write data to file
    call pr % write_data(FILETYPE_PARTICLE_RESTART, 'filetype')
    call pr % write_data(REVISION_PARTICLE_RESTART, 'revision')
    call pr % write_data(current_batch, 'current_batch')
    call pr % write_data(gen_per_batch, 'gen_per_batch')
    call pr % write_data(current_gen, 'current_gen')
    call pr % write_data(n_particles, 'n_particles')
    call pr % write_data(run_mode, 'run_mode')
    call pr % write_data(p % id, 'id')
    call pr % write_data(src % wgt, 'weight')
    call pr % write_data(src % E, 'energy')
    call pr % write_data(src % xyz, 'xyz', length = 3)
    call pr % write_data(src % uvw, 'uvw', length = 3)

    ! Close file
    call pr % file_close()
!$omp end critical (WriteParticleRestart)

  end subroutine write_particle_restart

end module particle_restart_write
