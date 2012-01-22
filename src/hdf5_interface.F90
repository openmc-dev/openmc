module hdf5_interface

  use constants
  use global

#ifdef HDF5
  use hdf5
  use h5lt
#endif

  implicit none

contains

#ifdef HDF5

!===============================================================================
! HDF5_CREATE_OUTPUT
!===============================================================================

  subroutine hdf5_create_output()

    character(9), parameter :: filename = "output.h5" ! File name
    integer :: error  ! Error flag

    ! Initialize FORTRAN interface.
    call h5open_f (error)

    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_output_file, error)

  end subroutine hdf5_create_output

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    integer          :: error
    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    character(8)     :: date_
    character(10)    :: time_
    character(19)    :: current_time
     
    ! Write version information
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_major", &
         rank, dims, (/ VERSION_MAJOR /), error)
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_minor", &
         rank, dims, (/ VERSION_MINOR /), error)
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_release", &
         rank, dims, (/ VERSION_RELEASE /), error)

    ! Write current date and time
    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)
    call h5ltmake_dataset_string_f(hdf5_output_file, "/date_and_time", &
         current_time, error)

    ! Write MPI information
    call h5ltmake_dataset_int_f(hdf5_output_file, "n_procs", &
         rank, dims, (/ n_procs /), error)
    call h5ltset_attribute_string_f(hdf5_output_file, "n_procs", &
         "description", "Number of MPI processes", error)

    ! Write criticality information
    if (problem_type == PROB_CRITICALITY) then
       ! Need to write integer(8)'s using double instead since there is no H5LT
       ! call for making a dataset of type long
       call h5ltmake_dataset_double_f(hdf5_output_file, "n_particles", &
            rank, dims, (/ real(n_particles,8) /), error)

       ! Use H5LT interface to write n_cycles, n_inactive, and n_active
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_cycles", &
            rank, dims, (/ n_cycles /), error)
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_inactive", &
            rank, dims, (/ n_inactive /), error)
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_active", &
            rank, dims, (/ n_cycles - n_inactive /), error)

       ! Add description of each variable
       call h5ltset_attribute_string_f(hdf5_output_file, "n_particles", &
            "description", "Number of particles per cycle", error)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_cycles", &
            "description", "Total number of cycles", error)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_inactive", &
            "description", "Number of inactive cycles", error)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_active", &
            "description", "Number of active cycles", error)
    end if

  end subroutine hdf5_write_summary

!===============================================================================
! HDF5_WRITE_TIMING
!===============================================================================

  subroutine hdf5_write_timing()

    integer          :: error
    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: timing_group
    integer(8)       :: total_particles
    real(8)          :: speed

    ! Create group for timing
    call h5gcreate_f(hdf5_output_file, "/timing", timing_group, error)

    ! Write timing data
    call h5ltmake_dataset_double_f(timing_group, "time_initialize", &
         rank, dims, (/ time_initialize % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_read_xs", &
         rank, dims, (/ time_read_xs % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_unionize", &
         rank, dims, (/ time_unionize % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_compute", &
         rank, dims, (/ time_compute % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_intercycle", &
         rank, dims, (/ time_intercycle % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_tallies", &
         rank, dims, (/ time_ic_tallies % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_sample", &
         rank, dims, (/ time_ic_sample % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_sendrecv", &
         rank, dims, (/ time_ic_sendrecv % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_rebuild", &
         rank, dims, (/ time_ic_rebuild % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_inactive", &
         rank, dims, (/ time_inactive % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_active", &
         rank, dims, (/ time_active % elapsed /), error)
    call h5ltmake_dataset_double_f(timing_group, "time_total", &
         rank, dims, (/ time_total % elapsed /), error)

    ! Add descriptions to timing data
    call h5ltset_attribute_string_f(timing_group, "time_initialize", &
         "description", "Total time elapsed for initialization (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_read_xs", &
         "description", "Time reading cross-section libraries (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_unionize", &
         "description", "Time unionizing energy grid (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_compute", &
         "description", "Total time in computation (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_intercycle", &
         "description", "Total time between cycles (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_tallies", &
         "description", "Time between cycles accumulating tallies (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_sample", &
         "description", "Time between cycles sampling source sites (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_sendrecv", &
         "description", "Time between cycles SEND/RECVing source sites (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_rebuild", &
         "description", "Time between cycles reconstructing source bank (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_inactive", &
         "description", "Total time in inactive cycles (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_active", &
         "description", "Total time in active cycles (s)", error)
    call h5ltset_attribute_string_f(timing_group, "time_total", &
         "description", "Total time elapsed (s)", error)

    ! Write calculation rate
    total_particles = n_particles * n_cycles
    speed = real(total_particles) / time_compute % elapsed
    call h5ltmake_dataset_double_f(timing_group, "neutrons_per_second", &
         rank, dims, (/ speed /), error)

    ! Close timing group
    call h5gclose_f(timing_group, error)

  end subroutine hdf5_write_timing

!===============================================================================
! HDF5_CLOSE_OUTPUT
!===============================================================================

  subroutine hdf5_close_output()

    integer :: error  ! Error flag

    ! Terminate access to the file.
    call h5fclose_f(hdf5_output_file, error)

    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine hdf5_close_output

#endif

end module hdf5_interface
