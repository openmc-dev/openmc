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

    ! Initialize FORTRAN interface.
    call h5open_f(hdf5_err)

    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_output_file, hdf5_err)

  end subroutine hdf5_create_output

!===============================================================================
! HDF5_WRITE_HEADER
!===============================================================================

  subroutine hdf5_write_header()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    character(8)     :: date_
    character(10)    :: time_
    character(19)    :: current_time
     
    ! Write version information
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_major", &
         rank, dims, (/ VERSION_MAJOR /), hdf5_err)
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_minor", &
         rank, dims, (/ VERSION_MINOR /), hdf5_err)
    call h5ltmake_dataset_int_f(hdf5_output_file, "version_release", &
         rank, dims, (/ VERSION_RELEASE /), hdf5_err)

    ! Write current date and time
    call date_and_time(DATE=date_, TIME=time_)
    current_time = date_(1:4) // "-" // date_(5:6) // "-" // date_(7:8) // &
         " " // time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)
    call h5ltmake_dataset_string_f(hdf5_output_file, "/date_and_time", &
         current_time, hdf5_err)

    ! Write MPI information
    call h5ltmake_dataset_int_f(hdf5_output_file, "n_procs", &
         rank, dims, (/ n_procs /), hdf5_err)
    call h5ltset_attribute_string_f(hdf5_output_file, "n_procs", &
         "description", "Number of MPI processes", hdf5_err)

  end subroutine hdf5_write_header

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
     
    ! Write criticality information
    if (problem_type == PROB_CRITICALITY) then
       ! Need to write integer(8)'s using double instead since there is no H5LT
       ! call for making a dataset of type long
       call h5ltmake_dataset_double_f(hdf5_output_file, "n_particles", &
            rank, dims, (/ real(n_particles,8) /), hdf5_err)

       ! Use H5LT interface to write n_cycles, n_inactive, and n_active
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_cycles", &
            rank, dims, (/ n_cycles /), hdf5_err)
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_inactive", &
            rank, dims, (/ n_inactive /), hdf5_err)
       call h5ltmake_dataset_int_f(hdf5_output_file, "n_active", &
            rank, dims, (/ n_cycles - n_inactive /), hdf5_err)

       ! Add description of each variable
       call h5ltset_attribute_string_f(hdf5_output_file, "n_particles", &
            "description", "Number of particles per cycle", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_cycles", &
            "description", "Total number of cycles", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_inactive", &
            "description", "Number of inactive cycles", hdf5_err)
       call h5ltset_attribute_string_f(hdf5_output_file, "n_active", &
            "description", "Number of active cycles", hdf5_err)
    end if

    call hdf5_write_geometry()
    call hdf5_write_materials()

  end subroutine hdf5_write_summary

!===============================================================================
! HDF5_WRITE_GEOMETRY
!===============================================================================

  subroutine hdf5_write_geometry()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: geometry_group
     
    ! Create group for geometry
    call h5gcreate_f(hdf5_output_file, "/geometry", geometry_group, hdf5_err)

    ! Use H5LT interface to write number of geometry objects
    call h5ltmake_dataset_int_f(geometry_group, "n_cells", &
         rank, dims, (/ n_cells /), hdf5_err)
    call h5ltmake_dataset_int_f(geometry_group, "n_surfaces", &
         rank, dims, (/ n_surfaces /), hdf5_err)
    call h5ltmake_dataset_int_f(geometry_group, "n_universes", &
         rank, dims, (/ n_universes /), hdf5_err)
    call h5ltmake_dataset_int_f(geometry_group, "n_lattices", &
         rank, dims, (/ n_lattices /), hdf5_err)

    ! Close geometry group
    call h5gclose_f(geometry_group, hdf5_err)

  end subroutine hdf5_write_geometry

!===============================================================================
! HDF5_WRITE_MATERIALS
!===============================================================================

  subroutine hdf5_write_materials()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: materials_group
     
    ! Create group for materials
    call h5gcreate_f(hdf5_output_file, "/materials", materials_group, hdf5_err)

    ! Use H5LT interface to write number of materials
    call h5ltmake_dataset_int_f(materials_group, "n_materials", &
         rank, dims, (/ n_materials /), hdf5_err)

    ! Close materials group
    call h5gclose_f(materials_group, hdf5_err)

  end subroutine hdf5_write_materials

!===============================================================================
! HDF5_WRITE_TIMING
!===============================================================================

  subroutine hdf5_write_timing()

    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
    integer(HID_T)   :: timing_group
    integer(8)       :: total_particles
    real(8)          :: speed

    ! Create group for timing
    call h5gcreate_f(hdf5_output_file, "/timing", timing_group, hdf5_err)

    ! Write timing data
    call h5ltmake_dataset_double_f(timing_group, "time_initialize", &
         rank, dims, (/ time_initialize % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_read_xs", &
         rank, dims, (/ time_read_xs % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_unionize", &
         rank, dims, (/ time_unionize % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_compute", &
         rank, dims, (/ time_compute % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_intercycle", &
         rank, dims, (/ time_intercycle % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_tallies", &
         rank, dims, (/ time_ic_tallies % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_sample", &
         rank, dims, (/ time_ic_sample % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_sendrecv", &
         rank, dims, (/ time_ic_sendrecv % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_rebuild", &
         rank, dims, (/ time_ic_rebuild % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_inactive", &
         rank, dims, (/ time_inactive % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_active", &
         rank, dims, (/ time_active % elapsed /), hdf5_err)
    call h5ltmake_dataset_double_f(timing_group, "time_total", &
         rank, dims, (/ time_total % elapsed /), hdf5_err)

    ! Add descriptions to timing data
    call h5ltset_attribute_string_f(timing_group, "time_initialize", &
         "description", "Total time elapsed for initialization (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_read_xs", &
         "description", "Time reading cross-section libraries (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_unionize", &
         "description", "Time unionizing energy grid (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_compute", &
         "description", "Total time in computation (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_intercycle", &
         "description", "Total time between cycles (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_tallies", &
         "description", "Time between cycles accumulating tallies (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_sample", &
         "description", "Time between cycles sampling source sites (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_sendrecv", &
         "description", "Time between cycles SEND/RECVing source sites (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_rebuild", &
         "description", "Time between cycles reconstructing source bank (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_inactive", &
         "description", "Total time in inactive cycles (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_active", &
         "description", "Total time in active cycles (s)", hdf5_err)
    call h5ltset_attribute_string_f(timing_group, "time_total", &
         "description", "Total time elapsed (s)", hdf5_err)

    ! Write calculation rate
    total_particles = n_particles * n_cycles
    speed = real(total_particles) / time_compute % elapsed
    call h5ltmake_dataset_double_f(timing_group, "neutrons_per_second", &
         rank, dims, (/ speed /), hdf5_err)

    ! Close timing group
    call h5gclose_f(timing_group, hdf5_err)

  end subroutine hdf5_write_timing

!===============================================================================
! HDF5_CLOSE_OUTPUT
!===============================================================================

  subroutine hdf5_close_output()

    ! Terminate access to the file.
    call h5fclose_f(hdf5_output_file, hdf5_err)

    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)

  end subroutine hdf5_close_output

#endif

end module hdf5_interface
