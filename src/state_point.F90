module state_point

!===============================================================================
! STATE_POINT -- This module handles writing and reading state point
! files. State points are contain complete tally results, source sites, and
! various other data. They can be used to restart a run or to reconstruct
! confidence intervals for tallies (this requires post-processing via Python
! scripts).
!
! State points can be written at any batch during a simulation, or at specified
! intervals, using the <state_point ... /> tag.
!===============================================================================

  use, intrinsic :: ISO_C_BINDING

  use bank_header,        only: Bank
  use constants
  use endf,               only: reaction_name
  use error,              only: fatal_error, warning, write_message
  use hdf5_interface
  use message_passing
  use mgxs_interface
  use nuclide_header,     only: nuclides
  use output,             only: time_stamp
  use random_lcg,         only: openmc_get_seed, openmc_set_seed
  use settings
  use simulation_header
  use string,             only: to_str, count_digits, zero_padded, to_f_string
  use tally_header
  use tally_filter_header
  use tally_derivative_header, only: tally_derivs
  use timer_header

  implicit none

  interface
    subroutine write_source_bank(group_id) bind(C)
      import HID_T
      integer(HID_T), value :: group_id
    end subroutine write_source_bank

    subroutine read_source_bank(group_id) bind(C)
      import HID_T
      integer(HID_T), value :: group_id
    end subroutine read_source_bank
  end interface

contains

!===============================================================================
! OPENMC_STATEPOINT_WRITE writes an HDF5 statepoint file to disk
!===============================================================================

  subroutine statepoint_write_f(file_id) bind(C)
    integer(HID_T), value :: file_id

    integer :: i, j
    integer :: i_xs
    integer, allocatable :: id_array(:)
    integer(HID_T) :: tallies_group, tally_group, &
                      filters_group, filter_group, derivs_group, &
                      deriv_group
    character(MAX_WORD_LEN), allocatable :: str_array(:)
    character(MAX_WORD_LEN, kind=C_CHAR) :: temp_name

    ! Open tallies group
    tallies_group = open_group(file_id, "tallies")

    ! Write information for derivatives.
    if (size(tally_derivs) > 0) then
      derivs_group = create_group(tallies_group, "derivatives")
      do i = 1, size(tally_derivs)
        associate(deriv => tally_derivs(i))
          deriv_group = create_group(derivs_group, "derivative " &
               // trim(to_str(deriv % id)))
          select case (deriv % variable)
          case (DIFF_DENSITY)
            call write_dataset(deriv_group, "independent variable", "density")
            call write_dataset(deriv_group, "material", deriv % diff_material)
          case (DIFF_NUCLIDE_DENSITY)
            call write_dataset(deriv_group, "independent variable", &
                 "nuclide_density")
            call write_dataset(deriv_group, "material", deriv % diff_material)
            call write_dataset(deriv_group, "nuclide", &
                 nuclides(deriv % diff_nuclide) % name)
          case (DIFF_TEMPERATURE)
            call write_dataset(deriv_group, "independent variable", &
                 "temperature")
            call write_dataset(deriv_group, "material", deriv % diff_material)
          case default
            call fatal_error("Independent variable for derivative " &
                 // trim(to_str(deriv % id)) // " not defined in &
                 &state_point.F90.")
          end select
          call close_group(deriv_group)
        end associate
      end do
      call close_group(derivs_group)
    end if

    ! Write number of filters
    filters_group = create_group(tallies_group, "filters")
    call write_attribute(filters_group, "n_filters", n_filters)

    if (n_filters > 0) then
      ! Write IDs of filters
      allocate(id_array(n_filters))
      do i = 1, n_filters
        id_array(i) = filters(i) % obj % id
      end do
      call write_attribute(filters_group, "ids", id_array)
      deallocate(id_array)

      ! Write filter information
      FILTER_LOOP: do i = 1, n_filters
        filter_group = create_group(filters_group, "filter " // &
             trim(to_str(filters(i) % obj % id)))
        call filters(i) % obj % to_statepoint(filter_group)
        call close_group(filter_group)
      end do FILTER_LOOP
    end if

    call close_group(filters_group)

      ! Write number of tallies
    call write_attribute(tallies_group, "n_tallies", n_tallies)

    if (n_tallies > 0) then
      ! Write array of tally IDs
      allocate(id_array(n_tallies))
      do i = 1, n_tallies
        id_array(i) = tallies(i) % obj % id
      end do
      call write_attribute(tallies_group, "ids", id_array)
      deallocate(id_array)

      ! Write all tally information except results
      TALLY_METADATA: do i = 1, n_tallies

        ! Get pointer to tally
        associate (tally => tallies(i) % obj)
        tally_group = create_group(tallies_group, "tally " // &
             trim(to_str(tally % id)))

        ! Write the name for this tally
        call write_dataset(tally_group, "name", tally % name)

        select case(tally % estimator)
        case (ESTIMATOR_ANALOG)
          call write_dataset(tally_group, "estimator", "analog")
        case (ESTIMATOR_TRACKLENGTH)
          call write_dataset(tally_group, "estimator", "tracklength")
        case (ESTIMATOR_COLLISION)
          call write_dataset(tally_group, "estimator", "collision")
        end select
        call write_dataset(tally_group, "n_realizations", &
             tally % n_realizations)

        call write_dataset(tally_group, "n_filters", size(tally % filter))
        if (size(tally % filter) > 0) then
          ! Write IDs of filters
          allocate(id_array(size(tally % filter)))
          do j = 1, size(tally % filter)
            id_array(j) = filters(tally % filter(j)) % obj % id
          end do
          call write_dataset(tally_group, "filters", id_array)
          deallocate(id_array)
        end if

        ! Set up nuclide bin array and then write
        allocate(str_array(tally % n_nuclide_bins))
        NUCLIDE_LOOP: do j = 1, tally % n_nuclide_bins
          if (tally % nuclide_bins(j) > 0) then
            if (run_CE) then
              i_xs = index(nuclides(tally % nuclide_bins(j)) % name, '.')
              if (i_xs > 0) then
                str_array(j) = nuclides(tally % nuclide_bins(j)) % name(1 : i_xs-1)
              else
                str_array(j) = nuclides(tally % nuclide_bins(j)) % name
              end if
            else
              call get_name_c(tally % nuclide_bins(j), len(temp_name), &
                              temp_name)
              i_xs = index(temp_name, '.')
              if (i_xs > 0) then
                str_array(j) = trim(temp_name(1 : i_xs-1))
              else
                str_array(j) = trim(temp_name)
              end if
            end if
          else
            str_array(j) = 'total'
          end if
        end do NUCLIDE_LOOP
        call write_dataset(tally_group, "nuclides", str_array)
        deallocate(str_array)

        ! Write derivative information.
        if (tally % deriv /= NONE) then
          call write_dataset(tally_group, "derivative", &
               tally_derivs(tally % deriv) % id)
        end if

        ! Write scores.
        call write_dataset(tally_group, "n_score_bins", tally % n_score_bins)
        allocate(str_array(size(tally % score_bins)))
        do j = 1, size(tally % score_bins)
          str_array(j) = reaction_name(tally % score_bins(j))
        end do
        call write_dataset(tally_group, "score_bins", str_array)

        deallocate(str_array)

        call close_group(tally_group)
        end associate
      end do TALLY_METADATA
    end if

    if (reduce_tallies) then
      ! Write global tallies
      call write_dataset(file_id, "global_tallies", global_tallies)

      ! Write tallies
      if (active_tallies % size() > 0) then
        ! Indicate that tallies are on
        call write_attribute(file_id, "tallies_present", 1)

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          associate (tally => tallies(i) % obj)
            ! Write sum and sum_sq for each bin
            tally_group = open_group(tallies_group, "tally " &
                 // to_str(tally % id))
            call tally % write_results_hdf5(tally_group)
            call close_group(tally_group)
          end associate
        end do TALLY_RESULTS
      else
        ! Indicate tallies are off
        call write_attribute(file_id, "tallies_present", 0)
      end if
    end if

    call close_group(tallies_group)
  end subroutine

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point() bind(C)

    integer :: i
    integer :: int_array(3)
    integer, allocatable :: array(:)
    integer(C_INT64_T) :: seed
    integer(HID_T) :: file_id
    integer(HID_T) :: tallies_group
    integer(HID_T) :: tally_group
    logical :: source_present
    character(MAX_WORD_LEN) :: word

    interface
      subroutine read_eigenvalue_hdf5(group) bind(C)
        import HID_T
        integer(HID_T), value :: group
      end subroutine

      subroutine restart_set_keff() bind(C)
      end subroutine
    end interface

    ! Write message
    call write_message("Loading state point " // trim(path_state_point) &
         // "...", 5)

    ! Open file for reading
    file_id = file_open(path_state_point, 'r', parallel=.true.)

    ! Read filetype
    call read_attribute(word, file_id, "filetype")
    if (word /= 'statepoint') then
      call fatal_error("OpenMC tried to restart from a non-statepoint file.")
    end if

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call read_attribute(array, file_id, "version")
    if (any(array /= VERSION_STATEPOINT)) then
      call fatal_error("State point version does not match current version &
           &in OpenMC.")
    end if

    ! Read and overwrite random number seed
    call read_dataset(seed, file_id, "seed")
    call openmc_set_seed(seed)

    ! It is not impossible for a state point to be generated from a CE run but
    ! to be loaded in to an MG run (or vice versa), check to prevent that.
    call read_dataset(word, file_id, "energy_mode")
    if (word == "multi-group" .and. run_CE) then
      call fatal_error("State point file is from multi-group run but &
                       & current run is continous-energy!")
    else if (word == "continuous-energy" .and. .not. run_CE) then
      call fatal_error("State point file is from continuous-energy run but &
                       & current run is multi-group!")
    end if

    ! Read and overwrite run information except number of batches
    call read_dataset(word, file_id, "run_mode")
    select case(word)
    case ('fixed source')
      run_mode = MODE_FIXEDSOURCE
    case ('eigenvalue')
      run_mode = MODE_EIGENVALUE
    end select
    call read_attribute(int_array(1), file_id, "photon_transport")
    if (int_array(1) == 1) then
      photon_transport = .true.
    else
      photon_transport = .false.
    end if
    call read_dataset(n_particles, file_id, "n_particles")
    call read_dataset(int_array(1), file_id, "n_batches")

    ! Take maximum of statepoint n_batches and input n_batches
    n_batches = max(n_batches, int_array(1))

    ! Read batch number to restart at
    call read_dataset(restart_batch, file_id, "current_batch")

    ! Check for source in statepoint if needed
    call read_attribute(int_array(1), file_id, "source_present")
    if (int_array(1) == 1) then
      source_present = .true.
    else
      source_present = .false.
    end if

    if (restart_batch > n_batches) then
      call fatal_error("The number batches specified in settings.xml is fewer &
           & than the number of batches in the given statepoint file.")
    end if

    ! Read information specific to eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
      call read_dataset(int_array(1), file_id, "n_inactive")
      call read_eigenvalue_hdf5(file_id)

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))
    end if

    ! Read number of realizations for global tallies
    call read_dataset(n_realizations, file_id, "n_realizations", indep=.true.)

    ! Set k_sum, keff, and current_batch based on whether restart file is part
    ! of active cycle or inactive cycle
    call restart_set_keff()
    current_batch = restart_batch

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      call fatal_error("Source bank must be contained in statepoint restart &
           &file")
    end if

    ! Read tallies to master. If we are using Parallel HDF5, all processes
    ! need to be included in the HDF5 calls.
#ifdef PHDF5
    if (.true.) then
#else
    if (master) then
#endif
      ! Read global tally data
      call read_dataset(global_tallies, file_id, "global_tallies")

      ! Check if tally results are present
      call read_attribute(int_array(1), file_id, "tallies_present")

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        tallies_group = open_group(file_id, "tallies")

        TALLY_RESULTS: do i = 1, n_tallies
          associate (t => tallies(i) % obj)
            ! Read sum, sum_sq, and N for each bin
            tally_group = open_group(tallies_group, "tally " // &
                 trim(to_str(t % id)))
            call t % read_results_hdf5(tally_group)
            call read_dataset(t % n_realizations, tally_group, &
                 "n_realizations")
            call close_group(tally_group)
          end associate
        end do TALLY_RESULTS

        call close_group(tallies_group)
      end if
    end if

    ! Read source if in eigenvalue mode
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file
        call file_close(file_id)

        ! Write message
        call write_message("Loading source file " // trim(path_source_point) &
             // "...", 5)

        ! Open source file
        file_id = file_open(path_source_point, 'r', parallel=.true.)
      end if

      ! Read source
      call read_source_bank(file_id)

    end if

    ! Close file
    call file_close(file_id)

  end subroutine load_state_point

end module state_point
