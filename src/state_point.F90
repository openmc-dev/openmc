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

  use constants
  use endf,               only: reaction_name
  use error,              only: fatal_error, warning, write_message
  use hdf5_interface
  use message_passing
  use mgxs_interface
  use nuclide_header,     only: nuclides
  use settings
  use simulation_header
  use string,             only: to_str
  use tally_header
  use tally_filter_header
  use tally_derivative_header

  implicit none

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
    type(TallyDerivative), pointer :: deriv 

    ! Open tallies group
    tallies_group = open_group(file_id, "tallies")

    ! Write information for derivatives.
    if (n_tally_derivs() > 0) then
      derivs_group = create_group(tallies_group, "derivatives")
      do i = 0, n_tally_derivs() - 1
        deriv => tally_deriv_c(i)
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
        id_array(i) = tallies(i) % obj % id()
      end do
      call write_attribute(tallies_group, "ids", id_array)
      deallocate(id_array)

      ! Write all tally information except results
      TALLY_METADATA: do i = 1, n_tallies

        ! Get pointer to tally
        associate (tally => tallies(i) % obj)
        tally_group = create_group(tallies_group, "tally " // &
             trim(to_str(tally % id())))

        ! Write the name for this tally
        call write_dataset(tally_group, "name", tally % name)

        select case(tally % estimator())
        case (ESTIMATOR_ANALOG)
          call write_dataset(tally_group, "estimator", "analog")
        case (ESTIMATOR_TRACKLENGTH)
          call write_dataset(tally_group, "estimator", "tracklength")
        case (ESTIMATOR_COLLISION)
          call write_dataset(tally_group, "estimator", "collision")
        end select
        call write_dataset(tally_group, "n_realizations", &
             tally % n_realizations)

        call write_dataset(tally_group, "n_filters", tally % n_filters())
        if (tally % n_filters() > 0) then
          ! Write IDs of filters
          allocate(id_array(tally % n_filters()))
          do j = 1, tally % n_filters() 
            id_array(j) = filters(tally % filter(j) + 1) % obj % id
          end do
          call write_dataset(tally_group, "filters", id_array)
          deallocate(id_array)
        end if

        ! Set up nuclide bin array and then write
        allocate(str_array(tally % n_nuclide_bins()))
        NUCLIDE_LOOP: do j = 1, tally % n_nuclide_bins()
          if (tally % nuclide_bins(j) >= 0) then
            if (run_CE) then
              i_xs = index(nuclides(tally % nuclide_bins(j)+1) % name, '.')
              if (i_xs > 0) then
                str_array(j) = nuclides(tally % nuclide_bins(j)+1) % name(1 : i_xs-1)
              else
                str_array(j) = nuclides(tally % nuclide_bins(j)+1) % name
              end if
            else
              call get_name_c(tally % nuclide_bins(j)+1, len(temp_name), &
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
        if (tally % deriv() /= C_NONE) then
          deriv => tally_deriv_c(tally % deriv())
          call write_dataset(tally_group, "derivative", deriv % id)
        end if

        ! Write scores.
        call write_dataset(tally_group, "n_score_bins", tally % n_score_bins())
        allocate(str_array(tally % n_score_bins()))
        do j = 1, tally % n_score_bins()
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
      if (active_tallies_size() > 0) then
        ! Indicate that tallies are on
        call write_attribute(file_id, "tallies_present", 1)

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          associate (tally => tallies(i) % obj)
            ! Write sum and sum_sq for each bin
            tally_group = open_group(tallies_group, "tally " &
                 // to_str(tally % id()))
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

  subroutine load_state_point_f(file_id) bind(C)
    integer(HID_T), value :: file_id

    integer :: i
    integer :: temp
    integer(HID_T) :: tallies_group
    integer(HID_T) :: tally_group

    ! Read global tally data
    call read_dataset(global_tallies, file_id, "global_tallies")

    ! Check if tally results are present
    call read_attribute(temp, file_id, "tallies_present")

    ! Read in sum and sum squared
    if (temp == 1) then
      tallies_group = open_group(file_id, "tallies")

      TALLY_RESULTS: do i = 1, n_tallies
        associate (t => tallies(i) % obj)
          ! Read sum, sum_sq, and N for each bin
          tally_group = open_group(tallies_group, "tally " // &
                trim(to_str(t % id())))
          call t % read_results_hdf5(tally_group)
          call read_dataset(t % n_realizations, tally_group, &
                "n_realizations")
          call close_group(tally_group)
        end associate
      end do TALLY_RESULTS

      call close_group(tallies_group)
    end if

  end subroutine load_state_point_f

end module state_point
