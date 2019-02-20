module tally

  use, intrinsic :: ISO_C_BINDING

  use bank_header
  use constants
  use material_header
  use message_passing
  use settings
  use simulation_header
  use tally_derivative_header
  use tally_filter
  use tally_header

  implicit none

contains

!===============================================================================
! ACCUMULATE_TALLIES accumulates the sum of the contributions from each history
! within the batch to a new random variable
!===============================================================================

  subroutine accumulate_tallies() bind(C)

    integer :: i
    real(C_DOUBLE) :: k_col ! Copy of batch collision estimate of keff
    real(C_DOUBLE) :: k_abs ! Copy of batch absorption estimate of keff
    real(C_DOUBLE) :: k_tra ! Copy of batch tracklength estimate of keff
    real(C_DOUBLE) :: val

#ifdef OPENMC_MPI
    interface
      subroutine reduce_tally_results() bind(C)
      end subroutine
    end interface

    ! Combine tally results onto master process
    if (reduce_tallies) call reduce_tally_results()
#endif

    ! Increase number of realizations (only used for global tallies)
    if (reduce_tallies) then
      n_realizations = n_realizations + 1
    else
      n_realizations = n_realizations + n_procs
    end if

    ! Accumulate on master only unless run is not reduced then do it on all
    if (master .or. (.not. reduce_tallies)) then
      if (run_mode == MODE_EIGENVALUE) then
        if (current_batch > n_inactive) then
          ! Accumulate products of different estimators of k
          k_col = global_tallies(RESULT_VALUE, K_COLLISION) / total_weight
          k_abs = global_tallies(RESULT_VALUE, K_ABSORPTION) / total_weight
          k_tra = global_tallies(RESULT_VALUE, K_TRACKLENGTH) / total_weight
          k_col_abs = k_col_abs + k_col * k_abs
          k_col_tra = k_col_tra + k_col * k_tra
          k_abs_tra = k_abs_tra + k_abs * k_tra
        end if
      end if

      ! Accumulate results for global tallies
      do i = 1, size(global_tallies, 2)
        val = global_tallies(RESULT_VALUE, i)/total_weight
        global_tallies(RESULT_VALUE, i) = ZERO

        global_tallies(RESULT_SUM, i) = global_tallies(RESULT_SUM, i) + val
        global_tallies(RESULT_SUM_SQ, i) = &
             global_tallies(RESULT_SUM_SQ, i) + val*val
      end do
    end if

    ! Accumulate results for each tally
    do i = 1, active_tallies_size()
      call tallies(active_tallies_data(i)+1) % obj % accumulate()
    end do

  end subroutine accumulate_tallies

!===============================================================================
! SETUP_ACTIVE_TALLIES
!===============================================================================

  subroutine setup_active_tallies() bind(C)

    integer :: i
    integer(C_INT) :: err
    logical(C_BOOL) :: active

    interface
      subroutine setup_active_tallies_c() bind(C)
      end subroutine
    end interface

    call setup_active_tallies_c()

    do i = 1, n_tallies
      associate (t => tallies(i) % obj)
        err = openmc_tally_get_active(i, active)
        if (active) then
          ! Check if tally contains depletion reactions and if so, set flag
          if (t % depletion_rx()) need_depletion_rx = .true.
        end if
      end associate
    end do

  end subroutine setup_active_tallies

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_tally_allocate(index, type) result(err) bind(C)
    ! Set the type of the tally
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    integer(C_INT32_T) :: empty(0)
    character(:), allocatable :: type_

    interface
      function tally_pointer(indx) bind(C) result(ptr)
        import C_INT, C_PTR
        integer(C_INT), value :: indx
        type(C_PTR)           :: ptr
      end function
    end interface

    ! Convert C string to Fortran string
    type_ = to_f_string(type)

    err = 0
    if (index >= 1 .and. index <= n_tallies) then
      if (allocated(tallies(index) % obj)) then
        err = E_ALLOCATE
        call set_errmsg("Tally type has already been set.")
      else
        select case (type_)
        case ('generic')
          allocate(TallyObject :: tallies(index) % obj)
        case default
          err = E_UNASSIGNED
          call set_errmsg("Unknown tally type: " // trim(type_))
        end select

        ! Assign the pointer to the C++ tally
        tallies(index) % obj % ptr = tally_pointer(index - 1)

        ! When a tally is allocated, set it to have 0 filters
        err = openmc_tally_set_filters(index, 0, empty)
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in tallies array is out of bounds.")
    end if
  end function openmc_tally_allocate

end module tally
