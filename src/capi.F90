module openmc_capi

  use, intrinsic :: ISO_C_BINDING

  use constants, only: K_BOLTZMANN
  use eigenvalue, only: k_sum
  use global

  private
  public :: openmc_cell_set_temperature
  public :: openmc_reset

contains

!===============================================================================
! OPENMC_CELL_SET_TEMPERATURE sets the temperature of a cell
!===============================================================================

  function openmc_cell_set_temperature(id, T) result(err) bind(C)
    integer(C_INT), value, intent(in) :: id  ! id of cell
    real(C_DOUBLE), value, intent(in) :: T
    integer(C_INT) :: err

    integer :: i

    err = -1
    if (allocated(cells)) then
      if (cell_dict % has_key(id)) then
        i = cell_dict % get_key(id)
        associate (c => cells(i))
          if (allocated(c % sqrtkT)) then
            c % sqrtkT(:) = sqrt(K_BOLTZMANN * T)
            err = 0
          end if
        end associate
      end if
    end if
  end function openmc_cell_set_temperature

!===============================================================================
! OPENMC_RESET resets all tallies
!===============================================================================

  subroutine openmc_reset() bind(C)
    integer :: i

    do i = 1, size(tallies)
      tallies(i) % n_realizations = 0
      if (allocated(tallies(i) % results)) then
        tallies(i) % results(:, :, :) = ZERO
      end if
    end do

    n_realizations = 0
    if (allocated(global_tallies)) then
      global_tallies(:, :) = ZERO
    end if
    k_col_abs = ZERO
    k_col_tra = ZERO
    k_abs_tra = ZERO
    k_sum(:) = ZERO

  end subroutine openmc_reset

end module openmc_capi
