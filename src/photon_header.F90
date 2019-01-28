module photon_header

  use, intrinsic :: ISO_C_BINDING

  use constants

  integer :: n_elements       ! Number of photon cross section tables

!===============================================================================
! ELEMENTMICROXS contains cached microscopic photon cross sections for a
! particular element at the current energy
!===============================================================================

  type, bind(C) :: ElementMicroXS
    integer(C_INT) :: index_grid      ! index on element energy grid
    real(C_DOUBLE) :: last_E = ZERO   ! last evaluated energy
    real(C_DOUBLE) :: interp_factor   ! interpolation factor on energy grid
    real(C_DOUBLE) :: total           ! microscropic total photon xs
    real(C_DOUBLE) :: coherent        ! microscopic coherent xs
    real(C_DOUBLE) :: incoherent      ! microscopic incoherent xs
    real(C_DOUBLE) :: photoelectric   ! microscopic photoelectric xs
    real(C_DOUBLE) :: pair_production ! microscopic pair production xs
  end type ElementMicroXS

  type(ElementMicroXS), allocatable, target :: micro_photon_xs(:) ! Cache for each element
!$omp threadprivate(micro_photon_xs)

contains

!===============================================================================
! FREE_MEMORY_PHOTON deallocates/resets global variables in this module
!===============================================================================

  subroutine free_memory_photon()
    interface
      subroutine free_memory_photon_c() bind(C)
      end subroutine
    end interface

    ! Deallocate photon cross section data
    n_elements = 0
  end subroutine free_memory_photon

  function micro_photon_xs_ptr() result(ptr) bind(C)
    type(C_PTR) :: ptr
    if (size(micro_photon_xs) > 0) then
      ptr = C_LOC(micro_photon_xs(1))
    else
      ptr = C_NULL_PTR
    end if
  end function

end module photon_header
