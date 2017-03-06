module URR_cross_sections

  use URR_constants, only: ZERO
  implicit none
  private
  public :: CrossSections,&
            xs_samples_tmp

!> Elastic, capture, fission, competitive, and total cross sections at a single
!! energy
  type CrossSections

    real(8) :: n
    real(8) :: g
    real(8) :: f
    real(8) :: x
    real(8) :: t

  contains

    procedure :: accum_resonance => accum_resonance
    procedure :: flush => flush_cross_sections

  end type CrossSections

  type(CrossSections), allocatable :: xs_samples_tmp(:,:)

contains


!> Accumulate contribution from an additional resonance to all partials
  subroutine accum_resonance(this, xs_contribution)

    class(CrossSections), intent(inout) :: this ! all partial xs values
    type(CrossSections) :: xs_contribution ! resonance contribution to xs

    this % n = this % n + xs_contribution % n
    this % g = this % g + xs_contribution % g
    this % f = this % f + xs_contribution % f
    this % x = this % x + xs_contribution % x
    this % t = this % t + xs_contribution % t

  end subroutine accum_resonance


!> Zero out partial xs values
  subroutine flush_cross_sections(this)

    class(CrossSections), intent(inout) :: this ! all partial xs values

    this % n = ZERO
    this % g = ZERO
    this % f = ZERO
    this % x = ZERO
    this % t = ZERO

  end subroutine flush_cross_sections


end module URR_cross_sections
