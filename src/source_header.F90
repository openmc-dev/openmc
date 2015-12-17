module source_header

  use distribution_univariate, only: Distribution
  use distribution_multivariate, only: UnitSphereDistribution

  implicit none

!===============================================================================
! EXTSOURCE describes an external source of neutrons for a fixed-source problem
! or for the starting source in a k eigenvalue problem
!===============================================================================

  type ExtSource
    integer :: type_space              ! spacial distribution, e.g. 'box' or 'point'
    real(8), allocatable :: params_space(:)  ! parameters for spatial distribution
    class(UnitSphereDistribution), allocatable :: angle  ! angle distribution
    class(Distribution), allocatable :: energy ! energy distribution
  end type ExtSource

end module source_header
