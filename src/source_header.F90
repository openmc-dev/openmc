module source_header

  use distribution_univariate, only: Distribution
  use distribution_multivariate, only: UnitSphereDistribution, SpatialDistribution

  implicit none

!===============================================================================
! SOURCEDISTRIBUTION describes an external source of particles for a
! fixed-source problem or for the starting source in a k eigenvalue problem
!===============================================================================

  type SourceDistribution
    real(8) :: strength ! source strength
    class(SpatialDistribution),    allocatable :: space  ! spatial distribution
    class(UnitSphereDistribution), allocatable :: angle  ! angle distribution
    class(Distribution),           allocatable :: energy ! energy distribution
  end type SourceDistribution

end module source_header
