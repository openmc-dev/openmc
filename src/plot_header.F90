module plot_header

  implicit none

!===============================================================================
! PLOT hold plot information
!===============================================================================

  type Plot
     integer :: id                    ! Unique ID
     integer :: type                  ! Type
     integer :: color                 ! quantity to color regions by
     real(8) :: origin(3)             ! xyz center of plot location
     real(8) :: aspect                ! spacing between rays in raytracer
     real(8) :: width(3)              ! xyz widths of plot
     integer :: basis                 ! direction of plot slice 
     integer :: pixels(2)             ! pixel width/height of plot slice
  end type Plot

  integer :: PLOT_TYPE_SLICE = 1
  integer :: PLOT_TYPE_POINTS = 2

  integer :: PLOT_BASIS_XY = 1
  integer :: PLOT_BASIS_XZ = 2
  integer :: PLOT_BASIS_YZ = 3

  integer :: PLOT_COLOR_CELLS = 1
  integer :: PLOT_COLOR_MATS = 2


end module plot_header
