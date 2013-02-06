module plot_header

  use constants

  implicit none

!===============================================================================
! ObjectColor holds color information for plotted objects
!===============================================================================

  type ObjectColor
    integer :: rgb(3)
  end type ObjectColor

!===============================================================================
! PLOTSLICE holds plot information
!===============================================================================

  type PlotSlice
    integer :: id                    ! Unique ID
    character(MAX_LINE_LEN) :: path_plot ! path for plot file
    integer :: type                  ! Type
    integer :: color_by              ! quantity to color regions by
    real(8) :: origin(3)             ! xyz center of plot location
    real(8) :: width(3)              ! xyz widths of plot
    integer :: basis                 ! direction of plot slice 
    integer :: pixels(2)             ! pixel width/height of plot slice
    type(ObjectColor) :: not_found   ! color for positions where no cell found
    type(ObjectColor), allocatable :: colors(:) ! colors of cells/mats
  end type PlotSlice

  ! Plot type -- note that only slice plots are implemented currently
  integer, parameter :: PLOT_TYPE_SLICE = 1
  integer, parameter :: PLOT_TYPE_POINTS = 2

  ! Plot basis plane
  integer, parameter :: PLOT_BASIS_XY = 1
  integer, parameter :: PLOT_BASIS_XZ = 2
  integer, parameter :: PLOT_BASIS_YZ = 3

  ! Indicate whether color refers to unique cell or unique material
  integer, parameter :: PLOT_COLOR_CELLS = 1
  integer, parameter :: PLOT_COLOR_MATS = 2

end module plot_header
