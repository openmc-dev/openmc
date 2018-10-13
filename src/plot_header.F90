module plot_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt

  implicit none

!===============================================================================
! ObjectColor holds color information for plotted objects
!===============================================================================

  type, bind(C) :: ObjectColor
    integer(C_INT) :: rgb(3)
  end type ObjectColor

!===============================================================================
! PLOTSLICE holds plot information
!===============================================================================

  type, bind(C) :: ObjectPlot
    integer(C_INT) :: id                              ! Unique ID
    integer(C_INT) :: type                            ! Type
    integer(C_INT) :: color_by                        ! quantity to color regions by
    real(C_DOUBLE) :: origin(3)                       ! xyz center of plot location
    real(C_DOUBLE) :: width(3)                        ! xyz widths of plot
    integer(C_INT) :: basis                           ! direction of plot slice
    integer(C_INT) :: pixels(3)                       ! pixel width/height of plot slice
    integer(C_INT) :: meshlines_width                 ! pixel width of meshlines
    integer(C_INT) :: level                           ! universe depth to plot the cells of
    integer(C_INT) :: index_meshlines_mesh = -1       ! index of  mesh to plot
    type(ObjectColor) :: meshlines_color              ! Color for meshlines
    type(ObjectColor) :: not_found                    ! color for positions where no cell found
    type(ObjectColor) :: colors(MAX_COORD)            ! colors of cells/mats
    character(MAX_WORD_LEN, kind=C_CHAR) :: path_plot ! path for plot file
  end type ObjectPlot

  ! Plot type
  integer, parameter :: PLOT_TYPE_SLICE = 1
  integer, parameter :: PLOT_TYPE_VOXEL = 2

  ! Plot level
  integer, parameter :: PLOT_LEVEL_LOWEST = -1

  ! Plot basis plane
  integer, parameter :: PLOT_BASIS_XY = 1
  integer, parameter :: PLOT_BASIS_XZ = 2
  integer, parameter :: PLOT_BASIS_YZ = 3

  ! Indicate whether color refers to unique cell or unique material
  integer, parameter :: PLOT_COLOR_CELLS = 1
  integer, parameter :: PLOT_COLOR_MATS = 2

  integer(C_INT32_T), bind(C) :: n_plots     ! # of plots

  type(ObjectPlot), allocatable, target :: plots(:)

  ! Dictionary that maps user IDs to indices in 'plots'
  type(DictIntInt) :: plot_dict

contains

!===============================================================================
! FREE_MEMORY_PLOT deallocates global arrays defined in this module
!===============================================================================

!===============================================================================
! RUN_PLOT controls the logic for making one or many plots
!===============================================================================

end module plot_header
