module plot_header

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header, only: DictIntInt
  use mesh_header, only: RegularMesh

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

  type ObjectPlot
    integer :: id                    ! Unique ID
    character(MAX_LINE_LEN) :: path_plot ! path for plot file
    integer :: type                  ! Type
    integer :: color_by              ! quantity to color regions by
    real(8) :: origin(3)             ! xyz center of plot location
    real(8) :: width(3)              ! xyz widths of plot
    integer :: basis                 ! direction of plot slice
    integer :: pixels(3)             ! pixel width/height of plot slice
    integer :: meshlines_width       ! pixel width of meshlines
    integer :: level                 ! universe depth to plot the cells of
    type(RegularMesh), pointer :: meshlines_mesh => null() ! mesh to plot
    type(ObjectColor) :: meshlines_color ! Color for meshlines
    type(ObjectColor) :: not_found   ! color for positions where no cell found
    type(ObjectColor), allocatable :: colors(:) ! colors of cells/mats
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

  subroutine free_memory_plot()
    n_plots = 0
    if (allocated(plots)) deallocate(plots)
    call plot_dict % clear()
  end subroutine free_memory_plot

end module plot_header
