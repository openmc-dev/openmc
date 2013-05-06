module plot_header

  use constants
  use mesh_header,             only: StructuredMesh

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
    type(ObjectColor) :: not_found   ! color for positions where no cell found
    type(ObjectColor), allocatable :: colors(:) ! colors of cells/mats
    type(StructuredMesh) :: pixmesh  ! pixmesh for reaction rate plots
    real(8), allocatable :: fisswgt(:)! Weights for rxn plot fissionable pixels
    real(8), allocatable :: fluxwgt(:)! Weights for rxn plot non-fiss pix
    integer              :: rrtype   ! particle event type to score
  end type ObjectPlot

  ! Plot type
  integer, parameter :: PLOT_TYPE_SLICE = 1
  integer, parameter :: PLOT_TYPE_VOXEL = 2
  integer, parameter :: PLOT_TYPE_RXNRATE = 3

  ! Plot basis plane
  integer, parameter :: PLOT_BASIS_XY = 1
  integer, parameter :: PLOT_BASIS_XZ = 2
  integer, parameter :: PLOT_BASIS_YZ = 3

  ! Indicate whether color refers to unique cell or unique material
  integer, parameter :: PLOT_COLOR_CELLS = 1
  integer, parameter :: PLOT_COLOR_MATS = 2

  ! Reaction rate plot types
  integer, parameter :: PLOT_RXN_FLUX_THERMAL = 1
  integer, parameter :: PLOT_RXN_FLUX_FAST = 2

  contains
  
end module plot_header
