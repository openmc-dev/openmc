module FoX_wkml

  use FoX_wxml
  use m_wkml_color
  use m_wkml_lowlevel
  use m_wkml_styling
  use m_wkml_core
  use m_wkml_contours
  use m_wkml_features
  use m_wkml_coverage
  use m_wkml_chart
 
  implicit none
  private

  public :: xmlf_t

! file management function
  public :: kmlBeginFile
  public :: kmlFinishFile
  public :: kmlOpenDocument
  public :: kmlCloseDocument
  public :: kmlOpenFolder
  public :: kmlCloseFolder

! feature function
  public :: kmlCreatePoints
  public :: kmlCreateLine
  public :: kmlStartRegion
  public :: kmlEndRegion
  public :: kmlAddInnerBoundary

! color handling functions and variables
  public :: color_t  
  public :: kmlGetCustomColor
  public :: kmlSetCustomColor

! style handling functions
  public :: kmlCreatePointStyle
  public :: kmlCreateLineStyle
  public :: kmlCreatePolygonStyle

! coverage functions
  public :: kmlCreateContours
  public :: kmlCreateCells
!  public :: kmlCreateRGBCells

end module FoX_wkml
