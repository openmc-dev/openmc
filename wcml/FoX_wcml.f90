module FoX_wcml

  use FoX_wxml, only: xmlf_t

  use m_wcml_core
  use m_wcml_coma
  use m_wcml_geometry
  use m_wcml_lattice
  use m_wcml_lists
  use m_wcml_metadata
  use m_wcml_molecule
  use m_wcml_parameter
  use m_wcml_property
  use m_wcml_stml
  use m_wcml_inputdec

  implicit none
  private

  public :: xmlf_t

  public :: cmlBeginFile
  public :: cmlFinishFile
  public :: cmlAddNamespace

  public :: cmlStartCml
  public :: cmlEndCml

  public :: cmlStartModule
  public :: cmlEndModule

  public :: cmlStartStep
  public :: cmlEndStep

!  public :: cmlAddLength
!  public :: cmlAddAngle
!  public :: cmlAddTorsion
!TOHW we don't use these ever

  public :: cmlAddCrystal
  public :: cmlAddLattice

  public :: cmlStartPropertyList
  public :: cmlEndPropertyList
  public :: cmlAddProperty

  public :: cmlStartParameterList
  public :: cmlEndParameterList
  public :: cmlAddParameter

  public :: cmlStartMetadataList
  public :: cmlEndMetadataList
  public :: cmlAddMetadata

  public :: cmlStartMolecule
  public :: cmlEndMolecule
  public :: cmlAddAtoms
  public :: cmlAddParticles
  public :: cmlAddMolecule

  public :: cmlStartBand
  public :: cmlEndBand

  public :: cmlAddBandList
  public :: cmlAddEigenValue
  public :: cmlAddEigenValueVector

  public :: cmlStartKPoint
  public :: cmlEndKPoint
  public :: cmlAddKPoint
  public :: cmlStartKPointList
  public :: cmlEndKPointList

  public :: cmlAddSymmetry

  public :: wcmlDumpDec
  public :: wcmlStartDecList, wcmlEndDecList
  public :: wcmlStartDec, wcmlEndDec
  public :: wcmlAddDecLine

end module FoX_wcml
