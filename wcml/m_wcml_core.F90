module m_wcml_core

  use FoX_wxml, only: xmlf_t

#ifndef DUMMYLIB
  use m_common_error, only: FoX_error

  use FoX_common, only: FoX_version
  use FoX_utils, only: generate_uuid
  use FoX_wxml, only: xml_OpenFile, xml_Close
  use FoX_wxml, only: xml_NewElement, xml_AddAttribute
  use FoX_wxml, only: xml_EndElement, xml_DeclareNamespace
  use FoX_wxml, only: xmlf_Name, xmlf_OpenTag
  use FoX_wxml, only: xmlf_GetExtendedData, xmlf_SetExtendedData

  use m_wcml_metadata, only: cmlAddMetadata
#endif

  implicit none
  private

  public :: cmlBeginFile
  public :: cmlFinishFile
  public :: cmlAddNamespace
  
  public :: cmlStartCml
  public :: cmlEndCml

contains

  subroutine cmlBeginFile(xf, filename, unit, replace)
    type(xmlf_t), intent(out) :: xf
    character(len=*), intent(in) :: filename
    integer, intent(in) :: unit
    logical, intent(in), optional :: replace

#ifndef DUMMYLIB
    if (unit==-1) then
      call xml_OpenFile(filename, xf, preserve_whitespace=.false., replace=replace)
    else
      call xml_OpenFile(filename, xf, preserve_whitespace=.false., unit=unit, replace=replace)
    endif
#endif

  end subroutine cmlBeginFile

  
  subroutine cmlFinishFile(xf)
    type(xmlf_t), intent(inout) :: xf

#ifndef DUMMYLIB
    call xml_Close(xf)
#endif

  end subroutine cmlFinishFile


  subroutine cmlAddNamespace(xf, prefix, URI)
    type(xmlf_t), intent(inout) :: xf  
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: URI

#ifndef DUMMYLIB
    if (xmlf_OpenTag(xf) /= "") &
      call FoX_error("Cannot do cmlAddNamespace after starting CML output")

    call xml_DeclareNamespace(xf, URI, prefix)
#endif

  end subroutine cmlAddNamespace


  subroutine cmlStartCml(xf, id, title, convention, dictref, fileId, version, compchem)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: dictref
    character(len=*), intent(in), optional :: fileId
    character(len=*), intent(in), optional :: version
    logical, intent(in), optional          :: compchem

#ifndef DUMMYLIB
    call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/schema')
    call xml_DeclareNamespace(xf, 'http://www.w3.org/2001/XMLSchema', 'xsd')
    call xml_DeclareNamespace(xf, 'http://www.uszla.me.uk/fpx', 'fpx')
    call xml_DeclareNamespace(xf, 'http://purl.org/dc/elements/1.1/', 'dc')
    call xml_DeclareNamespace(xf, 'http://www.uszla.me.uk/FoX/units', 'units')
    call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/units/units', 'cmlUnits')
    call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/units/siUnits', 'siUnits')
    call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/units/atomic', 'atomicUnits')
    if (present(compchem)) then
      if (compchem) then
        call xmlf_SetExtendedData(xf, 20)
        call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/convention/', 'convention')
        call xml_DeclareNamespace(xf, 'http://www.xml-cml.org/dictionary/compchem/', 'compchem')
      endif
    endif
! FIXME TOHW we may want other namespaces in here - particularly for units
! once PMR has stabilized that.

    call xml_NewElement(xf, 'cml')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(title)) call xml_AddAttribute(xf, 'title', title)
    if (present(dictref)) call xml_AddAttribute(xf, 'dictRef', dictref)
    if (present(convention)) then
      call xml_AddAttribute(xf, 'convention', convention)
    elseif (xmlf_GetExtendedData(xf).eq.20) then
      call xml_AddAttribute(xf, 'convention', 'convention:compchem')
    else
      call xml_AddAttribute(xf, 'convention', 'CMLComp')
    endif
    if (present(fileId)) then
      call xml_AddAttribute(xf, 'fileId', fileId)
    else
      call xml_AddAttribute(xf, 'fileId', xmlf_Name(xf))
    endif
    if (present(version)) then
      call xml_AddAttribute(xf, 'version', version)
    endif

    if (xmlf_GetExtendedData(xf).eq.20) then
      call cmlAddMetadata(xf, name='compchem:UUID', content=generate_uuid(1))
    else
      call cmlAddMetadata(xf, name='UUID', content=generate_uuid(1))
    endif
#endif

  end subroutine cmlStartCml


  subroutine cmlEndCml(xf)
    type(xmlf_t), intent(inout) :: xf

#ifndef DUMMYLIB
    call cmlAddMetadata(xf, name='dc:contributor', content='FoX-'//FoX_version//' (http://www1.gly.bris.ac.uk/~walker/FoX)')
    call xml_EndElement(xf, 'cml')
#endif

  end subroutine cmlEndCml

end module m_wcml_core
