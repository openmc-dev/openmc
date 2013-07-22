module FoX_wxml

  use m_wxml_core
  use m_wxml_overloads

  implicit none
  private

  public :: xmlf_t
  public :: xml_OpenFile
  public :: xml_Close
  public :: xml_NewElement
  public :: xml_EndElement
  public :: xml_AddXMLDeclaration
  public :: xml_AddXMLStylesheet
  public :: xml_AddXMLPI
  public :: xml_AddComment
  public :: xml_AddEntityReference

  public :: xml_AddCharacters
  public :: xml_AddNewline
  public :: xml_AddAttribute
  public :: xml_AddPseudoAttribute

  public :: xml_DeclareNamespace
  public :: xml_UndeclareNamespace

  public :: xml_AddDOCTYPE
  public :: xml_AddParameterEntity
  public :: xml_AddInternalEntity
  public :: xml_AddExternalEntity
  public :: xml_AddNotation
  public :: xml_AddElementToDTD
  public :: xml_AddAttlistToDTD
  public :: xml_AddPEreferenceToDTD

  public :: xmlf_GetPretty_print
  public :: xmlf_SetPretty_print

  public :: xmlf_Name
  public :: xmlf_OpenTag

end module FoX_wxml
