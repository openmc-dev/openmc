module FoX_sax

  use FoX_common
  use m_sax_operate

  implicit none
  private

  public :: open_xml_file
  public :: open_xml_string
  public :: close_xml_t
  public :: parse
  public :: stop_parser
  public :: SAX_OPEN_ERROR 

  public :: xml_t

  public :: dictionary_t
!SAX functions
  public :: getIndex 
  public :: getLength
  public :: getLocalName
  public :: getQName
  public :: getURI
  public :: getValue
  public :: getType
  public :: isSpecified
  public :: isDeclared
  public :: setSpecified
  public :: setDeclared
!For convenience
  public :: hasKey

end module FoX_sax
