module m_wcml_metadata

  use FoX_wxml, only: xmlf_t
#ifndef DUMMYLIB
  use FoX_wxml, only: xml_NewElement, xml_EndElement
  use FoX_wxml, only: xml_AddAttribute
#endif

  implicit none
  private

  public :: cmlAddMetadata

contains

   subroutine cmlAddMetadata(xf, name, content, convention, dictRef, id, title )
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: content
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: title

#ifndef DUMMYLIB
    call xml_NewElement(xf, "metadata")
    call xml_AddAttribute(xf, "name", name)
    call xml_AddAttribute(xf, name="content", value=content  )
    if (present(dictref)) call xml_AddAttribute(xf, "dictRef", dictref)
    if (present(id)) call xml_AddAttribute(xf, "id", trim(id))
    if (present(title)) call xml_AddAttribute(xf, "title", title)
    if (present(convention)) call xml_AddAttribute(xf, "convention", convention)
    call xml_EndElement(xf, "metadata")
#endif

   end subroutine cmlAddMetadata

end module m_wcml_metadata
