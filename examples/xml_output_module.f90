module example_xml_module

  use FoX_wxml

  implicit none
  private

  type(xmlf_t) :: xmlfile

  save :: xmlfile

  public :: initialize_xml_output
  public :: add_paragraph
  public :: end_xml_output

contains

  subroutine initialize_xml_output(title)
     character(len=*), intent(in) :: title

     call xml_OpenFile('output.xml', xmlfile)

     call xml_DeclareNamespace(xmlfile, 'http://www.w3.org/1999/xhtml')

     call xml_NewElement(xmlfile, 'html')


     call xml_NewElement(xmlfile, 'head')
     call xml_NewElement(xmlfile, 'title')
     call xml_AddCharacters(xmlfile, title)
     call xml_EndElement(xmlfile, 'title')
     call xml_EndElement(xmlfile, 'head')

     call xml_NewElement(xmlfile, 'body')

  end subroutine initialize_xml_output


  subroutine add_paragraph(para, italic)
    character(len=*), intent(in) :: para
    logical :: italic

    call xml_NewElement(xmlfile, 'p')
    if (italic) call xml_AddAttribute(xmlfile, 'style', 'font-style:italic')
    call xml_AddCharacters(xmlfile, para)
    call xml_EndElement(xmlfile, 'p')
  end subroutine add_paragraph


  subroutine end_xml_output

     call xml_EndElement(xmlfile, 'body')

     call xml_Close(xmlfile)
  end subroutine end_xml_output

end module example_xml_module
