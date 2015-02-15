program wxml_example

use FoX_wxml

type(xmlf_t) :: xf

integer :: age = 34
real, dimension(20)  :: x
real, dimension(4,4)  :: y

call xml_OpenFile("simple.xml",xf)
call xml_AddDOCTYPE(xf, "john", "hellodtd")
call xml_AddParameterEntity(xf, 'pe', '<!ENTITY def "what a load of nonsense">')
call xml_AddInternalEntity(xf, "abc", "A B C")
call xml_AddElementToDTD(xf, "br", "EMPTY")
call xml_AddAttlistToDTD(xf, "p", "class NMTOKENS #IMPLIED")
call xml_AddPEreferenceToDTD(xf, "pe")
call xml_AddXMLPI(xf, name="robots")
call xml_AddPseudoAttribute(xf, "index", "yes")
call xml_AddPseudoAttribute(xf, "follow", "no")
call xml_AddNotation(xf, name="GIF", system="http://lzw.org")
call xml_AddComment(xf, "a comment ...")
call xml_AddExternalEntity(xf, "def", "http://blah", public="h", notation="GIF")

call xml_AddXMLStylesheet(xf,href="simple.css",type="text/css",media="braille")
call xml_AddXMLPI(xf,name="ccode", data="{hello_world();}")

call xml_NewElement(xf,"john")
call xml_AddAttribute(xf,"age",age)


call xml_AddAttribute(xf,"with_markup","O'Reilly & Assoc is < OUP but > Wiley")
call xml_NewElement(xf,"peter")
call xml_AddComment(xf, "another comment ...")
call xml_NewElement(xf,"tim")
call xml_AddAttribute(xf,"age",37)
call xml_AddAttribute(xf,"weight",123.45d0,fmt="r3")
call xml_AddAttribute(xf,"cholesterol",167.0d0,fmt="r0")
call xml_AddCharacters(xf,"Ping-pong")
call xml_AddCharacters(xf,"champion")
call xml_EndElement(xf,"tim")

call xml_AddCharacters(xf," in years < 2004")

call xml_AddXMLPI(xf, name="robots2")
call xml_AddPseudoAttribute(xf, "index", "if you're nice")
call xml_AddEntityReference(xf, 'abc')

call xml_AddCharacters(xf, repeat("abcd   ",500))

call xml_NewElement(xf,"data")
call xml_AddAttribute(xf,"units","eV")
call random_number(x)
!call xml_AddArray(xf,x)
call xml_EndElement(xf,"data")
call xml_NewElement(xf,"data")
call xml_AddAttribute(xf,"units","Ryd")

call xml_AddEntityReference(xf, '#x2A9')

call xml_DeclareNamespace(xf, "http://www.w3.org/1999/xhtml", "h")
call xml_DeclareNamespace(xf, "http://www.w3.org/1999/svg", "svg")
call xml_NewElement(xf, "h:html")
call xml_NewElement(xf, "svg:svg")
call xml_EndElement(xf, "svg:svg")
call xml_NewElement(xf, "h:head")
call xml_DeclareNamespace(xf,"http://www.xml-cml.org/schema", "cml")
call xml_AddAttribute(xf, "cml:convention", "eMinerals")
call xml_EndElement(xf, "h:head")

!call xml_AddCharacters(xf,(/1, 2, 3, 4, 16 /))

! xml_Close will take care to close all outstanding elements

call xml_Close(xf)

! Equivalent code:
!
!!call xml_EndElement(xf,"data")
!!call xml_EndElement(xf,"peter")
!!call xml_EndElement(xf,"john")


end program wxml_example
