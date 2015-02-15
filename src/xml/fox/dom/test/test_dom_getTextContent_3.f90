program dom
 use FoX_dom
 implicit none

 integer :: i
 type(Node), pointer :: doc, name
 type(NodeList), pointer :: nameList
 character(200) :: name_text

 doc => parseFile("test_dom_getTextContent_3.xml_in")

 nameList => getElementsByTagname(doc, "name")

 do i = 0, getLength(nameList) - 1
   name_text = ''
   name => item(nameList,i)

   name_text = getTextContent(name)

   write(*,*) trim(name_text)
 enddo

 call destroy(doc)
end program dom
