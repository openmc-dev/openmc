program dom
 use FoX_dom
 implicit none

 integer :: i
 type(Node), pointer :: doc, name
 type(NodeList), pointer :: nameList
 character(200) :: name_text

 doc => parseString('<?xml version="1.0" encoding="ISO-8859-1"?> <name>[_tmp]:=somecommand(data, 0, 1)</name>')

 nameList => getElementsByTagname(doc, "name")

 do i = 0, getLength(nameList) - 1
   name_text = ''
   name => item(nameList,i)

   name_text = getTextContent(name)

   write(*,*) trim(name_text)
 enddo

 call destroy(doc)
end program dom
