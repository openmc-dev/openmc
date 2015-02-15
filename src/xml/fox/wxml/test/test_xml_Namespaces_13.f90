program test

  use FoX_wxml, only : xmlf_t, xml_OpenFile, xml_Close
  use FoX_wxml, only : xml_NewElement, xml_DeclareNamespace
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  call xml_OpenFile(filename, xf)
  call xml_DeclareNamespace(xf, "ns1", "ns1")
  call xml_NewElement(xf, "ns1:svg")
  call xml_DeclareNamespace(xf, "ns2", "ns2")
  call xml_NewElement(xf, "ns2:svg")
  call xml_DeclareNamespace(xf, "ns3", "ns3")
  call xml_NewElement(xf, "ns3:svg")
  call xml_DeclareNamespace(xf, "ns4", "ns4")
  call xml_NewElement(xf, "ns4:svg")
  call xml_DeclareNamespace(xf, "ns5", "ns5")
  call xml_NewElement(xf, "ns5:svg")
  call xml_DeclareNamespace(xf, "ns6", "ns6")
  call xml_NewElement(xf, "ns6:svg")
  call xml_DeclareNamespace(xf, "ns7", "ns7")
  call xml_NewElement(xf, "ns7:svg")
  call xml_DeclareNamespace(xf, "ns8", "ns8")
  call xml_NewElement(xf, "ns8:svg")
  call xml_DeclareNamespace(xf, "ns9", "ns9")
  call xml_NewElement(xf, "ns9:svg")
  call xml_DeclareNamespace(xf, "ns10", "ns10")
  call xml_NewElement(xf, "ns10:svg")
  call xml_DeclareNamespace(xf, "ns11", "ns11")
  call xml_NewElement(xf, "ns11:svg")

  call xml_Close(xf)

end program test
