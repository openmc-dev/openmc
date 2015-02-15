program test_sax_reader

  use m_sax_operate

  use m_handlers

  type(xml_t) :: xt

  call open_xml_string(xt, &
"<?xml version='1.0'?><!--comment--><?abc xyz?><pqr:aaa xmlns:pqr='http://www.example.com' "//&
"att='value' att2='valuejhg'> <!--c--><![CDATA[<>]]> <b/> jhg </pqr:aaa><")

  call parse(xt, &
       startDocument_handler=start_document_handler, &
       endDocument_handler=end_document_handler, &
       startElement_handler=begin_element_handler, &
       endElement_handler=end_element_handler, &
       startPrefixMapping_handler=start_prefix_handler, &
       endPrefixMapping_handler=end_prefix_handler, &
       characters_handler=characters_handler)

  call close_xml_t(xt)

end program test_sax_reader
