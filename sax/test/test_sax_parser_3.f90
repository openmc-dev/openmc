program test_sax_reader

  use m_sax_operate

  use m_handlers

  type(xml_t) :: xt
  integer :: i

  call open_xml_file(xt, "testin.xml", i)

  call parse(xt, &
       startDocument_handler=start_document_handler, &
       endDocument_handler=end_document_handler, &
       startElement_handler=begin_element_handler, &
       endElement_handler=end_element_handler, &
       startPrefixMapping_handler=start_prefix_handler, &
       endPrefixMapping_handler=end_prefix_handler)!, &
       !characters_handler=characters_handler)

  call close_xml_t(xt)

end program test_sax_reader
