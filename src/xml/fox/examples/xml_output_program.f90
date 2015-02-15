program example_xml_program

  use example_xml_module

  implicit none

  call initialize_xml_output('XML output by FoX')

  call add_paragraph('This output file was produced by FoX.', italic=.false.)

  call add_paragraph('FoX is an XML output library written in Fortran', italic=.true.)

  call end_xml_output

end program example_xml_program
