program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkEncName
  implicit none

  write(*,'(l1)') checkEncName('abcd')
  write(*,'(l1)') checkEncName('1abcd')
  write(*,'(l1)') checkEncName(':abcd')
  write(*,'(l1)') checkEncName('#abcd')
  write(*,'(l1)') checkEncName('e:abcd')
  write(*,'(l1)') checkEncName('e-abcd')

end program test
