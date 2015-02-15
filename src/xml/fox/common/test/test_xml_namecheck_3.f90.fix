program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkNCName
  implicit none

  write(*,'(l1)') checkNCName('abcd', XML1_0)
  write(*,'(l1)') checkNCName('1abcd', XML1_0)
  write(*,'(l1)') checkNCName(':abcd', XML1_0)
  write(*,'(l1)') checkNCName('#abcd', XML1_0)
  write(*,'(l1)') checkNCName('e:abcd', XML1_0)

  write(*,'(l1)') checkNCName('abcd', XML1_1)
  write(*,'(l1)') checkNCName('1abcd', XML1_1)
  write(*,'(l1)') checkNCName(':abcd', XML1_1)
  write(*,'(l1)') checkNCName('#abcd', XML1_1)
  write(*,'(l1)') checkNCName('e:abcd', XML1_1)

end program test
