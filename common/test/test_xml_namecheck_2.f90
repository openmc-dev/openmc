program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkQName
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkQName('abcd', xds)
  write(*,'(l1)') checkQName('1abcd', xds)
  write(*,'(l1)') checkQName(':abcd', xds)
  write(*,'(l1)') checkQName('#abcd', xds)
  write(*,'(l1)') checkQName('e:abcd', xds)

  xds%xml_version = XML1_1

  write(*,'(l1)') checkQName('abcd', xds)
  write(*,'(l1)') checkQName('1abcd', xds)
  write(*,'(l1)') checkQName(':abcd', xds)
  write(*,'(l1)') checkQName('#abcd', xds)
  write(*,'(l1)') checkQName('e:abcd', xds)

end program test
