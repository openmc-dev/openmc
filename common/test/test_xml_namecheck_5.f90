program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkPITarget
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkPITarget('abcd', xds)
  write(*,'(l1)') checkPITarget('1abcd', xds)
  write(*,'(l1)') checkPITarget(':abcd', xds)
  write(*,'(l1)') checkPITarget('#abcd', xds)
  write(*,'(l1)') checkPITarget('e:abcd', xds)
  write(*,'(l1)') checkPITarget('xmle:abcd', xds)
  write(*,'(l1)') checkPITarget('xMle:abcd', xds)

  xds%xml_version = XML1_1

  write(*,'(l1)') checkPITarget('abcd', xds)
  write(*,'(l1)') checkPITarget('1abcd', xds)
  write(*,'(l1)') checkPITarget(':abcd', xds)
  write(*,'(l1)') checkPITarget('#abcd', xds)
  write(*,'(l1)') checkPITarget('e:abcd', xds)
  write(*,'(l1)') checkPITarget('xmle:abcd', xds)
  write(*,'(l1)') checkPITarget('xMle:abcd', xds)

end program test
