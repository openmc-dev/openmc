program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkSystemID
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkSystemID('abcd')
  write(*,'(l1)') checkSystemID('1abcd /"()+,/=?;!*#@$%')
  write(*,'(l1)') checkSystemID('1abcd /"()+,/=?;!*#@$%'//achar(13))
  write(*,'(l1)') checkSystemID('1abcd /(")+,/=?;!*#@$%&')
  write(*,'(l1)') checkSystemID('1abcd /(")+,/=?;!*#@$%&'//"'")

end program test
