program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkCharacterEntityReference
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkCharacterEntityReference('#x00', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#x01', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#x0A', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#x0b', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#x20', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#xd7ff', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#xd800', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#00', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#01', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#10', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#11', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#32', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#55295', XML1_0)
  write(*,'(l1)') checkCharacterEntityReference('#55296', XML1_0)

  xds%xml_version = XML1_1

  write(*,'(l1)') checkCharacterEntityReference('#x00', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#x01', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#x0A', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#x0b', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#x20', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#xd7ff', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#xd800', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#00', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#01', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#10', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#11', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#32', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#55295', XML1_1)
  write(*,'(l1)') checkCharacterEntityReference('#55296', XML1_1)

end program test
