program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkPseudoAttValue
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkPseudoAttValue('abcd', xds)
  write(*,'(l1)') checkPseudoAttValue('1abcd', xds)
  write(*,'(l1)') checkPseudoAttValue('1abc&d', xds)
  write(*,'(l1)') checkPseudoAttValue('#abcd&;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a1;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a!;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&amp;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&&a;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &y;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &#x33;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &#x07;', xds)

  xds%xml_version = XML1_1

  write(*,'(l1)') checkPseudoAttValue('abcd', xds)
  write(*,'(l1)') checkPseudoAttValue('1abcd', xds)
  write(*,'(l1)') checkPseudoAttValue('1abc&d', xds)
  write(*,'(l1)') checkPseudoAttValue('#abcd&;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a1;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&a!;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&amp;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:abcd&&a;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &y;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &#x33;', xds)
  write(*,'(l1)') checkPseudoAttValue('e:a&amp;cd&gt; &#x07;', xds)

end program test
