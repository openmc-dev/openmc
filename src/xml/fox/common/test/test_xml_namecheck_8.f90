program test

  use m_common_charset, only : XML1_0, XML1_1
  use m_common_struct, only : xml_doc_state
  use m_common_namecheck, only : checkPEDef
  implicit none

  type(xml_doc_state) :: xds

  write(*,'(l1)') checkPEDef('abcd', xds)
  write(*,'(l1)') checkPEDef('1abcd', xds)
  write(*,'(l1)') checkPEDef('abcd%', xds)
  write(*,'(l1)') checkPEDef('1abc&d', xds)
  write(*,'(l1)') checkPEDef(':abcd%;', xds)
  write(*,'(l1)') checkPEDef('#abcd&;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a1;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a1;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a!;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a!;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%&a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&&a;', xds)
  write(*,'(l1)') checkPEDef('e:a%b;cd%a;', xds)
  write(*,'(l1)') checkPEDef('e:a&b;cd&a;', xds)
  write(*,'(l1)') checkPEDef('&e:;a%b;cd%a1;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; &#x33;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; %#x33;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; &#x07;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; %#x07;', xds)

  xds%xml_version = XML1_1

  write(*,'(l1)') checkPEDef('abcd', xds)
  write(*,'(l1)') checkPEDef('1abcd', xds)
  write(*,'(l1)') checkPEDef('abcd%', xds)
  write(*,'(l1)') checkPEDef('1abc&d', xds)
  write(*,'(l1)') checkPEDef(':abcd%;', xds)
  write(*,'(l1)') checkPEDef('#abcd&;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a1;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a1;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%a!;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&a!;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&%a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd%&a;', xds)
  write(*,'(l1)') checkPEDef('e:abcd&&a;', xds)
  write(*,'(l1)') checkPEDef('e:a%b;cd%a;', xds)
  write(*,'(l1)') checkPEDef('e:a&b;cd&a;', xds)
  write(*,'(l1)') checkPEDef('&e:;a%b;cd%a1;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; &#x33;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; %#x33;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; &#x07;', xds)
  write(*,'(l1)') checkPEDef('%e:;a&b;cd&a1; %#x07;', xds)

end program test
