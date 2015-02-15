function associated_bug(a) result(b)
  integer , pointer :: a
  integer , dimension(merge(1,2,associated(a))) :: b
  b = 0
end function associated_bug

