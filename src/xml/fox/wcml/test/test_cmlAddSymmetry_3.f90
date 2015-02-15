program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  real :: sym_ops(3,3,2)

  sym_ops(:,:,1) = reshape((/1.0, 0.0, 0.0, &
                             0.0, 1.0, 0.0, &
                             0.0, 0.0, 1.0/), (/3,3/))

  sym_ops(:,:,2) = reshape((/1.0, 0.0, 0.0, &
                             0.0, 1.0, 0.0, &
                             0.0, 0.0, 1.0/), (/3,3/))

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddSymmetry(xf, pointGroup="C2v", sym_ops=sym_ops)
  call cmlFinishFile(xf)

end program test
