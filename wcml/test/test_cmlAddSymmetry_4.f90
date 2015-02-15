program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  real :: sym_ops(3,3,2)
  real :: sym_disps(3,2)

  sym_ops(:,:,1) = reshape((/1.0, 0.0, 0.0, &
                             0.0, 1.0, 0.0, &
                             0.0, 0.0, 1.0/), (/3,3/))

  sym_ops(:,:,2) = reshape((/1.0, 0.0, 0.0, &
                             0.0, 1.0, 0.0, &
                             0.0, 0.0, 1.0/), (/3,3/))

  sym_disps(:,1) = (/0.5, 0.5, 0.5/)
  sym_disps(:,2) = (/-0.5, -0.5, -0.5/)

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  call cmlAddSymmetry(xf, spaceGroup="P -4 21 m", &
    sym_ops=sym_ops, sym_disps=sym_disps)
  call cmlFinishFile(xf)

end program test
