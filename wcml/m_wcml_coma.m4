dnl
include(`foreach.m4')`'dnl
dnl
include(`common.m4')`'dnl
dnl
dnl Below we only use arguments with a type of xsd:string
define(`TOHWM4_bandargs', `(dictRef,convention,title,id,ref,label)')dnl
dnl
define(`TOHWM4_eigenargs', `(dictRef,convention,title,id,type)')dnl
dnl
define(`TOHWM4_bandargslist', `dnl
m4_foreach(`x', TOHWM4_bandargs, `,x')dnl
')dnl
define(`TOHWM4_bandargsdecl',`dnl
m4_foreach(`x',TOHWM4_bandargs,`TOHWM4_dummyargdecl(x)')
')dnl
define(`TOHWM4_bandargsuse',`dnl
m4_foreach(`x',TOHWM4_bandargs,`TOHWM4_dummyarguse(x)')
')dnl
dnl
define(`TOHWM4_eigenargslist', `dnl
m4_foreach(`x', TOHWM4_eigenargs, `,x')dnl
')dnl
define(`TOHWM4_eigenargsdecl',`dnl
m4_foreach(`x',TOHWM4_eigenargs,`TOHWM4_dummyargdecl(x)')
')dnl
define(`TOHWM4_eigenargsuse',`dnl
m4_foreach(`x',TOHWM4_eigenargs,`TOHWM4_dummyarguse(x)')
')dnl
dnl
define(`TOHWM4_coma_real_subs', `dnl`'

  subroutine cmlStartKPoint$1(xf, coords, weight, kptfmt, wtfmt &
TOHWM4_bandargslist)
    type(xmlf_t), intent(inout)              :: xf
    real(kind=$1), dimension(3), intent(in)  :: coords
    real(kind=$1), intent(in), optional      :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
TOHWM4_bandargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "kpoint")
TOHWM4_bandargsuse

    call xml_AddAttribute(xf, "coords", coords, kptfmt)
    if (present(weight)) &
       call xml_AddAttribute(xf, "weight", weight, wtfmt)
#endif

  end subroutine cmlStartKPoint$1

  subroutine cmlAddKPoint$1(xf, coords, weight, kptfmt, wtfmt &
TOHWM4_bandargslist)
    type(xmlf_t), intent(inout)             :: xf
    real(kind=$1), dimension(3), intent(in) :: coords
    real(kind=$1), intent(in), optional     :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
TOHWM4_bandargsdecl

#ifndef DUMMYLIB
    call cmlStartKpoint(xf, coords, weight, kptfmt, wtfmt &
TOHWM4_bandargslist)
    call cmlEndKpoint(xf)
#endif

  end subroutine cmlAddKPoint$1

  subroutine cmlAddEigenValue$1(xf, value, units, fmt &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=$1), intent(in)              :: value
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: fmt

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
TOHWM4_eigenargsuse

    call stmAddValue(xf=xf, value=value, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValue$1

  subroutine cmlAddEigenValueVector$1(xf, eigval, eigvec, units, valfmt, vecfmt &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=$1), intent(in)              :: eigval
    real(kind=$1), intent(in)              :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: valfmt
    character(len=*), intent(in), optional :: vecfmt

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
TOHWM4_eigenargsuse

    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
!FIXME check 2nd dimension of matrix is 3
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVector$1

  subroutine cmlAddEigenValueVectorCmplx$1(xf, eigval, eigvec, units, vecfmt, valfmt &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=$1), intent(in)              :: eigval
    complex(kind=$1), intent(in)           :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: vecfmt
    character(len=*), intent(in), optional :: valfmt

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
TOHWM4_eigenargsuse

    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVectorCmplx$1

  subroutine cmlAddBandList$1(xf, values, fmt, units, spin &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=$1), intent(in)              :: values(:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: spin
    character(len=*), intent(in), optional :: fmt

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "bandList")
    if (present(spin)) then
      if (spin=="up".or.spin=="down") then
        call xml_AddAttribute(xf, "spin", spin)
      else
        !error
      endif
    endif
TOHWM4_eigenargsuse
    call xml_NewElement(xf, "eigen")
    call stmAddValue(xf=xf, value=values, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
    call xml_EndElement(xf, "bandList")
#endif

  end subroutine cmlAddBandList$1

  subroutine cmlAddSymmetry$1(xf, sym_ops, sym_disps, spaceGroup, pointGroup &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)                 :: xf
    real(kind=$1), intent(in)                   :: sym_ops(:,:,:)
    real(kind=$1), intent(in), optional         :: sym_disps(:,:)
    character(len=*), intent(in), optional      :: spaceGroup
    character(len=*), intent(in), optional      :: pointGroup

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB
    integer :: i, n
    real(kind=$1) :: seitzMatrix(4,4)

    call xml_NewElement(xf, "symmetry")
TOHWM4_eigenargsuse
    if (present(spaceGroup)) &
      call xml_AddAttribute(xf, "spaceGroup", spaceGroup)
    if (present(pointGroup)) &
      call xml_AddAttribute(xf, "pointGroup", pointGroup)

    if (present(sym_disps)) then
      if (size(sym_ops, 3)/=size(sym_disps, 2)) then
        ! FIXME error
      endif
    endif
    n = size(sym_ops, 3)

    do i = 1, n
      !Convert the 3x3 rotation and 1x3 translation into a 4x4 Seitz matrix
      if (.not.present(sym_disps)) then
        seitzMatrix = reshape((/sym_ops(:,1,i), 0.0_$1, &
                                sym_ops(:,2,i), 0.0_$1, &
                                sym_ops(:,3,i), 0.0_$1, &
                                0.0_$1, 0.0_$1, 0.0_$1, 1.0_$1/), (/4,4/))
      else
        seitzMatrix = reshape((/sym_ops(:,1,i), sym_disps(1,i), &
                                sym_ops(:,2,i), sym_disps(2,i), &
                                sym_ops(:,3,i), sym_disps(3,i), &
                                0.0_$1, 0.0_$1, 0.0_$1, 1.0_$1/), (/4,4/))
      endif
      call xml_NewElement(xf, "transform3")
      call xml_AddCharacters(xf, chars=seitzMatrix) 
      call xml_EndElement(xf, "transform3")
    end do
    call xml_EndElement(xf, "symmetry")
#endif

    end subroutine cmlAddSymmetry$1

')`'`'dnl
dnl
define(`TOHWM4_coma_subs', `dnl

  subroutine cmlStartBand(xf, spin &
TOHWM4_bandargslist)
    type(xmlf_t), intent(inout)            :: xf
    character(len=*), intent(in), optional :: spin
TOHWM4_bandargsdecl

#ifndef DUMMYLIB
    call xml_NewElement(xf, "band")
    if (present(spin)) then
      if (spin=="up".or.spin=="down") then
        call xml_AddAttribute(xf, "spin", spin)
      else
        !error
      endif
    endif
TOHWM4_bandargsuse
#endif

  end subroutine cmlStartBand

  subroutine cmlEndKpoint(xf)
    type(xmlf_t), intent(inout) :: xf
#ifndef DUMMYLIB
    call xml_EndElement(xf, "kpoint")
#endif
  end subroutine cmlEndKpoint

  subroutine cmlEndBand(xf)
    type(xmlf_t), intent(inout) :: xf
#ifndef DUMMYLIB
    call xml_EndElement(xf, "band")
#endif
  end subroutine cmlEndBand

  subroutine cmlAddSymmetryNoOps(xf, spaceGroup, pointGroup &
TOHWM4_eigenargslist)
    type(xmlf_t), intent(inout)                 :: xf
    character(len=*), intent(in), optional      :: spaceGroup
    character(len=*), intent(in), optional      :: pointGroup

TOHWM4_eigenargsdecl

#ifndef DUMMYLIB

    call xml_NewElement(xf, "symmetry")
TOHWM4_eigenargsuse
    if (present(spaceGroup)) &
      call xml_AddAttribute(xf, "spaceGroup", spaceGroup)
    if (present(pointGroup)) &
      call xml_AddAttribute(xf, "pointGroup", pointGroup)
    call xml_EndElement(xf, "symmetry")
#endif

    end subroutine cmlAddSymmetryNoOps

')`'dnl
dnl
module m_wcml_coma
  ! Implements routines relating to electronic structure

  use fox_m_fsys_realtypes, only: sp, dp
  use FoX_wxml, only: xmlf_t
#ifndef DUMMYLIB
  use FoX_wxml, only: xml_NewElement, xml_AddAttribute
  use FoX_wxml, only: xml_EndElement
  use m_wcml_stml, only: stmAddValue

! Fix for pgi, requires this explicitly:
  use m_wxml_overloads
#endif

  implicit none
  private

  public :: cmlStartKpoint
  public :: cmlEndKpoint
  public :: cmlAddKpoint

  public :: cmlStartBand
  public :: cmlEndBand
  public :: cmlAddBandList
  public :: cmlAddEigenValue
  public :: cmlAddEigenValueVector

  public :: cmlAddSymmetry

  interface cmlAddEigenValue
    module procedure cmlAddEigenValueSP
    module procedure cmlAddEigenValueDP
  end interface

  interface cmlAddEigenValueVector
    module procedure cmlAddEigenValueVectorSP
    module procedure cmlAddEigenValueVectorDP
    module procedure cmlAddEigenValueVectorCmplxSP
    module procedure cmlAddEigenValueVectorCmplxDP
  end interface

  interface cmlStartKpoint
    module procedure cmlStartKpointSP
    module procedure cmlStartKpointDP
  end interface

  interface cmlAddKpoint
    module procedure cmlAddKpointSP
    module procedure cmlAddKpointDP
  end interface

  interface cmlAddBandList
    module procedure cmlAddBandListSP
    module procedure cmlAddBandListDP
  end interface

  interface cmlAddSymmetry
    module procedure cmlAddSymmetrySP
    module procedure cmlAddSymmetryDP
    module procedure cmlAddSymmetryNoOps
  end interface

contains

TOHWM4_coma_real_subs(`sp')`'dnl

TOHWM4_coma_real_subs(`dp')`'dnl

TOHWM4_coma_subs


end module m_wcml_coma
