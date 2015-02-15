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


  subroutine cmlStartKPointsp(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)              :: xf
    real(kind=sp), dimension(3), intent(in)  :: coords
    real(kind=sp), intent(in), optional      :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label



#ifndef DUMMYLIB
    call xml_NewElement(xf, "kpoint")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(ref)) call xml_addAttribute(xf, "ref", ref)
    if (present(label)) call xml_addAttribute(xf, "label", label)



    call xml_AddAttribute(xf, "coords", coords, kptfmt)
    if (present(weight)) &
       call xml_AddAttribute(xf, "weight", weight, wtfmt)
#endif

  end subroutine cmlStartKPointsp

  subroutine cmlAddKPointsp(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)             :: xf
    real(kind=sp), dimension(3), intent(in) :: coords
    real(kind=sp), intent(in), optional     :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label



#ifndef DUMMYLIB
    call cmlStartKpoint(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    call cmlEndKpoint(xf)
#endif

  end subroutine cmlAddKPointsp

  subroutine cmlAddEigenValuesp(xf, value, units, fmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: value
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: fmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=value, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValuesp

  subroutine cmlAddEigenValueVectorsp(xf, eigval, eigvec, units, valfmt, vecfmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: eigval
    real(kind=sp), intent(in)              :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: valfmt
    character(len=*), intent(in), optional :: vecfmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
!FIXME check 2nd dimension of matrix is 3
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVectorsp

  subroutine cmlAddEigenValueVectorCmplxsp(xf, eigval, eigvec, units, vecfmt, valfmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: eigval
    complex(kind=sp), intent(in)           :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: vecfmt
    character(len=*), intent(in), optional :: valfmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVectorCmplxsp

  subroutine cmlAddBandListsp(xf, values, fmt, units, spin &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=sp), intent(in)              :: values(:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: spin
    character(len=*), intent(in), optional :: fmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "bandList")
    if (present(spin)) then
      if (spin=="up".or.spin=="down") then
        call xml_AddAttribute(xf, "spin", spin)
      else
        !error
      endif
    endif
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)


    call xml_NewElement(xf, "eigen")
    call stmAddValue(xf=xf, value=values, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
    call xml_EndElement(xf, "bandList")
#endif

  end subroutine cmlAddBandListsp

  subroutine cmlAddSymmetrysp(xf, sym_ops, sym_disps, spaceGroup, pointGroup &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)                 :: xf
    real(kind=sp), intent(in)                   :: sym_ops(:,:,:)
    real(kind=sp), intent(in), optional         :: sym_disps(:,:)
    character(len=*), intent(in), optional      :: spaceGroup
    character(len=*), intent(in), optional      :: pointGroup

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    integer :: i, n
    real(kind=sp) :: seitzMatrix(4,4)

    call xml_NewElement(xf, "symmetry")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)


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
        seitzMatrix = reshape((/sym_ops(:,1,i), 0.0_sp, &
                                sym_ops(:,2,i), 0.0_sp, &
                                sym_ops(:,3,i), 0.0_sp, &
                                0.0_sp, 0.0_sp, 0.0_sp, 1.0_sp/), (/4,4/))
      else
        seitzMatrix = reshape((/sym_ops(:,1,i), sym_disps(1,i), &
                                sym_ops(:,2,i), sym_disps(2,i), &
                                sym_ops(:,3,i), sym_disps(3,i), &
                                0.0_sp, 0.0_sp, 0.0_sp, 1.0_sp/), (/4,4/))
      endif
      call xml_NewElement(xf, "transform3")
      call xml_AddCharacters(xf, chars=seitzMatrix) 
      call xml_EndElement(xf, "transform3")
    end do
    call xml_EndElement(xf, "symmetry")
#endif

    end subroutine cmlAddSymmetrysp



  subroutine cmlStartKPointdp(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)              :: xf
    real(kind=dp), dimension(3), intent(in)  :: coords
    real(kind=dp), intent(in), optional      :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label



#ifndef DUMMYLIB
    call xml_NewElement(xf, "kpoint")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(ref)) call xml_addAttribute(xf, "ref", ref)
    if (present(label)) call xml_addAttribute(xf, "label", label)



    call xml_AddAttribute(xf, "coords", coords, kptfmt)
    if (present(weight)) &
       call xml_AddAttribute(xf, "weight", weight, wtfmt)
#endif

  end subroutine cmlStartKPointdp

  subroutine cmlAddKPointdp(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)             :: xf
    real(kind=dp), dimension(3), intent(in) :: coords
    real(kind=dp), intent(in), optional     :: weight
    character(len=*), intent(in), optional   :: kptfmt
    character(len=*), intent(in), optional   :: wtfmt
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label



#ifndef DUMMYLIB
    call cmlStartKpoint(xf, coords, weight, kptfmt, wtfmt &
,dictRef,convention,title,id,ref,label)
    call cmlEndKpoint(xf)
#endif

  end subroutine cmlAddKPointdp

  subroutine cmlAddEigenValuedp(xf, value, units, fmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: value
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: fmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=value, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValuedp

  subroutine cmlAddEigenValueVectordp(xf, eigval, eigvec, units, valfmt, vecfmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: eigval
    real(kind=dp), intent(in)              :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: valfmt
    character(len=*), intent(in), optional :: vecfmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
!FIXME check 2nd dimension of matrix is 3
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVectordp

  subroutine cmlAddEigenValueVectorCmplxdp(xf, eigval, eigvec, units, vecfmt, valfmt &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: eigval
    complex(kind=dp), intent(in)           :: eigvec(:,:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: vecfmt
    character(len=*), intent(in), optional :: valfmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "eigen")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)



    call stmAddValue(xf=xf, value=eigval, fmt=valfmt, units=units)
    call stmAddValue(xf=xf, value=eigvec, fmt=vecfmt, units="units:dimensionless")
    call xml_EndElement(xf, "eigen")
#endif

  end subroutine cmlAddEigenValueVectorCmplxdp

  subroutine cmlAddBandListdp(xf, values, fmt, units, spin &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)            :: xf
    real(kind=dp), intent(in)              :: values(:)
    character(len=*), intent(in)           :: units
    character(len=*), intent(in), optional :: spin
    character(len=*), intent(in), optional :: fmt

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    call xml_NewElement(xf, "bandList")
    if (present(spin)) then
      if (spin=="up".or.spin=="down") then
        call xml_AddAttribute(xf, "spin", spin)
      else
        !error
      endif
    endif
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)


    call xml_NewElement(xf, "eigen")
    call stmAddValue(xf=xf, value=values, fmt=fmt, units=units)
    call xml_EndElement(xf, "eigen")
    call xml_EndElement(xf, "bandList")
#endif

  end subroutine cmlAddBandListdp

  subroutine cmlAddSymmetrydp(xf, sym_ops, sym_disps, spaceGroup, pointGroup &
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)                 :: xf
    real(kind=dp), intent(in)                   :: sym_ops(:,:,:)
    real(kind=dp), intent(in), optional         :: sym_disps(:,:)
    character(len=*), intent(in), optional      :: spaceGroup
    character(len=*), intent(in), optional      :: pointGroup

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB
    integer :: i, n
    real(kind=dp) :: seitzMatrix(4,4)

    call xml_NewElement(xf, "symmetry")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)


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
        seitzMatrix = reshape((/sym_ops(:,1,i), 0.0_dp, &
                                sym_ops(:,2,i), 0.0_dp, &
                                sym_ops(:,3,i), 0.0_dp, &
                                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/4,4/))
      else
        seitzMatrix = reshape((/sym_ops(:,1,i), sym_disps(1,i), &
                                sym_ops(:,2,i), sym_disps(2,i), &
                                sym_ops(:,3,i), sym_disps(3,i), &
                                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/4,4/))
      endif
      call xml_NewElement(xf, "transform3")
      call xml_AddCharacters(xf, chars=seitzMatrix) 
      call xml_EndElement(xf, "transform3")
    end do
    call xml_EndElement(xf, "symmetry")
#endif

    end subroutine cmlAddSymmetrydp



  subroutine cmlStartBand(xf, spin &
,dictRef,convention,title,id,ref,label)
    type(xmlf_t), intent(inout)            :: xf
    character(len=*), intent(in), optional :: spin
    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: ref
    character(len=*), intent(in), optional :: label



#ifndef DUMMYLIB
    call xml_NewElement(xf, "band")
    if (present(spin)) then
      if (spin=="up".or.spin=="down") then
        call xml_AddAttribute(xf, "spin", spin)
      else
        !error
      endif
    endif
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(ref)) call xml_addAttribute(xf, "ref", ref)
    if (present(label)) call xml_addAttribute(xf, "label", label)


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
,dictRef,convention,title,id,type)
    type(xmlf_t), intent(inout)                 :: xf
    character(len=*), intent(in), optional      :: spaceGroup
    character(len=*), intent(in), optional      :: pointGroup

    character(len=*), intent(in), optional :: dictRef
    character(len=*), intent(in), optional :: convention
    character(len=*), intent(in), optional :: title
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: type



#ifndef DUMMYLIB

    call xml_NewElement(xf, "symmetry")
    if (present(dictRef)) call xml_addAttribute(xf, "dictRef", dictRef)
    if (present(convention)) call xml_addAttribute(xf, "convention", convention)
    if (present(title)) call xml_addAttribute(xf, "title", title)
    if (present(id)) call xml_addAttribute(xf, "id", id)
    if (present(type)) call xml_addAttribute(xf, "type", type)


    if (present(spaceGroup)) &
      call xml_AddAttribute(xf, "spaceGroup", spaceGroup)
    if (present(pointGroup)) &
      call xml_AddAttribute(xf, "pointGroup", pointGroup)
    call xml_EndElement(xf, "symmetry")
#endif

    end subroutine cmlAddSymmetryNoOps




end module m_wcml_coma
