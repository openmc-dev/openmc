module URR_resonance

  implicit none
  private
  public :: BWResonanceListVec, BWResonanceVec, BWResonanceVecVec,&
            RMResonanceVec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCE is a Breit-Wigner resonance object
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type :: BWResonance

    real(8) :: E_lam ! resonance energy
    real(8) :: AJ    ! total angular momentum, J
    real(8) :: GN    ! neutron width
    real(8) :: GG    ! gamma width
    real(8) :: GF    ! fission width
    real(8) :: GX    ! competitive width
    real(8) :: GT    ! total width

  end type BWResonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RMRESONANCE is a Reich-Moore resonance object
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type, extends(BWResonance) :: RMResonance

    real(8) :: GFA ! fission width A
    real(8) :: GFB ! fission width B

  end type RMResonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCEVEC is a vector of Breit-Wigner resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type :: BWResonanceVec

    type(BWResonance), allocatable :: res(:)

  ! type-bound procedures
  contains

    ! allocate vector of Breit-Wigner resonances
    procedure :: alloc => alloc_bwresonance_vec

    ! deallocate vector of Breit-Wigner resonances
    procedure :: clear => clear_bwresonance_vec

  end type BWResonanceVec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RMRESONANCEVEC is a vector of Reich-Moore resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type :: RMResonanceVec

    type(RMResonance), allocatable :: res(:)

  ! type-bound procedures
  contains

    ! allocate vector of Reich-Moore resonances
    procedure :: alloc => alloc_rmresonance_vec

    ! deallocate vector of Reich-Moore resonances
    procedure :: clear => clear_rmresonance_vec

  end type RMResonanceVec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCEVECVEC is a vector of length NJS(l) containing a vector of NRS(l,J)
! Breit-Wigner resonances for each J(l)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type :: BWResonanceVecVec

    type(BWResonanceVec), allocatable :: J(:)

  end type BWResonanceVecVec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCEELEM contains one element of a linked list of Breit-Wigner
! resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type, extends(BWResonance) :: BWResonanceElem

    type(BWResonanceElem), pointer :: next => null() ! next resonance in list

  end type BWResonanceElem

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCELIST is a linked list object containing Breit-Wigner resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type BWResonanceList

    type(BWResonanceElem) :: first ! first resonance
    type(BWResonanceElem), pointer :: res => null() ! current resonance

  ! type-bound procedures
  contains

    ! deallocate linked list of Breit-Wigner resonances
    procedure :: clear => clear_bwresonance_list

  end type BWResonanceList

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BWRESONANCELISTVEC is a vector of length NJS(l) containing a linked list of
! NRS(l,J) Breit-Wigner resonances for each J(l)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type BWResonanceListVec

    type(BWResonanceList), allocatable :: J(:)

  end type BWResonanceListVec

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_BWRESONANCE_VEC allocates a vector of Breit-Wigner resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_bwresonance_vec(this, NRS)

    class(BWResonanceVec), intent(inout) :: this ! resonance vector
    integer :: NRS ! number of resonances for this l-wave

    allocate(this % res(NRS))

  end subroutine alloc_bwresonance_vec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_RMRESONANCE_VEC allocates a vector of Reich-Moore resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_rmresonance_vec(this, NRS)

    class(RMResonanceVec), intent(inout) :: this ! resonance vector
    integer :: NRS ! number of resonances for this l-wave

    allocate(this % res(NRS))

  end subroutine alloc_rmresonance_vec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CLEAR_BWRESONANCE_VEC deallocates a vector of Breit-Wigner resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine clear_bwresonance_vec(this)

    class(BWResonanceVec), intent(inout) :: this ! resonance vector

    deallocate(this % res)

  end subroutine clear_bwresonance_vec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CLEAR_RMRESONANCE_VEC deallocates a vector of Reich-Moore resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine clear_rmresonance_vec(this)

    class(RMResonanceVec), intent(inout) :: this ! resonance vector

    deallocate(this % res)

  end subroutine clear_rmresonance_vec

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CLEAR_BWRESONANCE_LIST deallocates a linked list of Breit-Wigner resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine clear_bwresonance_list(this)

    class(BWResonanceList), target, intent(inout) :: this ! resonance element
    type(BWResonanceElem), pointer :: next => null() ! next resonance element

    ! start at beginning of list
    this % res => this % first

    ! while the current resonance element is associated
    do while(associated(this % res))

      ! stash the next resonance element
      next => this % res % next

      ! if the current resonance element is the first, don't need to deallocate
      if (.not.associated(this % res, target=this % first))&
           deallocate(this % res)

      ! point to the next resonance element
      this % res => next

    end do

    ! finish it off
    nullify(next)

  end subroutine clear_bwresonance_list

end module URR_resonance
