module URR_resonance

  use URR_constants,      only: PI,&
                                FOUR,&
                                SLBW,&
                                MLBW,&
                                REICH_MOORE,&
                                MNBW
  use URR_cross_sections, only: CrossSections
  implicit none
  private
  public :: BreitWignerResonance,&
            BreitWignerResonanceListVector1D,&
            BreitWignerResonanceVector1D,&
            BreitWignerResonanceVector2D,&
            ReichMooreResonanceVector1D,&
            Resonance,&
            wigner_level_spacing


!> Data for a resonance contributing to URR cross sections
  type Resonance

    ! sampled unresolved resonance parameters
    real(8) :: E_lam ! sampled resonance energy
    real(8) :: Gam_t ! sampled total width
    real(8) :: Gam_n ! sampled neutron width
    real(8) :: Gam_g ! sampled radiative width
    real(8) :: Gam_f ! sampled fission width
    real(8) :: Gam_x ! sampled competitive width

    ! counter for the number of resonances added for a given spin sequence
    integer :: i_res

    ! partial and total contributions from resonance to xs values at E_n
    type(CrossSections) :: xs_contribution

  end type Resonance


!> Breit-Wigner resonance parameters
  type :: BreitWignerResonance

    real(8) :: E_lam ! resonance energy
    real(8) :: AJ    ! total angular momentum, J
    real(8) :: GN    ! neutron width
    real(8) :: GG    ! gamma width
    real(8) :: GF    ! fission width
    real(8) :: GX    ! competitive width
    real(8) :: GT    ! total width

  end type BreitWignerResonance


!> Reich-Moore resonance parameters
  type, extends(BreitWignerResonance) :: ReichMooreResonance

    real(8) :: GFA ! fission width A
    real(8) :: GFB ! fission width B

  end type ReichMooreResonance


!> Vector of Breit-Wigner resonances
  type :: BreitWignerResonanceVector1D

    type(BreitWignerResonance), allocatable :: res(:)

  contains

    procedure :: alloc   => alloc_bwresonance_vec
    procedure :: dealloc => dealloc_bwresonance_vec

  end type BreitWignerResonanceVector1D


!> Vector of Reich-Moore resonances
  type :: ReichMooreResonanceVector1D

    type(ReichMooreResonance), allocatable :: res(:)

  contains

    procedure :: alloc   => alloc_rmresonance_vec
    procedure :: dealloc => dealloc_rmresonance_vec

  end type ReichMooreResonanceVector1D


!> Vector of length NJS(l) containing a vector of NRS(l,J) Breit-Wigner
!! resonances for each J(l)
  type :: BreitWignerResonanceVector2D

    type(BreitWignerResonanceVector1D), allocatable :: J(:)

  end type BreitWignerResonanceVector2D


!> One element of a linked list of Breit-Wigner resonances
  type, extends(BreitWignerResonance) :: BreitWignerResonanceElement

    type(BreitWignerResonanceElement), pointer :: next => null() ! next resonance in list

  end type BreitWignerResonanceElement


!> Linked list object containing Breit-Wigner resonances
  type BreitWignerResonanceList

    type(BreitWignerResonanceElement) :: first ! first resonance
    type(BreitWignerResonanceElement), pointer :: res => null() ! current resonance

  contains

    procedure :: dealloc => dealloc_bwresonance_list

  end type BreitWignerResonanceList


!> Vector of length NJS(l) containing a linked list of NRS(l,J) Breit-Wigner
!! resonances for each J(l)
  type BreitWignerResonanceListVector1D

    type(BreitWignerResonanceList), allocatable :: J(:)

  end type BreitWignerResonanceListVector1D

contains


!> Allocates a vector of Breit-Wigner resonances
  subroutine alloc_bwresonance_vec(this, NRS)

    class(BreitWignerResonanceVector1D), intent(inout) :: this ! resonance vector
    integer :: NRS ! number of resonances for this l-wave

    allocate(this % res(NRS))

  end subroutine alloc_bwresonance_vec


!> Allocates a vector of Reich-Moore resonances
  subroutine alloc_rmresonance_vec(this, NRS)

    class(ReichMooreResonanceVector1D), intent(inout) :: this ! resonance vector
    integer :: NRS ! number of resonances for this l-wave

    allocate(this % res(NRS))

  end subroutine alloc_rmresonance_vec


!> Deallocates a vector of Breit-Wigner resonances
  subroutine dealloc_bwresonance_vec(this)

    class(BreitWignerResonanceVector1D), intent(inout) :: this ! resonance vector

    deallocate(this % res)

  end subroutine dealloc_bwresonance_vec


!> Deallocates a vector of Reich-Moore resonances
  subroutine dealloc_rmresonance_vec(this)

    class(ReichMooreResonanceVector1D), intent(inout) :: this ! resonance vector

    deallocate(this % res)

  end subroutine dealloc_rmresonance_vec


!> Deallocates a linked list of Breit-Wigner resonances
  subroutine dealloc_bwresonance_list(this)

    class(BreitWignerResonanceList), target, intent(inout) :: this ! resonance element
    type(BreitWignerResonanceElement), pointer :: next => null() ! next resonance element

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

  end subroutine dealloc_bwresonance_list


!> Samples the Wigner distribution for level spacings
  function wigner_level_spacing(D_avg, prn) result(D_samp)

    real(8) :: D_avg  ! mean level spacing
    real(8) :: prn    ! pseudo-random number
    real(8) :: D_samp ! sampled level spacing

    ! sample a level spacing by directly inverting the Wigner distribution CDF
    D_samp = D_avg * sqrt(-FOUR * log(prn) / PI)

  end function wigner_level_spacing


end module URR_resonance
