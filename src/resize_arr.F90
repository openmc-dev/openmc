module resize_arr

  use bank_header,        only: Bank
  use particle_header,    only: Particle, ParticleBuffer

  implicit none

!===============================================================================
! EXTEND_ARRAY resizes a 1D array to a larger size as specified, keeping the
! contents of the original array
!===============================================================================

  interface extend_array
    module procedure extend_array_bank, &
                     extend8_array_bank, &
                     extend_array_particlebuffer, &
                     extend8_array_particlebuffer
  end interface

!===============================================================================
! RESIZE_ARRAY resizes a 1D array to a larger size as specified, discarding the
! contects of the original array
!===============================================================================

  interface resize_array
    module procedure resize_array_bank, &
                     resize8_array_bank, &
                     resize_array_particle, &
                     resize8_array_particle, &
                     resize_array_particlebuffer, &
                     resize8_array_particlebuffer
  end interface

contains

!===============================================================================
! EXTEND_ARRAY_BANK
!===============================================================================

  subroutine extend_array_bank(array, current_size, requested_size, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer, intent(inout)                     :: current_size
    integer, intent(in)                        :: requested_size
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      tmp(1:current_size) = array
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine extend_array_bank

!===============================================================================
! EXTEND8_ARRAY_BANK
!===============================================================================

  subroutine extend8_array_bank(array, current_size, requested_size, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer(8), intent(inout)                  :: current_size
    integer(8), intent(in)                     :: requested_size
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      tmp(1:current_size) = array
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine extend8_array_bank

!===============================================================================
! EXTEND_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine extend_array_particlebuffer(array, current_size, requested_size, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer, intent(inout)                           :: current_size
    integer, intent(in)                              :: requested_size
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      tmp(1:current_size) = array
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine extend_array_particlebuffer

!===============================================================================
! EXTEND8_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine extend8_array_particlebuffer(array, current_size, requested_size, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer(8), intent(inout)                        :: current_size
    integer(8), intent(in)                           :: requested_size
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      tmp(1:current_size) = array
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine extend8_array_particlebuffer

!===============================================================================
! RESIZE_ARRAY_BANK
!===============================================================================

  subroutine resize_array_bank(array, current_size, requested_size, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer, intent(inout)                     :: current_size
    integer, intent(in)                        :: requested_size
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize_array_bank

!===============================================================================
! RESIZE8_ARRAY_BANK
!===============================================================================

  subroutine resize8_array_bank(array, current_size, requested_size, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer(8), intent(inout)                  :: current_size
    integer(8), intent(in)                     :: requested_size
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize8_array_bank

!===============================================================================
! RESIZE_ARRAY_PARTICLE
!===============================================================================

  subroutine resize_array_particle(array, current_size, requested_size, alloc_err)

    type(Particle), allocatable, intent(inout)     :: array(:)
    integer, intent(inout)                         :: current_size
    integer, intent(in)                            :: requested_size
    integer, intent(out)                           :: alloc_err

    type(PArticle), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize_array_particle

!===============================================================================
! RESIZE8_ARRAY_PARTICLE
!===============================================================================

  subroutine resize8_array_particle(array, current_size, requested_size, alloc_err)

    type(Particle), allocatable, intent(inout)     :: array(:)
    integer(8), intent(inout)                      :: current_size
    integer(8), intent(in)                         :: requested_size
    integer, intent(out)                           :: alloc_err

    type(PArticle), allocatable                    :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize8_array_particle

!===============================================================================
! RESIZE_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine resize_array_particlebuffer(array, current_size, requested_size, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer, intent(inout)                           :: current_size
    integer, intent(in)                              :: requested_size
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize_array_particlebuffer

!===============================================================================
! RESIZE8_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine resize8_array_particlebuffer(array, current_size, requested_size, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer(8), intent(inout)                        :: current_size
    integer(8), intent(in)                           :: requested_size
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (requested_size > current_size) then

      allocate(tmp(requested_size), STAT=alloc_err)
      call move_alloc(FROM=tmp, TO=array)

      current_size = requested_size

    end if

  end subroutine resize8_array_particlebuffer

end module resize_arr
