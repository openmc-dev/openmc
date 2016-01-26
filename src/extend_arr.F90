module extend_arr

  use bank_header,        only: Bank
  use particle_header,    only: Particle, ParticleBuffer

  implicit none

!===============================================================================
! EXTEND_ARRAY resizes a 1D array to a larger size as specified, keeping the
! contents of the original array if required
!===============================================================================

  interface extend_array
    module procedure extend_array_bank, &
                     extend8_array_bank, &
                     extend_array_particle, &
                     extend8_array_particle, &
                     extend_array_particlebuffer, &
                     extend8_array_particlebuffer
  end interface

contains

!===============================================================================
! EXTEND_ARRAY_BANK
!===============================================================================

  subroutine extend_array_bank(array, new_size, keep_data, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer, intent(in)                        :: new_size
    logical, intent(in)                        :: keep_data
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend_array_bank

!===============================================================================
! EXTEND8_ARRAY_BANK
!===============================================================================

  subroutine extend8_array_bank(array, new_size, keep_data, alloc_err)

    type(Bank), allocatable, intent(inout)     :: array(:)
    integer(8), intent(in)                     :: new_size
    logical, intent(in)                        :: keep_data
    integer, intent(out)                       :: alloc_err

    type(Bank), allocatable                    :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend8_array_bank

!===============================================================================
! EXTEND_ARRAY_PARTICLE
!===============================================================================

  subroutine extend_array_particle(array, new_size, keep_data, alloc_err)

    type(Particle), allocatable, intent(inout)     :: array(:)
    integer, intent(in)                            :: new_size
    logical, intent(in)                            :: keep_data
    integer, intent(out)                           :: alloc_err

    type(PArticle), allocatable                    :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend_array_particle

!===============================================================================
! EXTEND8_ARRAY_PARTICLE
!===============================================================================

  subroutine extend8_array_particle(array, new_size, keep_data, alloc_err)

    type(Particle), allocatable, intent(inout)     :: array(:)
    integer(8), intent(in)                         :: new_size
    logical, intent(in)                            :: keep_data
    integer, intent(out)                           :: alloc_err

    type(PArticle), allocatable                    :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend8_array_particle

!===============================================================================
! EXTEND_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine extend_array_particlebuffer(array, new_size, keep_data, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer, intent(in)                              :: new_size
    logical, intent(in)                              :: keep_data
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend_array_particlebuffer

!===============================================================================
! EXTEND8_ARRAY_PARTICLEBUFFER
!===============================================================================

  subroutine extend8_array_particlebuffer(array, new_size, keep_data, alloc_err)

    type(ParticleBuffer), allocatable, intent(inout) :: array(:)
    integer(8), intent(in)                           :: new_size
    logical, intent(in)                              :: keep_data
    integer, intent(out)                             :: alloc_err

    type(ParticleBuffer), allocatable                :: tmp(:)

    if (new_size > size(array)) then

      allocate(tmp(new_size), STAT=alloc_err)
      if (keep_data) tmp(1:size(array)) = array
      call move_alloc(FROM=tmp, TO=array)

    end if

  end subroutine extend8_array_particlebuffer

end module extend_arr
