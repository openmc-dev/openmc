module ppmlib

  implicit none
 
!===============================================================================
! Image holds RGB information for output PPM image
!===============================================================================

  type Image
    integer, dimension(:,:), pointer :: red, green, blue
    integer                          :: width, height
  end type Image

contains
 
  subroutine init_image(img)

    type(Image), intent(out) :: img

    nullify(img % red)
    nullify(img % green)
    nullify(img % blue)

    img % width = 0
    img % height = 0

  end subroutine init_image

  subroutine allocate_image(img, w, h)

    type(Image), intent(inout) :: img
    integer,     intent(in)    :: w, h
 
    allocate(img % red(w, h))
    allocate(img % green(w, h))
    allocate(img % blue(w, h))

    img % width = w
    img % height = h

  end subroutine allocate_image
 
  subroutine deallocate_image(img)

    type(Image) :: img
 
    if ( associated(img % red) )   deallocate(img % red)
    if ( associated(img % green) ) deallocate(img % green)
    if ( associated(img % blue) )  deallocate(img % blue)
  end subroutine deallocate_image

 
  function inside_image(img, x, y) result(inside)

    type(Image), intent(in) :: img
    integer,     intent(in) :: x, y
    logical                 :: inside
 
    inside = .false.

    if ( ( x < img % width)   .and. &
         ( y < img % height ) .and. &
         ( x >= 0 )           .and. &
         ( y >= 0 ) ) inside = .true.

  end function inside_image
 

  function valid_image(img) result(valid)

    type(Image), intent(in) :: img
    logical                 :: valid

    valid = .false.

    if ( img % width  == 0 ) return
    if ( img % height == 0 ) return
    if ( .not. associated(img % red)   .or. &
         .not. associated(img % green) .or. &
         .not. associated(img % blue) ) return

    valid = .true.

  end function valid_image

 
  subroutine set_pixel(img, x, y, r, g, b)

    type(Image), intent(inout) :: img
    integer,     intent(in)    :: x, y
    integer,     intent(in)    :: r, g, b

    if ( inside_image(img, x, y) .and. valid_image(img)) then
      img % red(x+1,y+1)    = mod(abs(r), 256)
      img % green(x+1, y+1) = mod(abs(g), 256)
      img % blue(x+1, y+1)  = mod(abs(b), 256)
    end if
  end subroutine set_pixel

 
end module ppmlib
