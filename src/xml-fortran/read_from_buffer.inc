! Part of XML-Fortran library:
!
! $Id: read_from_buffer.inc,v 1.2 2006/03/26 19:05:48 arjenmarkus Exp $
!
   character(len=*), intent(in)                  :: buffer
   integer, intent(inout)                        :: ierror

   integer                                       :: n
   integer                                       :: i
   integer                                       :: step
   integer                                       :: ierr
   !
   ! First allocate an array that is surely large enough
   ! Note:
   ! This is not completely failsafe: with list-directed
   ! input you can also use repeat counts (10000*1.0 for
   ! instance).
   !
   allocate( work(len(buffer)/2+1) )

   !
   ! NOTE:
   ! This is not portable!!
   !
   ! read( buffer, *, iostat = ierror ) (work(n), n=1,size(work))
   !
   ! So, use a different strategy: a binary search
   ! First: establish that we have at least one item to read
   ! Second: do the binary search
   !
!   read( buffer, *, iostat = ierr ) work(1)
!   if ( ierr /= 0 ) then
!      n = 0
!   else
      n = 1
      do while ( n <= size(work) )
         n = 2 * n
      enddo
      n    = n / 2
      step = n / 2
!      step = n / 2

      do while ( step > 0 )
         read( buffer, *, iostat = ierr ) (work(i), i = 1,n)
         if ( ierr /= 0 ) then
            ierror = ierr       ! Store the error code for later use
            n = n - step
         else
            n = n + step
         endif
         step = step / 2
      enddo
!   endif

   !
   ! Then allocate an array of the actual size needed
   ! and copy the data
   !
   !
   if ( associated( var ) ) then
      deallocate( var )
   endif
   !
   ! One complication: we may have one too many
   ! (consequence of the binary search)
   !
   read( buffer, *, iostat = ierr ) (work(i), i = 1,n)
   if ( ierr < 0 ) then
      n = n - 1
   endif

   allocate( var(n) )
   var(1:n) = work(1:n)
   deallocate( work )

   if ( ierror .lt. 0 ) then
      ierror = 0
   endif

