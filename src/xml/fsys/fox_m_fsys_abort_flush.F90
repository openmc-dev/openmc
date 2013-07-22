module fox_m_fsys_abort_flush

  implicit none

  public :: pxfflush
  public :: pxfabort
  public :: pure_pxfabort

  ! pxf.F90 - assortment of Fortran wrappers to various
  ! unix-y system calls.
  !
  ! Copyright Toby White, <tow21@cam.ac.uk> 2005

  ! The name pxf is intended to be reminiscent of the POSIX
  ! fortran interfaces defined by IEEE 1003.9-1992, although
  ! in fact I don't think that either flush or abort were
  ! covered by said standard.

  ! Of the preprocessor defines used here, only xlC is
  ! automatically defined by the appropriate compiler. All
  ! others must be defined by some other mechanism - I
  ! recommend autoconf.


  ! FLUSH: flushes buffered output on a given unit. Not guaranteed
  ! to do anything at all (particularly under MPI when even FLUSHed 
  ! buffers may not make it to the host cpu after an abort.
  !
  ! IMPLEMENTATION: in F2003, this is a native operation called by the
  ! FLUSH statement.
  ! In almost every compiler, there is a FLUSH intrinsic subroutine 
  ! available which takes one argument, the unit to be flushed.
  ! (some will flush all open units when no argument is given - this 
  !  facility is not used here.
  ! NAG complicates matters by having to USE a module to get at flush.
  ! 
  ! If no FLUSH is available, the subroutine is a no-op.


CONTAINS

  subroutine pxfflush(unit)
#ifdef __NAG__
    use f90_unix_io, only : flush
#endif
    integer, intent(in) :: unit

#if defined(F2003)
    flush(unit)
#elif defined(xlC)
    call flush_(unit)
#elif defined (FC_HAVE_FLUSH)
    call flush(unit)
#else
    continue
#endif

  end subroutine pxfflush

  ! ABORT: terminates program execution in such a way that a backtrace
  ! is produced. (see abort() in the C Standard Library). No arguments.
  !
  ! IMPLEMENTATION: In F2003, the C interoperability features mean that
  ! the abort in stdlib.h is available to be linked against.
  ! In several other compilers an ABORT intrinsic subroutine is available.
  ! Again, NAG complicates matters by needing to USE a module.
  !
  ! In the case where no native ABORT can be found, we emulate one
  ! by dereferencing a null pointer. This has reliably produced coredumps
  ! on every platform I've tried it with. Just in case it doesn't (it need
  ! not even stop execution) a stop is given to force termination.

  subroutine pxfabort()
#ifdef __NAG__
    use f90_unix_proc, only : abort
#endif
#ifdef F2003
    interface
      subroutine abort() bind(c)
      end subroutine abort
    end interface
#define FC_HAVE_ABORT
#endif
#ifndef FC_HAVE_ABORT
    integer, pointer :: i
#endif

    call pxfflush(6)

#ifdef FC_HAVE_ABORT
#ifdef FC_ABORT_UNDERSCORE
    call abort_()
#elif defined(FC_ABORT_ARG)
    call abort("")
#else
    call abort()
#endif
#else
    i=>null()
    Print*,i
#endif
    stop

  end subroutine pxfabort

  ! For where we need a pxfabort that is pure, we have
  ! this below. I am less sure of its working everywhe
  ! also it must be used as a function not a subroutine
  ! (otherwise it would be optimized away as side-effect
  ! free

  pure function pure_pxfabort() result (crash)
    integer :: crash
    integer, pointer :: i
    nullify(i)
    crash = i
  end function pure_pxfabort

end module fox_m_fsys_abort_flush
