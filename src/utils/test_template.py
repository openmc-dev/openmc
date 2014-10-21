#!/usr/bin/env python2

import os
import inspect

template = """module test_{name}

  use global,           only: master, message
  use error,            only: warning
  use testing_header,   only: TestSuiteClass, TestClass

  implicit none
  private

  type, extends(TestClass) :: test
    contains
      procedure         :: init     => test_init
      procedure, nopass :: setup    => test_setup
      procedure, nopass :: execute  => test_execute
      procedure, nopass :: check    => test_check
      procedure, nopass :: teardown => test_teardown
  end type test

  type(test), public :: {name}_test

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(self)

    class(test), intent(inout) :: self

    self % name = "test_{name}"
    
    ! Do any initialization specific to unit testing here

  end subroutine test_init

!===============================================================================
! SETUP
!===============================================================================

  subroutine test_setup(suite)

    class(TestSuiteClass), intent(inout) :: suite
    
    ! Do any openmc setup needed before calling the subroutine or function to
    ! check here.  For example, here is where you might hardcode in some
    ! materials, geometry, and particles

    if (master) then
      message = "TEST NOT IMPLEMENTED"
      call warning(force=.true.)
      call suite % skip()
    end if

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    ! Call the subroutines or functions to check here
    
  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    ! Add test checks here.  For example, check that the results are what you
    ! expect given whatever setup was hardcoded into test_setup
    
    ! If success, do:  call suite % pass()
    ! If failure, do:  call suite % fail()
    
  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    ! Add teardown code here.  For example, deallocate all memory that was used.
    
  end subroutine test_teardown

end module test_{name}"""


################################################################################
def parse_options():
  """Process command line arguments"""

  from optparse import OptionParser
  usage = r"""%prog [options] <name>"""
  p = OptionParser(usage=usage)
  p.add_option('-o', '--output', action='store', dest='output',
             default=None, help='Path of test file to create.  By default the '
                                'test is created in src/unittests using the '
                                'provided name')
  parsed = p.parse_args()
  if not parsed[1]:
    p.print_help()
    return parsed
  return parsed

################################################################################
def main(name, o):

  print(name)
  self_ = os.path.realpath(inspect.getfile(inspect.currentframe()))
  dir_ = os.path.join(os.path.dirname(self_), '..')
  dir_ = os.path.join(dir_, 'unittests')
  file_ = os.path.join(dir_, 'test_{0}.F90'.format(name))
  if os.path.exists(file_):
    proceed = raw_input('Overwrite {0}? ([y]/n) '.format(file_))
    if proceed.lower() == 'N'.lower() or proceed.lower() == 'NO'.lower(): return
  print('Creating file:\n{0}'.format(file_))
  with open(file_, 'w') as fh:
    fh.write(template.format(name=name))
  print('Be sure to add the tests to testing.F90')

################################################################################
if __name__ == '__main__':
  (options, args) = parse_options()
  if args:
    main(args[0], options)
