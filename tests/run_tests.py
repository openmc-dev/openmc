#!/bin/sh env python

# load in packages
import os
import sys
import nose
from subprocess import call

pwd = os.getcwd()

# run compile test
result = nose.run(argv=['run_tests.py','-v','test_compile'])
if not result:
    sys.exit('Did not pass compile tests.')

# run gfortran tests
print '\n--------------'
print 'gfortran tests'
print '--------------\n'
os.chdir(pwd)
os.rename(pwd + '/../src/openmc-gfortran', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile','test_basic'])
if not result:
    sys.exit('Did not pass gfortran tests')

# run gfortran-dbg tests
print '\n--------------'
print 'gfortran-dbg tests'
print '--------------\n'
os.rename(pwd + '/../src/openmc-gfortran-dbg', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile','test_basic'])
if not result:
    sys.exit('Did not pass gfortran-dbg tests')

# run gfortran-opt tests
print '\n--------------'
print 'gfortran-opt tests'
print '--------------\n'
os.rename(pwd + '/../src/openmc-gfortran-opt', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile','test_basic'])
if not result:
    sys.exit('Did not pass gfortran-opt tests')

# run gfortran-hdf5 tests
print '\n--------------'
print 'gfortran-hdf5 tests'
print '--------------\n'
#os.rename(pwd + '/../src/openmc-gfortran-hdf5', pwd + '/../src/openmc')
#result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile','test_basic'])
#if not result:
#    sys.exit('Did not pass gfortran-hdf5 tests')

# run gfortran-mpi tests
print '\n--------------'
print 'gfortran-mpi tests'
print '--------------\n'
os.rename(pwd + '/../src/openmc-gfortran-mpi', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec','test_basic'])
if not result:
    sys.exit('Did not pass gfortran-mpi tests')

# run gfortran-phdf5 tests
print '\n--------------'
print 'gfortran-phdf5 tests'
print '--------------\n'
#result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
#                        '--mpi-np','3','--mpi-exec',
#                        '/opt/mpich/3.0.4-gnu/bin/mpiexec','test_basic'])
#if not result:
#    sys.exit('Did not pass gfortran-phdf5 tests')

# run gfortran-petsc tests
print '\n--------------'
print 'gfortran-petsc tests'
print '--------------\n'
os.rename(pwd + '/../src/openmc-gfortran-petsc', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec','test_basic'])
if not result:
    sys.exit('Did not pass gfortran-petsc tests')

# run gfortran-phdf5-petsc tests
print '\n--------------'
print 'gfortran-phdf5-petsc tests'
print '--------------\n'
os.rename(pwd + '/../src/openmc-gfortran-phdf5-petsc', pwd + '/../src/openmc')
result = nose.run(argv=['run_tests.py','-v','--exclude','test_compile',
                        '--mpi-np','3','--mpi-exec',
                        '/opt/mpich/3.0.4-gnu/bin/mpiexec','test_basic'])
if not result:
    sys.exit('Did not pass gfortran-phdf5-petsc tests')
