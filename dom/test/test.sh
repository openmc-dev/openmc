#!/bin/sh

# Core files off, because they are huge on a Mac
ulimit -c 0

# NB Note that we ensure all locally-produced files 
# have Unix line endings only by using 'tr', in
# order to compare properly to our canonical versions.

passed=no
if ! make $1.exe; then
  echo $1 >> failed.out
  echo "------------" >> failed.out
  echo Cannot compile $1 >> failed.out
  echo "------------" >> failed.out
else
  ./$1.exe 2>&1 | tr -d '\15' | grep -v 'STOP' > test.out
  if test -f $1.xml
  then
    tr -d '\15' < test.xml | grep -v UUID > test.xml.tmp; mv test.xml.tmp test.xml
    if test -f test.xml
    then
      if diff test.xml $1.xml > /dev/null; then
        passed=yes
      else
        echo $1 >> failed.out
        echo "------------" >> failed.out
        diff -u test.xml $1.xml >> failed.out
        echo "------------" >> failed.out
      fi
    else
      echo $1 >> failed.out
      echo " -----------"
      echo "test.xml not produced"
      echo " -----------"
    fi
  elif test -f $1.out
  then
    # Note that for most of the test.sh files we only check
    # that the DIFFerences are in one direction. Here, we need
    # to check both directions as the "test_input" file reports
    # errors in a non-standard way (by adding lines). 
    # FIXME: Better to correct test_input?
    if diff -B -b test.out $1.out | grep "^[><]" > /dev/null; then
      echo $1 >> failed.out
      echo "------------" >> failed.out
      diff -u test.out $1.out >> failed.out
      echo "------------" >> failed.out
    else
      passed=yes
    fi
  else
    echo $1 >> failed.out
    echo "------------" >> failed.out
    echo No test output found for $1
    echo "------------" >> failed.out
  fi
fi

if [ $passed = yes ]; then
  echo 'PASSED: ' $1 
  echo 'PASSED: ' $1 >> tests.out
  echo '1' >> passed.score
else
  echo 'FAILED: ' $1 
  echo 'FAILED: ' $1 >> tests.out
  echo '1' >> failed.score
fi
