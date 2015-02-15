#!/bin/sh

INCFLAGS=`../../FoX-config --fcflags`
export INCFLAGS
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

for t in test_dom_*.sh
do
  ./$t
done

echo RESULT dom/ Test Results:
echo RESULT dom/ Passed: `grep -c PASSED tests.out`
echo RESULT dom/ Failed: `grep -c FAILED tests.out`

echo RESULT dom/ See dom/test/failed.out for details of failed tests.
