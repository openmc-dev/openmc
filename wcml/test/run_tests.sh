#!/bin/sh

INCFLAGS=`../../FoX-config --fcflags`
export INCFLAGS
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

for t in test_cml?*.sh
do
  ./$t
done

echo RESULT wcml/ Test Results:
echo RESULT wcml/ Passed: `wc -l passed.score| cut -f 1 -d 'p'`
echo RESULT wcml/ Failed: `wc -l failed.score| cut -f 1 -d 'f'`

echo RESULT wcml/ See wcml/test/failed.out for details of failed tests.
