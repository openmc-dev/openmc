#!/bin/sh -e

export INCFLAGS=`../../FoX-config --fcflags --wkml`
make clean
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

for t in test_kml*.sh
do
  ./$t
done

echo RESULT wkml/ Test Results:
echo RESULT wkml/ Passed: `wc -l passed.score| cut -f 1 -d 'p'`
echo RESULT wkml/ Failed: `wc -l failed.score| cut -f 1 -d 'f'`

echo RESULT wkml/ See wkml/test/failed.out for details of failed tests.

