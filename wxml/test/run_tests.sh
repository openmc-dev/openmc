#!/bin/sh

INCFLAGS=`../../FoX-config --fcflags`
export INCFLAGS
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

for t in test_xml_*.sh
do
  ./$t
done

echo RESULT wxml/ Test Results:
echo RESULT wxml/ Passed: `grep -c PASSED tests.out`
echo RESULT wxml/ Failed: `grep -c FAILED tests.out`

echo RESULT wxml/ See wxml/test/failed.out for details of failed tests.
