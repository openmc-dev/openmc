#!/bin/sh -e

INCFLAGS=`../../FoX-config --fcflags`
export INCFLAGS
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

./test.sh test_UUID
./test.sh test_nullURI
./test.sh test_baseURI
./test.sh test_rfc2396
./test.sh test_rfc2396b
./test.sh test_URI

echo RESULT utils/ Test Results:
echo RESULT utils/ Passed: `wc -l passed.score| cut -f 1 -d 'p'`
echo RESULT utils/ Failed: `wc -l failed.score| cut -f 1 -d 'f'`

echo RESULT utils/ See utils/test/failed.out for details of failed tests.
