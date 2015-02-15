#!/bin/sh -e

INCFLAGS=`../../FoX-config --fcflags`
export INCFLAGS
rm -f passed.score failed.score
rm -f tests.out failed.out
touch passed.score failed.score

# I don't know why make won't do this ...

make m_handlers.o
for t in test_sax*.sh
do
  ./$t
done

echo RESULT sax/ Test Results:
echo RESULT sax/ Passed: `wc -l passed.score| cut -f 1 -d 'p'`
echo RESULT sax/ Failed: `wc -l failed.score| cut -f 1 -d 'f'`

echo RESULT sax/ See sax/test/failed.out for details of failed tests.
