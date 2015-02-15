#!/bin/sh

if [ $1 = yes ]; then
  echo 'PASSED: ' $2
  echo 'PASSED: ' $2 >> tests.out
  echo -n '1' >> passed.score
else
  echo 'FAILED: ' $2
  echo 'FAILED: ' $2 >> tests.out
  echo -n '1' >> failed.score
fi
