#!/bin/sh

if [ $1 = yes ]; then
  echo 'PASSED: ' $2
  echo 'PASSED: ' $2 >> tests.out
else
  echo 'FAILED: ' $2
  echo 'FAILED: ' $2 >> tests.out
fi
