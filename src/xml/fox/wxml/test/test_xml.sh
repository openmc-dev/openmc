#!/bin/sh -e

for t in ${0%.sh}_*.f90
do
  TEST=${t%.f90}
  ./test.sh $TEST
done
