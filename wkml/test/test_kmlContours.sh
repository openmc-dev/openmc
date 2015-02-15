#!/bin/sh -e

rm -f ${0%.sh}_*_real.*
for t in ${0%.sh}_*.f90
do
  TEST=${t%.f90}
  cat m_contours_test_data_sp.f90 $TEST.f90 > "$TEST"_real.f90
  if [ -f $TEST.xml ]; then ln -s $TEST.xml "$TEST"_real.xml; fi
  if [ -f $TEST.out ]; then ln -s $TEST.out "$TEST"_real.out; fi
  ./test.sh "$TEST"_real
  
done
