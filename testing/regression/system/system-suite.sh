#!/bin/sh

###########################################
#
# script for running jobs on Kilkenny head
#
###########################################

export TEST_DIR=/home/walshjon/dev/shutmc/testing/regression/system
cd $TEST_DIR

export LD_LIBRARY_PATH=/opt/gcc/4.9.0/lib64:/opt/mpich/3.1.2-gnu/lib

export SHUTXS=/home/walshjon/dev/shutmc/src/build/bin/shutmc

cd $TEST_DIR/case-1
$SHUTXS
cd $TEST_DIR/case-2
$SHUTXS
cd $TEST_DIR/case-3
$SHUTXS
cd $TEST_DIR/case-4
$SHUTXS
cd $TEST_DIR/case-5
$SHUTXS
cd $TEST_DIR/case-6
$SHUTXS
cd $TEST_DIR/case-7
$SHUTXS
cd $TEST_DIR/case-8
$SHUTXS
cd $TEST_DIR/case-9
$SHUTXS
cd $TEST_DIR/case-10
$SHUTXS
cd $TEST_DIR/case-11
$SHUTXS
cd $TEST_DIR/case-12
$SHUTXS
cd $TEST_DIR/case-13
$SHUTXS
cd $TEST_DIR/case-14
$SHUTXS
cd $TEST_DIR/case-15
$SHUTXS
cd $TEST_DIR/case-16
$SHUTXS
cd $TEST_DIR/case-17
$SHUTXS
cd $TEST_DIR/case-18
$SHUTXS
cd $TEST_DIR/case-19
$SHUTXS
cd $TEST_DIR/case-20
$SHUTXS
cd $TEST_DIR/case-21
$SHUTXS
cd $TEST_DIR/case-22
$SHUTXS
cd $TEST_DIR/case-23
$SHUTXS
cd $TEST_DIR/case-24
$SHUTXS
cd $TEST_DIR/case-25
$SHUTXS
