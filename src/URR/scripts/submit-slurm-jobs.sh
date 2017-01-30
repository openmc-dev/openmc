#!/bin/bash

export JOB_SCRIPT=purxs.slurm

for dir in */; do
  cp $JOB_SCRIPT $dir
  cd $dir
  echo "Submitting slurm job $PWD/$JOB_SCRIPT"
  sbatch $JOB_SCRIPT
  cd ..
done
