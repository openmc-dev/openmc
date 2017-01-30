#!/bin/bash -l

export EXE=$HOME/dev/shutmc/src/build/bin/shutmc

##SBATCH -N 1
#SBATCH -n 1
#SBATCH -p pdebug
#SBATCH -t 12:00:00

srun -n1 $EXE $PWD --purxs
