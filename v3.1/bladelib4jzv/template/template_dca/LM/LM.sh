#! /bin/bash -l

source modules

sleep 30

export OMP_NUM_THREADS=8
mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS ./LM "$@"
