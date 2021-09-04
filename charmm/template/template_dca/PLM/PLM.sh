#! /bin/bash -l

source modules

sleep 30

export OMP_NUM_THREADS=1
# time mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS nvprof -o profile_$SLURM_JOBID.out ./PLM "$@"
time mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS ./PLM "$@"
