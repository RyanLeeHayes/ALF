#! /bin/bash -l

source modules

sleep 30

# mpiexec -x OMP_NUM_THREADS=10 ./GenerateMoments ../PLM/h.dat ../PLM/J.dat

# export OMP_NUM_THREADS=10
# ./GenerateMoments ../PLM/h.dat ../PLM/J.dat

export OMP_NUM_THREADS=8
mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS ./GenerateMoments $1 $2 $3 $4
