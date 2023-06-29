#! /bin/bash -l

source modules

ncentral=`cat ../../ncentral`

sleep 30

export OMP_NUM_THREADS=1
mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS ./Filter $2 $3/Lambda.$1.$ncentral.dat $4/Filter.$1.dat
