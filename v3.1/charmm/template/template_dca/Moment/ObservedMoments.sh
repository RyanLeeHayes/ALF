#! /bin/bash -l

source modules

ncentral=`cat ../../ncentral`

sleep 30

export OMP_NUM_THREADS=1
mpirun -n 1 -bynode --bind-to none -x OMP_NUM_THREADS ./ObservedMoments $2/Filter.$1.dat $2/m1.$1.obs.dat $2/m2.$1.obs.dat
