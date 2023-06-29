#! /bin/bash -l

# make clean

# OPTLEVEL="-g"
OPTLEVEL="-O3"
# OPTLEVEL="-g -O3"
# OPTLEVEL="-g -G"

# module load pgi/15.1 cuda/7.0
# CC="nvcc"
# CFLAGS="-x cu -arch=sm_20"
# LDFLAGS=""

source modules

CC=`which mpicc`
CFLAGS="-qopenmp"
# LDFLAGS="-lm -fopenmp" # gcc syntax?
LDFLAGS="-lm -qopenmp"
  
CFLAGS="$OPTLEVEL $CFLAGS"

export CC CFLAGS LDLIBS LDFLAGS

make $1 GenerateMoments
