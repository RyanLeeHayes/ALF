#! /bin/bash -l

# make clean

OPTLEVEL="-O3"
# OPTLEVEL="-g"
# OPTLEVEL="-g -O3"
# OPTLEVEL="-g -G"

source modules
CC=`which mpicc`
CFLAGS="-qopenmp"
# LDFLAGS="-lm -fopenmp" # gcc syntax?
LDFLAGS="-lm -qopenmp"

CFLAGS="$OPTLEVEL $CFLAGS"

export CC CFLAGS LDLIBS LDFLAGS

make $1 ObservedMoments
