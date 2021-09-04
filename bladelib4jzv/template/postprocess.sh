#! /bin/bash

module load anaconda/3.5.3.0

# i=64
# eqS=5
# S=20
# N=5

sysname=`cat name`
NREPS=`cat nreps`

DIR=`pwd`
PANADIR=$DIR/analysis$(( $i - 1 ))
ANADIR=$DIR/analysis${i}

if [ ! -d $ANADIR ]; then
  mkdir $ANADIR
fi
if [ ! -d $ANADIR/data ]; then
  mkdir $ANADIR/data
fi

# cp -r analysisPhase4/* $ANADIR/
cp $PANADIR/b_sum.dat $ANADIR/b_prev.dat
cp $PANADIR/c_sum.dat $ANADIR/c_prev.dat
cp $PANADIR/x_sum.dat $ANADIR/x_prev.dat
cp $PANADIR/s_sum.dat $ANADIR/s_prev.dat
cd $ANADIR

../ALF/GetLambdas.py $i $N $eqS $S

../ALF/GetEnergy.py $i $i $skipE
../ALF/RunWham.sh $(( $N * $NREPS )) `cat ../ntersiteprod`
../ALF/GetFreeEnergy5.py `cat ../ntersiteprod`
../ALF/SetVars.py $(( $i + 1 ))

if [ `wc -w ../nsubs | awk '{print $1}'` -le 5 ]; then
  ../ALF/GetVariance.py $N
fi
