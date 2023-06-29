#! /bin/bash

module load anaconda/3.5.3.0

i=33
NF=5
FREQ=1
ANADIR=`pwd`/../analysis$i
# DEBUG=

cp $ANADIR/b_prev.dat b_prev.dat
cp $ANADIR/c_prev.dat c_prev.dat
cp $ANADIR/x_prev.dat x_prev.dat
cp $ANADIR/s_prev.dat s_prev.dat
cp $ANADIR/../n* ../

export DATDIR=`pwd`/data
mkdir $DATDIR

cd Filter

./Filter.sub $NF $FREQ $ANADIR/data $DATDIR

cd ../Moment

./Moment.sub $NF $DATDIR

cd ../PLM

./PLM.sub

if [[ -v DEBUG ]]
then
  cd ../GetZ
  ./GenerateMoments.sub
fi

cd ../

./GetVariancePLM.py $NF $DATDIR
./GetModel.py $NF $DATDIR

# awk '{print $5,$7}' Result.txt > Result.dat
