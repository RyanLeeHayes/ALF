#! /bin/bash

module load python

for Ndim in `seq 2 10`
do
  echo $Ndim
  ./Entropy1.py $Ndim
  ./Entropy12.py $Ndim
  ./Entropy2.py $Ndim
done

./Entropy11.py 10
