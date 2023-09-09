#! /bin/bash

mkdir multisite
time $( dirname -- "${BASH_SOURCE[0]}" )/whamweight $1 $2 $3 $4

rm Lambda/Lambda.dat
for i in `seq $1`
do
  cat Lambda/Lambda$i.dat >> Lambda/Lambda.dat
done

time $( dirname -- "${BASH_SOURCE[0]}" )/lmalf $2 $3 $4 Lambda/Lambda.dat weight.dat OUT.dat $CRITERIA
