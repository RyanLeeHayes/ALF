#! /bin/bash

source env-slurm

for p in a b c d e
do

export step=211
export p
export nitt=5
sbatch --time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL --array=1-1%1 $DEPEND $NICE ./runset4.sh

done
