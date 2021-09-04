#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`
export nitt=1

# DEPEND="--dependency=afterok:"
# NICE="--nice=500"

for p in a b c d e
do

export ini=111
export i=${ini}$p
sbatch --time=1440 --ntasks=$(($nnodes * $nreps)) --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL --array=1-40%1 $DEPEND $NICE ./runset4.sh

done
