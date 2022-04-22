#! /bin/bash

module load slurm

for p in a b c d e
do

export step=111
export p
export nitt=10
sbatch --time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL --array=1-10%1 $DEPEND $NICE ./runset4.sh

done
