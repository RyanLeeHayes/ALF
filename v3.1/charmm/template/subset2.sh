#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`

# DEPEND="--dependency=afterok:"
sbatch --time=2880 --ntasks=$(($nnodes * $nreps)) --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset2.sh
