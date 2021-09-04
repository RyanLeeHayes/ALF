#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`

export i=113
export eqS=10
export S=40
export N=5

# DEPEND="--dependency=afterok:590822"
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postvolume.sh
