#! /bin/bash

module load slurm

export i=111
export eqS=25
export S=100
export N=5
export skipE=10

# DEPEND="--dependency=afterok:"
# --mem=48G
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh
