#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`

export i=111
export eqS=25
export S=100
export N=5
export skipE=100

# DEPEND="--dependency=afterok:"
# --mem=48G
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh
