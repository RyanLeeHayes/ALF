#! /bin/bash

source env-slurm

export i=211
export eqS=1
export S=5
export N=5
export skipE=1

# DEPEND="--dependency=afterok:"
# --mem=48G
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL $DEPEND ./postprocess.sh
