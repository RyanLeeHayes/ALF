#! /bin/bash

source env-slurm

# DEPEND="--dependency=afterok:"
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL $DEPEND ./runset2.sh
