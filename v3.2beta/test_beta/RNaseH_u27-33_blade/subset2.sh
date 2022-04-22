#! /bin/bash

module load slurm

# DEPEND="--dependency=afterok:"
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset2.sh
