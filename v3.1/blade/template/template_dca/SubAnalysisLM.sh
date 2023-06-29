#! /bin/bash

module load slurm

# DEPEND="--dependency=afterok:590822"
sbatch --time=4320 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p brooks --export=ALL $DEPEND ./RunAnalysisLM.sh
