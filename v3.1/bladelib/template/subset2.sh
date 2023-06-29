#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`

# DEPEND="--dependency=afterok:"
sbatch --time=2880 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL $DEPEND ./runset2.sh
