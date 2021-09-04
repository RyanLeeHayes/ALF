#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`
export nitt=10

# DEPEND="--dependency=afterok:"
# NICE="--nice=500"

for p in a b c d e
do

export ini=111
export i=${ini}$p
sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-10%1 $DEPEND $NICE ./runset4.sh

done
