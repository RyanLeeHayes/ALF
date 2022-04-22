#! /bin/bash

module load slurm

export i=211
export NF=5
export BS=50
export FREQ=10

export PHASE=1
# DEPEND="--dependency=afterok:"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --export=ALL $DEPEND ./postprocessLM.sh | awk '{print $4}'`

export PHASE=2
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --export=ALL --array=0-$(( $NF - 1 )) $DEPEND ./postprocessLM.sh | awk '{print $4}'`

export PHASE=3
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --export=ALL $DEPEND ./postprocessLM.sh | awk '{print $4}'`

export PHASE=4
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=8 --cpus-per-task=1 -p gpu --export=ALL --array=0-$BS $DEPEND ./postprocessLM.sh | awk '{print $4}'`

export PHASE=5
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --export=ALL $DEPEND ./postprocessLM.sh | awk '{print $4}'`
