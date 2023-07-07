#! /bin/bash

export SLURMOPTSGPU1="--time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p free-gpu --gres=gpu:1 --export=ALL"
export SLURMOPTSCPU1="--time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p free --export=ALL"
export SLURMOPTSCPU8="--time=1440 --ntasks=1 --tasks-per-node=8 --cpus-per-task=1 -p free --export=ALL"

export i=213
export NF=5
export BS=50
export FREQ=10

export PHASE=1
DEPEND=""
PID=`sbatch --parsable $SLURMOPTSCPU1 $DEPEND ./postprocessLM.sh`

export PHASE=2
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --parsable $SLURMOPTSCPU1 --array=0-$(( $NF - 1 )) $DEPEND ./postprocessLM.sh`

export PHASE=3
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --parsable $SLURMOPTSCPU1 $DEPEND ./postprocessLM.sh`

export PHASE=4
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --parsable $SLURMOPTSCPU8 --array=0-$BS $DEPEND ./postprocessLM.sh`

export PHASE=5
DEPEND="--dependency=afterok:$PID"
PID=`sbatch --parsable $SLURMOPTSCPU1 $DEPEND ./postprocessLM.sh`
