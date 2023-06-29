#! /bin/bash

source env-slurm

# DEPEND="--dependency=afterok:"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL $DEPEND ./runset2.sh | awk '{print $4}'`

DEPEND="--dependency=afterok:$PID"

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export step=211
export p
export nitt=5
PID=`sbatch --time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL --array=1-1%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","

done

DEPEND=$NEWDEPEND

export i=211
export eqS=1
export S=5
export N=5
export skipE=1

# DEPEND="--dependency=afterok:"
# --mem=48G
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 $SLURMOPTS --export=ALL $DEPEND ./postprocess.sh
