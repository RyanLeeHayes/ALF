#! /bin/bash

module load slurm

# DEPEND="--dependency=afterok:"
PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./runset2.sh | awk '{print $4}'`

DEPEND="--dependency=afterok:$PID"

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export step=211
export p
export nitt=10
PID=`sbatch --time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=4 -p gpu --gres=gpu:1 --export=ALL --array=1-10%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","

done

DEPEND=$NEWDEPEND

export i=211
export eqS=25
export S=100
export N=5
export skipE=10

# DEPEND="--dependency=afterok:"
# --mem=48G
sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh
