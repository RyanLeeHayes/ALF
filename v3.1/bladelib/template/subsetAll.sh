#! /bin/bash

module load slurm
export nnodes=`cat nnodes`
export nreps=`cat nreps`

# DEPEND="--dependency=afterok:"
PID=`sbatch --time=2880 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL $DEPEND ./runset2.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
PID=`sbatch --time=2880 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL $DEPEND ./runset3.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"



export nitt=1

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export ini=111
export i=$ini$p
eval PID=`sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-100%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

DEPEND=$NEWDEPEND

export i=111
export eqS=25
export S=100
export N=5
export skipE=10

PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"
