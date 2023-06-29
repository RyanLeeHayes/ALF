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

export ini=61
export i=$ini$p
eval PID=`sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-5%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

DEPEND=$NEWDEPEND

export i=61
export eqS=1
export S=5
export N=5
export skipE=1

PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"



export nitt=1

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export ini=62
export i=$ini$p
eval PID=`sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-20%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

DEPEND=$NEWDEPEND

export i=62
export eqS=5
export S=20
export N=5
export skipE=1

PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"



export nitt=10

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export ini=63
export i=$ini$p
eval PID=`sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-10%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

DEPEND=$NEWDEPEND

export i=63
export eqS=25
export S=100
export N=5
export skipE=5

PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"



export nitt=10

NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do

export ini=64
export i=$ini$p
eval PID=`sbatch --time=1440 --ntasks=$nreps --tasks-per-node=1 --cpus-per-task=$nnodes -p gpu --gres=gpu:gtx1080:$nnodes --export=ALL --array=1-50%1 $DEPEND $NICE ./runset4.sh | awk '{print $4}'`
NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
COMMA=","
echo "$p-$PID"

done

DEPEND=$NEWDEPEND

export i=64
export eqS=125
export S=500
export N=5
export skipE=25

PID=`sbatch --time=2880 --ntasks=1 --tasks-per-node=1 --cpus-per-task=2 -p gpu --gres=gpu:1 --export=ALL $DEPEND ./postprocess.sh | awk '{print $4}'`
echo $PID

DEPEND="--dependency=afterok:$PID"
# NICE="--nice=500"
