#!/bin/bash

export CHARMMDIR=/home/rhaye/CHARMM/chv1/charmm_exe
source $CHARMMDIR/modules
export CHARMMEXEC=$CHARMMDIR/gnu/charmm
source ../setupenv # python with alf module

iend=$(( $SLURM_ARRAY_TASK_ID * $nitt ))
python -c "import alf; alf.runprod($step,'$p',$iend)"
