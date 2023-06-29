#!/bin/bash

source env-blade
source ../../setupenv # python with alf module

iend=$(( $SLURM_ARRAY_TASK_ID * $nitt ))
python -c "import alf; alf.runprod($step,'$p',$iend,engine='blade')"
