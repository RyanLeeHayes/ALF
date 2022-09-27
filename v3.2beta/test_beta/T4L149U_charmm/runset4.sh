#!/bin/bash

source env-charmm
source ../setupenv # python with alf module

iend=$(( $SLURM_ARRAY_TASK_ID * $nitt ))
python -c "import alf; alf.runprod($step,'$p',$iend)"
