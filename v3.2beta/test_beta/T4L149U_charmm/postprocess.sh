#! /bin/bash

export CHARMMDIR=/home/rhaye/CHARMM/chv1/charmm_exe
source $CHARMMDIR/modules
export CHARMMEXEC=$CHARMMDIR/gnu/charmm
source ../setupenv # python with alf module

python -c "import alf; alf.postprocess($i,$eqS,$S,$N,$skipE,True)"
