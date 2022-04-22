#! /bin/bash

export BLADEDIR=/home/rhaye/BLaDE
source $BLADEDIR/modules
export BLADEEXEC=$BLADEDIR/build/blade
source ../setupenv # python with alf module

python -c "import alf; alf.postprocess($i,$eqS,$S,$N,$skipE,True,engine='blade')"
